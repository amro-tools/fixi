from ase.constraints import FixConstraint
from ase import Atoms
from ase.geometry import find_mic
import pyfixi.fixicpp
import numpy as np
import numpy.typing as npt
from typing import Optional

class FixBondLengths(FixConstraint):
    """Constrain bond lengths with RATTLE."""

    def __init__(
        self, pairs, tolerance=1e-13, bondlengths=None, maxiter=500, atoms=None
    ):
        self.tolerance = tolerance
        self.maxiter = maxiter
        self.pairs = np.asarray(pairs)
        self.bondlengths = bondlengths
        self._scaled_virial : Optional[npt.NDArray] = None

        if bondlengths is None and not atoms is None:
            self.bondlengths = self.compute_initial_bond_lengths(atoms, self.pairs)

        # Try to construct the fixi object
        self.construct_fixi_objects(pairs, bondlengths)

    def get_removed_dof(self, atoms: Atoms):
        return len(self.pairs)

    def compute_initial_bond_lengths(self, atoms: Atoms, pairs):
        dij = []
        for pair in pairs:
            i = pair[0]
            j = pair[1]

            v, vlen = find_mic(
                atoms.positions[i] - atoms.positions[j], atoms.cell, atoms.pbc
            )

            dij.append(vlen)
        return dij

    def construct_fixi_objects(self, pairs, bondlengths):
        if bondlengths is None:
            self.fixi_pairs = None
            self.rattle = None
            return

        fixi_pairs = []
        for pair, dij in zip(pairs, bondlengths):
            fixi_pair = pyfixi.fixicpp.FixedBondLengthPair()
            i_idx = pair[0]
            j_idx = pair[1]
            fixi_pair.i = i_idx
            fixi_pair.j = j_idx
            fixi_pair.dij = dij
            fixi_pairs.append(fixi_pair)

        self.fixi_pairs = fixi_pairs
        self.rattle = pyfixi.fixicpp.Rattle(
            self.maxiter, self.tolerance, self.fixi_pairs
        )

    def adjust_positions(self, atoms: Atoms, newpositions:npt.ArrayLike):
        """Calculates the new adjusted positions using RATTLE. Also calculates the "unscaled virial" 
        from the constraint forces applied. To get the virial, simply multiply by (timestep)^2. 

        Args:
            atoms (Atoms): _description_
            newpositions (_type_): _description_
        """
        if self.rattle is None:
            self.bondlengths = self.compute_initial_bond_lengths(atoms, self.pairs)
            self.construct_fixi_objects(self.pairs, self.bondlengths)

        unadjusted_positions = np.array(newpositions)

        cell_lengths = atoms.cell.cellpar()[:3]
        pbc = atoms.get_pbc()
        masses = atoms.get_masses()

        self.hij_max = self.rattle.adjust_positions(
            unadjusted_positions, newpositions, masses, cell_lengths, pbc
        )

        self._scaled_virial = self.rattle.get_scaled_virial()

    def get_virial(self, dt):
        return self._scaled_virial * dt**2

    def adjust_momenta(self, atoms: Atoms, momenta):
        if self.rattle is None:
            self.bondlengths = self.compute_initial_bond_lengths(atoms, self.pairs)
            self.construct_fixi_objects(self.pairs, self.bondlengths)

        positions = atoms.positions
        masses = atoms.get_masses()

        masses3 = np.vstack((masses, masses, masses)).T

        adjusted_velocities = momenta / masses3

        cell_lengths = atoms.cell.cellpar()[:3]
        pbc = atoms.get_pbc()

        self.hij_max_v = self.rattle.adjust_velocities(
            positions, adjusted_velocities, masses, cell_lengths, pbc
        )

        momenta[:] = adjusted_velocities * masses3

    def adjust_forces(self, atoms, forces):
        """
        Used in energy minimization in ASE
        """
        self.constraint_forces = -forces
        self.adjust_momenta(atoms, forces)
        self.constraint_forces += forces

    def get_indices(self):
        return np.unique(self.pairs.ravel())

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""
        map = np.zeros(len(atoms), int)
        map[ind] = 1
        n = map.sum()
        map[:] = -1
        map[ind] = range(n)
        pairs = map[self.pairs]
        self.pairs = pairs[(pairs != -1).all(1)]

        self.fixi_pairs = self.construct_bond_length_pairs(self.pairs, self.bondlengths)
        self.rattle = pyfixi.fixicpp.rattle(
            self.maxiter, self.tolerance, self.fixi_pairs
        )

        if len(self.pairs) == 0:
            raise IndexError("Constraint not part of slice")

    def todict(self):
        return {
            "name": "FixBondLengths",
            "kwargs": {
                "pairs": self.pairs.tolist(),
                "tolerance": self.tolerance,
                "bondlengths": self.bondlengths,
                "maxiter": self.maxiter,
            },
        }
