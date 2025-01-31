from ase.constraints import FixConstraint
from ase import Atoms
from ase.geometry import find_mic
import pyfixi.fixicpp
import numpy as np


class FixBondLengths(FixConstraint):
    """Constrain bond lengths with RATTLE."""

    def __init__(self, pairs, tolerance=1e-13, bondlengths=None, maxiter=500):
        self.tolerance = tolerance
        self.maxiter = maxiter
        self.pairs = np.asarray(pairs)
        self.bondlengths = bondlengths

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

    def adjust_positions(self, atoms: Atoms, newpositions):
        if self.rattle is None:
            self.bondlengths = self.compute_initial_bond_lengths(atoms, self.pairs)
            self.construct_fixi_objects(self.pairs, self.bondlengths)

        unadjusted_positions = atoms.positions

        cell_lengths = atoms.cell.cellpar()[:3]
        pbc = atoms.get_pbc()
        masses = atoms.get_masses()

        self.hij_max = self.rattle.adjust_positions(
            unadjusted_positions, newpositions, masses, cell_lengths, pbc
        )

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
