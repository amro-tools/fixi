from ase.units import fs
from ase.md import VelocityVerlet
from ase import Atoms
from ase.optimize.fire2 import FIRE2
from pyfixi import check_constraints
import numpy as np
import numpy.typing as npt


def velocity_verlet_step(
    atoms: Atoms, forces: npt.NDArray, dt: float, apply_constraint: bool
):

    p = atoms.get_momenta()
    p += 0.5 * dt * forces
    masses = atoms.get_masses()[:, None]
    r = atoms.get_positions()

    unadjusted_positions = r + dt * p / masses

    # if we have constraints then this will do the first part of the
    # RATTLE algorithm:
    atoms.set_positions(unadjusted_positions, apply_constraint=apply_constraint)

    return unadjusted_positions


def check(atoms):
    fixi_constraint = atoms._get_constraints()[0]

    hij_max, hij_max_v = check_constraints(
        fixi_constraint.fixi_pairs,
        atoms.get_positions(),
        atoms.get_velocities(),
        atoms.cell.cellpar()[:3],
        atoms.get_pbc(),
    )

    assert hij_max < fixi_constraint.tolerance
    assert hij_max_v < fixi_constraint.tolerance


def test_constraint_forces(water_system):
    atoms = water_system
    # atoms.positions += 1e-1 * np.random.uniform(size=atoms.positions.shape)
    atoms.set_momenta(np.zeros(shape=atoms.positions.shape), apply_constraint=False)

    original_positions = atoms.get_positions()

    DT = 5 * fs

    # also save the constraint forces calculated by fixi
    fixi_constraint = atoms._get_constraints()[0]

    forces = atoms.get_forces()

    unadjusted_positions = velocity_verlet_step(
        atoms, forces, dt=DT, apply_constraint=True
    )

    cell_lengths = atoms.cell.cellpar()[:3]
    pbc = atoms.get_pbc()
    virial = fixi_constraint.rattle.get_virial(
        DT, unadjusted_positions, cell_lengths, pbc
    )
    constraint_forces = fixi_constraint.rattle.get_constraint_forces(len(atoms))
    adjusted_positions = atoms.get_positions()

    # Now we want to re-run the velcocity verlet algorithm
    atoms.set_positions(original_positions, apply_constraint=False)
    atoms.set_momenta(np.zeros(shape=atoms.positions.shape), apply_constraint=False)

    forces_potential = atoms.get_forces()
    forces = forces_potential + constraint_forces

    unadjusted_positions = velocity_verlet_step(
        atoms, forces, dt=DT, apply_constraint=False
    )
    positions_no_constraint = atoms.get_positions()

    print(virial)
    print(adjusted_positions)
    print(positions_no_constraint)

    max_diff = np.max(np.abs(adjusted_positions - positions_no_constraint))
    print(f"{max_diff = }")


if __name__ == "__main__":
    import conftest

    test_constraint_forces(conftest.water_system())
