import pytest
from pathlib import Path
from ase.io import read
from pyfixi.constraints import FixBondLengths


def constrain_water(atoms):
    n_atoms = len(atoms)
    n_molecules = int(n_atoms / 3)

    pairs = []
    for i_molecule in range(n_molecules):
        iO = 3 * i_molecule
        iH1 = iO + 1
        iH2 = iO + 2

        pairs.append([iO, iH1])
        pairs.append([iO, iH2])
        pairs.append([iH1, iH2])

    atoms.set_constraint(FixBondLengths(pairs, tolerance=1e-5))


def construct_calculator(atoms):
    from ase.calculators.tip4p import TIP4P

    return TIP4P()


@pytest.fixture
def water_system():
    input_xyz = Path(__file__).parent / "resources/system.xyz"

    # Read the system using ASE
    with open(input_xyz, "r") as f:
        atoms = read(f, format="xyz")

    atoms.set_cell([20.0, 20.0, 20.0])
    constrain_water(atoms)
    atoms.calc = construct_calculator(atoms)
    return atoms


@pytest.fixture
def lj_dimer():
    import numpy as np
    from ase import Atoms
    from ase.calculators.lj import LennardJones

    positions = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    atoms = Atoms(positions=positions)
    atoms.set_masses([1.0, 1.0])
    atoms.set_pbc(False)
    atoms.calc = LennardJones()
    atoms.set_constraint(
        FixBondLengths(pairs=[(0, 1)], bondlengths=[2.0], tolerance=1e-5)
    )

    return atoms
