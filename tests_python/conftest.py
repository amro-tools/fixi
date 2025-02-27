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


@pytest.fixture
def water_system():
    input_xyz = Path(__file__).parent / "resources/system.xyz"

    # Read the system using ASE
    with open(input_xyz, "r") as f:
        atoms = read(f, format="xyz")

    atoms.set_cell([20.0, 20.0, 20.0])
    constrain_water(atoms)
    return atoms