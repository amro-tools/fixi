import pytest 
from pathlib import Path
from ase.io import read

@pytest.fixture
def water_system():
    input_xyz = Path(__file__).parent / "resources/system.xyz"

    # Read the system using ASE
    with open(input_xyz, "r") as f:
        atoms = read(f, format="xyz")

    atoms.set_cell([20.0, 20.0, 20.0])
    return atoms