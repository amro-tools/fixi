from ase.io import read, write, Trajectory
from ase.units import fs
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.constraints import FixBondLengths
from pyfixi.constraints import FixBondLengths
from pathlib import Path


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

    atoms.set_constraint(FixBondLengths(pairs))


def construct_calculator(atoms):
    from ase.calculators.tip4p import TIP4P

    return TIP4P()


def test_ase_constraints():
    n_iter = 1000
    input_xyz = Path(__file__).parent / "resources/system.xyz"
    output_xyz = Path(__file__).parent / "final.xyz"
    temperature = 300
    logfile = "log.txt"
    friction = 0.5

    # Read the system using ASE
    with open(input_xyz, "r") as f:
        atoms = read(f, format="xyz")

    atoms.set_cell([20.0, 20.0, 20.0])
    atoms.calc = construct_calculator(atoms)

    print(atoms.cell.cellpar())

    dt = 1 * fs

    dyn = Langevin(
        atoms,
        timestep=dt,
        temperature_K=temperature,
        friction=friction,
        logfile=logfile,
    )

    constrain_water(atoms)

    dyn.run(steps=n_iter)
    dyn.close()

    if not output_xyz is None:
        with open(output_xyz, "w") as f:
            write(f, atoms)


test_ase_constraints()
