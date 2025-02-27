from ase.units import fs
from ase.md.langevin import Langevin
from ase.optimize.fire2 import FIRE2
from pyfixi import check_constraints


def construct_calculator(atoms):
    from ase.calculators.tip4p import TIP4P

    return TIP4P()


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


def test_ase_constraints_nvt(water_system):
    n_iter = 100
    temperature = 300
    logfile = "log.txt"
    friction = 0.5

    atoms = water_system
    atoms.calc = construct_calculator(atoms)

    dt = 1 * fs

    dyn = Langevin(
        atoms,
        timestep=dt,
        temperature_K=temperature,
        friction=friction,
        logfile=logfile,
    )

    dyn.run(steps=n_iter)
    dyn.close()

    check(atoms)


def test_ase_constraints_min(water_system):
    fmax = 5e-1
    logfile = "log_min.txt"

    atoms = water_system
    atoms.calc = construct_calculator(atoms)

    dt = 0.1 * fs

    dyn = FIRE2(
        atoms,
        dt=dt,
        logfile=logfile,
    )

    dyn.run(fmax=fmax)
    dyn.close()

    check(atoms)


if __name__ == "__main__":
    test_ase_constraints_nvt()
    test_ase_constraints_min()
