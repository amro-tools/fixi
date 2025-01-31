# Fixi

`Fixi`, so far, implements bond length constraints using the RATTLE algorithm [1]. `Fixi` is written in `C++20` and has optional support for multi-threading with `OpenMP`. `Fixi` also provides an interface to `ASE` [2,3], so that that it can be used as a drop-in replacement for `ASE`'s `FixedBondLengths` class (see also the example `tests_python/test_water.py`). 

## Installation

To install, clone the repository and run `pip install`:
```bash
pip install .
```

*NOTE:* The `pip install` command will try to compile the C++ source code and build a binary wheel. 
Hence, this step depends on your system being equipped with a working `C++` compiler and the necessary build tools. 
We recommend to obtain them from the shipped conda `environment.yml` file. 
See below:

## Obtaining build dependencies
We use micromamba here, but `conda` or `mamba` will work just as well.
First create the environment, then activate it.

```bash
micromamba create -f environment.yml # only once
micromamba activate fixienv # every time 
```
After you have activated the environment it shoudl be possible to install fixi with `pip install .`

## Advanced building and running the tests

We use `meson` as the build system: 

```bash
meson setup build
meson compile -C build
```

To run the `C++` unit tests:
```bash
meson test -C build
```

To run the tests in `Python` (which are in `pytest`): 

```bash
pytest -v
```

## PYPI

TODO

## Usage

`fixi.FixBondLengths` is a drop-in replacement for `ase.constraints.FixBondLengths`:

```python
from pyfixi.constraints import FixBondLengths
from ase import Atoms
# ...

# Fix the pairwwise bond lengths between the first three atoms
atoms.set_constraint( FixBondLengths( [ [0,1], [0,2], [1,2] ], tolerance=1e-5) )
# ...
```
The `tolerance` and `maxiter` are optional, with default values equal to `1e-13` and `500`, respectively. 

The number of threads is controlled via the `OMP_NUM_THREADS` environment variable, e.g:
```bash
export OMP_NUM_THREADS=2
```

## References
[1] Hans C Andersen, Rattle: A “velocity” version of the shake algorithm for molecular dynamics calculations, Journal of Computational Physics, Volume 52, Issue 1, 1983

[2] ASE (Atomic Simulation Environment), https://gitlab.com/ase/ase

[3] Ask Hjorth Larsen et al 2017 J. Phys.: Condens. Matter 29 2730022