import pyfixi.fixicpp

import numpy as np

n_atoms = 27
n_molecules = 9

mO = 16.0
mH = 1.0

r_OH = 1
r_HH = 1.5

maxiter = 5000
tolerance = 1e-6
cell_lengths = [20.0, 20.0, 20.0]

pbc = [False, False, False]


pairs = []
masses = np.zeros(n_atoms)
masses[::3] = mO
masses[1::3] = mH
masses[2::3] = mH


# Define pairs and masses
for i in range(n_molecules):
    idx_O = 3 * i
    idx_H1 = 3 * i + 1
    idx_H2 = 3 * i + 2

    pair = pyfixi.fixicpp.FixedBondLengthPair()
    pair.i = idx_O
    pair.j = idx_H1
    pair.dij = r_OH
    pairs.append(pair)

    pair = pyfixi.fixicpp.FixedBondLengthPair()
    pair.i = idx_O
    pair.j = idx_H2
    pair.dij = r_OH
    pairs.append(pair)

    pair = pyfixi.fixicpp.FixedBondLengthPair()
    pair.i = idx_H1
    pair.j = idx_H2
    pair.dij = r_HH
    pairs.append(pair)

# Define the positions
positions = []
for i in range(n_molecules):
    rO = 2.0 * np.random.uniform(size=(3))
    rH1 = rO + 0.4 * np.random.uniform(size=(3))
    rH2 = rO + 0.4 * np.random.uniform(size=(3))
    positions.append(rO)
    positions.append(rH1)
    positions.append(rH2)
positions = np.array(positions)

velocities = 0.1 * np.random.uniform(size=(n_atoms, 3))


rattle = pyfixi.fixicpp.Rattle(maxiter, tolerance, pairs)


adjusted_positions = np.array(positions)


rattle.adjust_positions(positions, adjusted_positions, masses, cell_lengths, pbc)
rattle.adjust_velocities(adjusted_positions, velocities, masses, cell_lengths, pbc)

check = pyfixi.check_constraints(
    pairs, adjusted_positions, velocities, cell_lengths, [False, False, False]
)
print(check)

import pickle

with open("pickle_test", "wb") as f:
    pickle.dump(rattle, f)

with open("pickle_test", "rb") as file:
    r = pickle.load(file)
    print(r.maxiter)
    print(r.tolerance)
    print(r.get_pairs()[2])
