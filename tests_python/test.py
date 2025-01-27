import pyfixi.fixicpp
import numpy as np

pair = pyfixi.fixicpp.FixedBondLengthPair()
pair.i = 0
pair.j = 1
pair.dij = 2.0

pairs = [pair]

maxiter = 500
tolerance = 1e-6
rattle = pyfixi.fixicpp.Rattle(maxiter, tolerance, pairs, [True, True, True])

positions = np.random.uniform(size=(10, 3))
velocities = np.random.uniform(size=(10, 3))

adjusted_positions = np.array(positions)
rattle.adjust_positions(positions, adjusted_positions, [10.0, 10.0, 10.0])
