from arrow import StochasticSystem
import numpy as np

seed = 807952948
stoich = np.load('stoich.npy')
mol = np.load('complex-counts.npy')
rates = np.load('rates.npy')

system = StochasticSystem(stoich, random_seed=seed)
for i in range(10000):
	if i % 100 == 0:
		print(i)

	result = system.evolve(1, mol, rates)
