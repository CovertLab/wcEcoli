#!/usr/bin/env python

"""A numpy benchmark,
tweaked from https://gist.github.com/markus-beuckelmann/8bc25531b11158431a5b09a45abd6276
which was based on: http://stackoverflow.com/questions/11443302/compiling-numpy-with-openblas-integration
"""

import pprint
from time import time

import numpy as np

# Take the randomness out of random numbers, for reproducibility.
np.random.seed(0)

size = 1024
A, B = np.random.random((size, size)), np.random.random((size, size))
C, D = np.random.random((size * 128,)), np.random.random((size * 128,))
E = np.random.random((int(size / 2), int(size / 4)))
F = np.random.random((int(size / 2), int(size / 2)))
F = np.dot(F, F.T)
G = np.random.random((int(size / 2), int(size / 2)))


def t1() -> None:  # Matrix multiplication
	N = 20
	t = time()
	for i in range(N):
		np.dot(A, B)
	delta = time() - t
	print('Dotted two %dx%d matrices in %0.2f ms.' % (size, size, 1e3 * delta / N))


def t2() -> None:  # Vector multiplication
	N = 5000
	t = time()
	for i in range(N):
		np.dot(C, D)
	delta = time() - t
	print('Dotted two vectors of length %d in %0.2f µs.' % (size * 128, 1e6 * delta / N))


def t3() -> None:  # Singular Value Decomposition (SVD)
	N = 3
	t = time()
	for i in range(N):
		np.linalg.svd(E, full_matrices = False)
	delta = time() - t
	print("SVD of a %dx%d matrix in %0.2f ms." % (size / 2, size / 4, 1e3 * delta / N))


def t4() -> None:  # Cholesky Decomposition
	N = 3
	t = time()
	for i in range(N):
		np.linalg.cholesky(F)
	delta = time() - t
	print("Cholesky decomposition of a %dx%d matrix in %0.2f ms."
		  % (size / 2, size / 2, 1e3 * delta / N))


def t5() -> None:  # Eigendecomposition
	N = 3
	t = time()
	for i in range(N):
		np.linalg.eig(G)
	delta = time() - t
	print("Eigendecomposition of a %dx%d matrix in %0.2f s." % (size / 2, size / 2, delta / N))


def print_np_dependencies() -> None:
	print('\nNumpy Build Dependencies:')
	if np.lib.NumpyVersion(np.__version__) >= '1.25.0':
		# noinspection PyNoneFunctionAssignment,PyArgumentList,PyTypeChecker
		d: dict = np.show_config(mode='dicts')  # type: ignore
		pprint.pp(d['Build Dependencies'])
	else:
		np.show_config()


if __name__ == '__main__':
	t1()
	t2()
	t3()
	t4()
	t5()

	print_np_dependencies()
