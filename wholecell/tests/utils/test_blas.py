#!/usr/bin/env python
"""
Test a Numpy dot product. Repeat with different values of the
OPENBLAS_NUM_THREADS environment variable which should vary the parallelism
**if** Numpy is installed to use OpenBLAS.

Check that the result using 1 thread is as expected to test if the wcEcoli
simulation results might be consistent across platforms.

The number of BLAS threads can change the dot product computation,
presumably by changing evaluation order and thus floating point rounding.

Any of these ways works to run this test. The first one is quiet if the test
passes. The others always show stdout:
	pytest wholecell/tests/utils/test_blas.py
	python -m wholecell.tests.utils.test_blas
	wholecell/tests/utils/test_blas.py
"""

import os
import sys
import time
import unittest
import warnings

import numpy as np

import wholecell.utils.filepath as fp
from wholecell.utils import parallelization


THIS_DIR_PATH = os.path.dirname(os.path.realpath(__file__))


def dot_product() -> tuple[float, float]:
	"""Return (a dot product of 2 data arrays, the time in ns to compute it).
	"""
	diff = np.load(os.path.join(THIS_DIR_PATH, 'diff.npy'))
	mass = np.load(os.path.join(THIS_DIR_PATH, 'mass.npy'))

	elapsed_start = time.monotonic_ns()
	dot = diff.dot(mass)
	elapsed_end = time.monotonic_ns()
	elapsed_nanoseconds = elapsed_end - elapsed_start

	return dot, elapsed_nanoseconds


class Test_blas(unittest.TestCase):
	def test_dot_product(self):
		"""Compare a dot product running with various numbers of OpenBLAS
		threads -- which matters **if** Numpy is installed to use OpenBLAS.
		"""
		products = []
		total_nanoseconds = 0.0
		thread_range = [str(c) for c in range(1, parallelization.cpus() + 1)] + ['']
		print('OPENBLAS')
		print('{:>7}  {:>26} {:>26} {:>11}'.format('THREADS', 'DOT PRODUCT', 'DIFF FROM 1 THREAD', 'NANOSECS'))

		# NOTE: Setting OPENBLAS_NUM_THREADS affects OpenBLAS but it's unlikely
		# to affect macOS Accelerate or other BLAS libraries.
		for num_threads in thread_range:
			env = dict(os.environ, OPENBLAS_NUM_THREADS=num_threads)
			command = 'python -m wholecell.tests.utils.test_blas DOT'.split()
			output = fp.run_cmd(command, env=env).split(',')
			dot = float(output[0])
			nanoseconds = float(output[1])
			total_nanoseconds += nanoseconds

			products.append(dot)
			diff = dot - products[0]
			print('{:>7}  {:26.17g} {:26.17g} {:11,g}'.format(num_threads, dot, diff, nanoseconds))

		print(f'{"":62} {11 * "-"}')
		print(f'{"":62} {round(total_nanoseconds):11,}')
		print()

		# Issue #931: The expected value came from Numpy's copy of openblas on
		# Intel CPUs on Sherlock, Mac, and Linux. But:
		#
		#   * Apple M1 CPU on Mac computed 0.01668380558411259
		#   * Intel CPU on WSL on Windows computed 0.016683805584112667
		#
		# Reproducible simulations require reproducible floating point results,
		# but that might be unachievable across platforms.
		expected = 0.016683805584112754
		if products[0] != expected:
			warnings.warn(f"Didn't reproduce the expected dot product using 1"
						  f" OpenBLAS thread, {products[0]} != {expected}, so"
						  f" simulation results aren't portable.")

		# Check that multi-threaded results are within some tolerance.
		np_products = np.array(products)
		reference = np.full(len(products), products[0])
		np.testing.assert_allclose(np_products, reference, rtol=1.5e-14)


if __name__ == '__main__':
	if sys.argv[1:] == ['DOT']:
		product, nanoseconds = dot_product()
		print(f'{product}, {nanoseconds}')
	else:
		unittest.main()
