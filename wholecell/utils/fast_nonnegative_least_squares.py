"""
Faster implementation of nonnegative least squares.
"""

import numpy as np
from scipy.optimize import nnls

def fast_nnls(A, b):
	"""
	Faster implementation of the nonnegative least squares algorithm, which
	returns a nonnegative vector x that minimizes ||Ax - b||_2. This function
	utilizes the property that both matrix A and vector b can be divided into
	matrices and vectors that each form a smaller nonnegative least squares
	problem, which can each be solved independently and the solutions later
	concatenated to yield the full vector x.

	Args:
		A: numpy.ndarray of size (N, M)
		b: numpy.ndarray of size (N, )
	Returns:
		x: numpy.ndarray of size (M, ), the solution to the NNLS problem.
		rnorm: float, the residual (||Ax - b||_2) of the NNLS problem.
	"""
	# Check matrix dimensions
	if A.ndim != 2 or b.ndim != 1:
		raise Exception('A must be a 2-dimensional matrix and b must be a vector.')

	# Divide matrix A into smaller submatrices
	visited_row_indexes = set()
	visited_column_indexes = set()
	submatrix_indexes = []

	def column_DFS(index, all_row_indexes, all_column_indexes):
		"""
		Recursive function to look for columns and rows in matrix A that should
		be grouped into the same NNLS problem.
		"""
		visited_column_indexes.add(index)
		all_column_indexes.append(index)

		for i in np.where(A[:, index] != 0)[0]:
			if i not in visited_row_indexes:
				row_DFS(i, all_row_indexes, all_column_indexes)

	def row_DFS(index, all_row_indexes, all_column_indexes):
		"""
		Recursive function to look for columns and rows in matrix A that should
		be grouped into the same NNLS problem.
		"""
		visited_row_indexes.add(index)
		all_row_indexes.append(index)

		for i in np.where(A[index, :] != 0)[0]:
			if i not in visited_column_indexes:
				column_DFS(i, all_row_indexes, all_column_indexes)

	# Loop through each column of matrix A
	for column_index in range(A.shape[1]):
		# Search for columns and rows that can be grouped into a single NNLS
		# problem as the given column
		if column_index not in visited_column_indexes:
			submatrix_row_indexes = []
			submatrix_column_indexes = []
			column_DFS(column_index, submatrix_row_indexes, submatrix_column_indexes)

			submatrix_indexes.append((
				np.array(submatrix_row_indexes), np.array(submatrix_column_indexes)
				))

	# Initialize x
	x = np.zeros(A.shape[1])

	# Solve NNLS for each subproblem identified above
	for (row_indexes, column_indexes) in submatrix_indexes:
		if len(row_indexes) == 1 and len(column_indexes) == 1:
			x[column_indexes] = b[row_indexes]
		else:
			x_subproblem, _ = nnls(
				A[np.ix_(row_indexes, column_indexes)],
				b[row_indexes]
				)
			x[column_indexes] = x_subproblem

	assert np.all(x >= 0)
	rnorm = np.linalg.norm(A.dot(x) - b)

	return x, rnorm
