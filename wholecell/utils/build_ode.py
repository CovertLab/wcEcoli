"""
Utilities to compile functions, esp. from Sympy-constructed Matrix math.
"""

import numpy as np
from sympy import Matrix
from typing import Any, Callable, Tuple


def build_functions(arguments, expression):
	# type: (str, str) -> Tuple[Callable, Callable]
	"""Build a function from its arguments and source code expression.
	(This USED TO return a second function that had Numba set up to JIT-compile
	the code on demand, but the compiled code time savings don't pay back the
	compilation time in Python 3.9+.)

	This still returns two functions so the JIT compiler could be restored
	someday.

	Args:
		arguments (str): comma-separated lambda argument names
		expression (str): expression to compile

	Returns:
		a function(arguments),
		the same function(arguments) [was a JIT-enabled version]
	"""
	local_dict: dict[str, Any] = {}
	expression = f'def f({arguments}):\n' + expression
	exec(expression, globals(), local_dict)
	f = local_dict['f']

	return f, f


def _matrix_to_array(matrix):
	# type: (Matrix) -> str
	"""Convert a sympy Matrix expression to a function body."""
	rows, cols = matrix.shape
	_ = np  # So the tools won't warn about unused np import.

	function_str = f'\tarr = np.zeros(({rows}, {cols}))\n'
	for i in range(rows):
		for j in range(cols):
			if matrix[i, j] != 0:
				function_str += f'\tarr[{i}, {j}] = {matrix[i, j]}\n'

	function_str += '\treturn arr'
	return function_str


def derivatives(matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized derivatives ODE function(y, t)."""
	return build_functions('y, t',
		_matrix_to_array(matrix) + '.reshape(-1)')

def derivatives_jacobian(jacobian_matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized derivatives ODE Jacobian function(y, t)."""
	return build_functions('y, t', _matrix_to_array(jacobian_matrix))

def rates(matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized rates function(t, y, kf, kr)."""
	return build_functions('t, y, kf, kr',
		_matrix_to_array(matrix) + '.reshape(-1)')

def rates_jacobian(jacobian_matrix):
	# type: (Matrix) -> Tuple[Callable, Callable]
	"""Build an optimized rates Jacobian function(t, y, kf, kr)."""
	return build_functions('t, y, kf, kr',
		_matrix_to_array(jacobian_matrix))
