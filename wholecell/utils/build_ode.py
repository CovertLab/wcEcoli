from __future__ import absolute_import, division, print_function

import numpy as np
from numba import njit
from typing import Callable


def build_function(arguments, expression, jit=True):
	# type: (str, str, bool) -> Callable
	"""Build a function from its arguments and source code expression, give it
	access to `Numpy as np`, and set up Numba to JIT-compile it on demand.

	Numba will optimize expressions like 1.0*y[2]**1.0 while compiling it
	to machine code.

	Args:
		arguments (str): comma-separated lambda argument names
		expression (str): expression to compile
		jit (bool): whether to JIT-compile the function; this option is just to
			work around Numba compiler bugs

	Returns:
		a Numba Dispatcher function(arguments)
	"""
	f = eval('lambda {}: {}'.format(arguments, expression), {'np': np}, {})

	if jit:
		# Too bad cache=True doesn't work with string source code.
		f_jit = njit(f, error_model='numpy')
		return f_jit
	return f


# TODO(jerry): Surely we can extract the argument of "Matrix([arg])" using
#  sympy calls more reliably than str(expr)[7:-1].

def derivatives(derivatives, jit=True):
	"""Build an optimized derivatives ODE function(y, t)."""
	return build_function('y, t',
		'np.array(' + str(derivatives)[7:-1] + ').reshape(-1)',
		jit)

def derivatives_jacobian(derivatives_jacobian, jit=True):
	"""Build an optimized derivatives ODE Jacobian function(y, t)."""
	return build_function('y, t',
		'np.array(' + str(derivatives_jacobian)[7:-1] + ')',
		jit)

def derivatives_with_rates(derivatives, jit=True):
	"""Build an optimized derivatives ODE function(y, t, kf, kr)."""
	return build_function('y, t, kf, kr',
		'np.array(' + str(derivatives)[7:-1] + ').reshape(-1)',
		jit)

def derivatives_jacobian_with_rates(derivatives_jacobian, jit=True):
	"""Build an optimized derivatives ODE Jacobian function(y, t, kf, kr)."""
	return build_function('y, t, kf, kr',
		'np.array(' + str(derivatives_jacobian)[7:-1] + ')',
		jit)
