# cython: language_level=3str
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
Test case for #1399, where passing a float from Cython code (but not from Python
code) to a function that takes a np.npy_intp (like weird_func() below or
np.random.RandomState.multinomial()) raises
    TypeError: 'float' object cannot be interpreted as an integer
"""

import numpy as np
cimport numpy as np

def func(int length):
    # func(1.1) returns 1
    return length

def weird_func(np.npy_intp length):
    # weird_func(1.1) raises TypeError: 'float' object cannot be interpreted as an integer
    return length
