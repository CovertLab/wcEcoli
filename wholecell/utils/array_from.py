"""
Defines a function array_from, that takes in a dictionary and
returns an array of the values in the given dictionary
"""
import numpy as np

def array_from(d):
    return np.array(list(d.values()))