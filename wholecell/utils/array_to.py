"""
Defines a function array_to, that takes an array of keys and
a second array of data to create and return one dictionary.
"""

def array_to(keys, array):
    return {
        key: array[index]
        for index, key in enumerate(keys)}