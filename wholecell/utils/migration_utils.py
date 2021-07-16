"""
Defines two helper functions to be used in migration.
"""
import numpy as np

def arrays_to(n, attrs):
    ds = []
    for index in np.arange(n):
        d = {}
        for attr in attrs.keys():
            d[attr] = attrs[attr][index]
        ds.append(d)

    return ds

def add_elements(elements, id):
    return {
        '_add': [{
            'key': element[id],
            'state': element}
            for element in elements]}