import json
import unum
import os
import numpy as np

def write_json(path, numpy_dict):
    INFINITY = float('inf')

    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif obj == INFINITY:
                return '__INFINITY__'
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, unum.Unum):
                return float(obj)
            elif isinstance(obj, np.bool_):
                return bool(obj)
            else:
                return super(NpEncoder, self).default(obj)

    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, 'w') as outfile:
        json.dump(numpy_dict, outfile, cls=NpEncoder)