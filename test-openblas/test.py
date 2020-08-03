#! /usr/bin/env python

import numpy as np


diff = np.load('diff.npy')
mass = np.load('mass.npy')

print(diff.dot(mass))
