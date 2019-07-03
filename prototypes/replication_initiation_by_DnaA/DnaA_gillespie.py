#! /usr/bin/python
"""
Runs a stochastic Gillespie simulation of DnaA-mediated replication initiation
in E. coli.
"""

# Imports
from __future__ import division, print_function, absolute_import

import os

import numpy as np
import matplotlib.pyplot as plt

# Parameters
# TODO: Add parameters from literature and model runs


class DnaA_gillespie(object):
	def __init__(self):
		pass


if __name__ == '__main__':
	sim = DnaA_gillespie()
	sim.run()
	sim.plot_results()
