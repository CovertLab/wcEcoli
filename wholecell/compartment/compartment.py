from __future__ import absolute_import
from __future__ import division

from wholecell.containers.bulk_objects_container import BulkObjectsContainer

class Compartment(object):
	def __init__(self):
		self.subpartments = []
		self.environment = {}
		self.states = {}
		self.processes = []
