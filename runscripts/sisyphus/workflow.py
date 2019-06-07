"""Workflow builder.

TODO(jerry): Move this module to a Sisyphus Python client library.
"""

from __future__ import absolute_import, division, print_function

import json
import os
from typing import Any, Dict, List, Optional, Set


class Task(object):
	"""An AST element."""

	def __init__(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> None
		"""Construct a Workflow Task from the Task args: `name`, Docker `image`,
		`inputs` (dependencies) and `outputs` path mappings, and `commands`.

		The `upstream_tasks` and the `>>` operator are just convenient ways to
		add `inputs`.

		If present, `kwargs['storage_prefix']` is a storage bucket:path to
		prepend to all the remote input and output paths.
		TODO(jerry): Keep the storage_prefix feature?
		"""
		prefix = kwargs.get('storage_prefix', '')
		def fix(subpath):
			return os.path.join(prefix, subpath)

		self.name = kwargs['name']
		self.image = kwargs['image']
		self.inputs  = {fix(k): v for k, v in kwargs.get('inputs',  {}).viewitems()}
		self.outputs = {fix(k): v for k, v in kwargs.get('outputs', {}).viewitems()}
		self.commands = kwargs['commands']

		# Sisyphus uses inputs and outputs to determine data flow dependencies.
		# self.upstream_tasks is just for visualizing the DAG.
		# TODO(jerry): Skip this? Keep a dict of input paths to Tasks instead?
		self.upstream_tasks = set()  # type: Set[Task]

		for task in upstream_tasks:
			task >> self

	def __rshift__(self, t2):
		# type: (Task) -> Task
		"""Set downstream: `t1 >> t2` makes `t2` depend on `t1`'s outputs.
		Return `t2` for chaining.
		"""
		t2.upstream_tasks.add(self)
		t2.inputs.update(self.outputs)
		return t2

	def _set_path_mapping(self, io_map, storage_prefix, *path_elements):
		# type: (Dict[str, str], str, *str) -> None
		path = os.path.join(*path_elements)
		io_map[os.path.join(storage_prefix, path)] = path

	def set_input_mapping(self, storage_prefix, *path_elements):
		# type: (str, *str) -> None
		"""Set an input mapping storage_prefix:path --> path.
		NOTE: A path ending with '/' is treated as an entire directory.
		"""
		self._set_path_mapping(self.inputs, storage_prefix, *path_elements)

	def set_output_mapping(self, storage_prefix, *path_elements):
		# type: (str, *str) -> None
		"""Set an input mapping storage_prefix:path <-- path.
		NOTE: A path ending with '/' is treated as an entire directory.
		"""
		self._set_path_mapping(self.outputs, storage_prefix, *path_elements)

	def task_doc(self):
		# type: () -> Dict[str, Any]
		"""Return the Sisyphus task document to run this Task."""
		fields = vars(self)
		return {key: fields[key] for key in (
			'name', 'image', 'inputs', 'outputs', 'commands')}


class Workflow(object):
	"""A Sisyphus workflow builder. It's essentially an AST."""

	def __init__(self, verbose_logging=True):
		# type: (bool) -> None
		self.verbose_logging = verbose_logging
		self._tasks = {}  # type: Dict[str, Task]

	def log_info(self, message):
		if self.verbose_logging:
			print(message)

	def get(self, task_name, default=None):
		# type: (str, Optional[Task]) -> Optional[Task]
		return self._tasks.get(task_name, default)

	def __getitem__(self, task_name):
		# type: (str) -> Task
		return self._tasks[task_name]

	def add_task(self, task):
		# type: (Task) -> Task
		"""Add a task object. Return it for chaining."""
		self._tasks[task.name] = task
		self.log_info('    Added task: {}'.format(task.name))
		return task

	def as_dag(self):
		# type: () -> List[Dict[str, Any]]
		"""Return the workflow DAG as JSON-serializable data."""
		return [task.task_doc() for task in self._tasks.viewvalues()]

	def as_json(self):
		# type: () -> str
		"""Return the workflow DAG in JSON format."""
		dag = self.as_dag()
		result = json.dumps(dag, ensure_ascii=False, indent=4,
			separators=(',', ': '), sort_keys=True) + '\n'
		return result

	def enqueue(self):  # TODO(jerry): "start()"? "run()"? "send()"?
		# type: () -> None
		"""Construct a workflow and enqueue it on the work servers."""
		json_dag = self.as_json()

		# TODO(jerry): *** POST it to the workflow manager via HTTP. ***
		print()
		print(json_dag)  # meanwhile, just to see something happen
