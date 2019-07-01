"""Generic Sisyphus/Gaia workflow builder."""

from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import os
from typing import Any, Dict, Iterable, List, Optional

from wholecell.utils import filepath as fp


def _rebase(path, old_prefix, new_prefix):
	# type: (str, str, str) -> str
	"""Return a path starting with new_prefix in place of old_prefix."""
	new_path = os.path.join(new_prefix, os.path.relpath(path, old_prefix))
	if path.endswith(os.sep) and not new_path.endswith(os.sep):
		new_path = os.path.join(new_path, '')
	assert '..' not in new_path, "Rebasing a path that doesn't start with old_prefix?"
	return new_path

def _keyify(paths):
	# type: (Iterable[str]) -> Dict[str, str]
	"""Map sequential keys to the given paths."""
	return {str(i): path for i, path in enumerate(paths)}

def _re_keyify(paths, old_prefix, new_prefix):
	# type: (Iterable[str], str, str) -> Dict[str, str]
	"""Map sequential keys to rebased paths."""
	return _keyify([_rebase(path, old_prefix, new_prefix) for path in paths])

def _copy_list(value):
	# type: (List[str]) -> List[str]
	"""Copy a list. Also check that it's a list to catch a goof like
	`outputs=plot_dir` instead of `outputs=[plot_dir]`.
	"""
	if isinstance(value, list):
		return list(value)
	raise ValueError('Expected a list, not {}'.format(repr(value)))


class Task(object):
	"""A workflow task builder."""

	def __init__(self, upstream_tasks=(), **kwargs):
		# type: (Iterable[Task], **Any) -> None
		"""Construct a Workflow Task from the kwargs: key, image, commands,
		inputs, and outputs.

		upstream_tasks and the `>>` operator are convenient ways to add inputs.
		"""
		self.key = kwargs['key']  # type: str  # the task name
		self.image = kwargs['image']  # type: str
		self.commands = kwargs['commands']  # type: List[Dict[str, List[str]]]
		self.inputs  = _copy_list(kwargs.get('inputs',  []))  # type: List[str]
		self.outputs = _copy_list(kwargs.get('outputs', []))  # type: List[str]
		self.storage_prefix = kwargs.get('storage_prefix', '')
		self.local_prefix = kwargs.get('local_prefix', '')

		for task in upstream_tasks:
			task >> self

	def __rshift__(self, t2):
		# type: (Task) -> Task
		"""Set downstream: `t1 >> t2` adds `t1`'s outputs to `t2`'s inputs.
		Return `t2` for chaining.
		"""
		t2.inputs.extend(self.outputs)
		return t2

	def get_command(self):
		# type: () -> Dict[str, Any]
		"""Return a Sisyphus Command to run this Task."""
		return dict(
			key=self.key,
			image=self.image,
			commands=self.commands,
			inputs=_keyify(self.inputs),
			outputs=_keyify(self.outputs),
			vars={})

	def get_process(self):
		# type: () -> Dict[str, Any]
		"""Return a Sisyphus Process to run this Task."""
		return dict(
			key=self.key,
			command=self.key,
			inputs=_re_keyify(self.inputs, self.local_prefix, self.storage_prefix),
			outputs=_re_keyify(self.outputs, self.local_prefix, self.storage_prefix))


class Workflow(object):
	"""A workflow builder."""

	def __init__(self, namespace, verbose_logging=True):
		# type: (str, bool) -> None
		self.namespace = namespace
		self.verbose_logging = verbose_logging
		self._tasks = OrderedDict()  # type: Dict[str, Task]

	def log_info(self, message):
		if self.verbose_logging:
			print(message)

	def get(self, task_key, default=None):
		# type: (str, Optional[Task]) -> Optional[Task]
		return self._tasks.get(task_key, default)

	def __getitem__(self, task_key):
		# type: (str) -> Task
		return self._tasks[task_key]

	def add_task(self, task):
		# type: (Task) -> Task
		"""Add a task object. Return it for chaining."""
		# TODO(jerry): Workaround until we can set a namespace on all commands and processes.
		task.key = self.namespace + '-' + task.key

		self._tasks[task.key] = task
		self.log_info('    Added task: {}'.format(task.key))
		return task

	def get_commands(self):
		# type: () -> List[dict]
		"""Build this workflow's Commands."""
		return [task.get_command() for task in self._tasks.itervalues()]

	def get_processes(self):
		# type: () -> List[dict]
		"""Build this workflow's Processes."""
		return [task.get_process() for task in self._tasks.itervalues()]

	def send(self):
		# type: () -> None
		"""Build the workflow and send it to the work server."""
		commands = self.get_commands()
		processes = self.get_processes()

		self.log_info('\nWriting wcm-commands.json, wcm-processes.json')
		fp.write_json_file('wcm-commands.json', commands)
		fp.write_json_file('wcm-processes.json', processes)

		# TODO(jerry): *** Upload w/the namespace to the workflow manager. ***
