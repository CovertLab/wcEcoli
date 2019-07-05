"""Generic Sisyphus/Gaia workflow builder."""

from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import multiprocessing
import os
from typing import Any, Dict, Iterable, List, Optional

from gaia.client import Gaia
from requests import ConnectionError

from wholecell.utils import filepath as fp


# Config details to pass to Gaia.
# ASSUMES: An ssh tunnel is open. See runscripts/sisyphus/ssh-tunnel.sh.
# ASSUMES: /etc/hosts has the line:
#   127.0.0.1   zookeeper-prime
GAIA_CONFIG = {'gaia_host': 'localhost:24442', 'kafka_host': 'localhost:9092'}


def _rebase(path, old_prefix, new_prefix):
	# type: (str, str, str) -> str
	"""Return a path starting with new_prefix in place of old_prefix."""
	new_path = os.path.join(new_prefix, os.path.relpath(path, old_prefix))

	### Should new_prefix end with os.sep iff path does? relpath() drops the
	### trailing os.sep. This `if` statement will reattach it but disable it
	### until we figure out what to do.
	# if path.endswith(os.sep) and not new_path.endswith(os.sep):
	# 	new_path = os.path.join(new_path, '')

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

def _copy_path_list(value):
	# type: (List[str]) -> List[str]
	"""Copy a list and check that it's a list of absolute paths to catch goofs
	like `outputs=plot_dir` (instead of `outputs=[plot_dir]`) or
	`outputs=['out/metadata.json']`. Sisyphus needs absolute paths to mount
	into the Docker container.
	"""
	assert isinstance(value, list), 'Expected a list, not {}'.format(repr(value))
	for path in value:
		assert os.path.isabs(path), 'Expected a absolute path, not {}'.format(path)
	return list(value)

def _launch_sisyphus(worker_name):
	# type: (str) -> None
	os.system("runscripts/sisyphus/launch-sisyphus.sh {}".format(worker_name))


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
		self.inputs  = _copy_path_list(kwargs.get('inputs',  []))  # type: List[str]
		self.outputs = _copy_path_list(kwargs.get('outputs', []))  # type: List[str]
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
		# type: (str) -> None
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

	def write(self):
		# type: () -> None
		"""Build the workflow and write it as JSON files for debugging that can
		be manually sent to the Gaia server.
		"""
		commands = self.get_commands()
		processes = self.get_processes()

		fp.makedirs('out')
		commands_path = os.path.join('out', 'wcm-commands.json')
		processes_path = os.path.join('out', 'wcm-processes.json')

		self.log_info('\nWriting {}, {}'.format(commands_path, processes_path))
		fp.write_json_file(commands_path, commands)
		fp.write_json_file(processes_path, processes)

	def launch_workers(self, count):
		# type: (int) -> None
		if count <= 0:
			return

		self.log_info('\nLaunching {} worker node(s).'.format(count))
		user = os.environ['USER']
		names = ['{}-{}'.format(user, i) for i in range(count)]

		pool = multiprocessing.Pool(10)
		pool.map(_launch_sisyphus, names)

	def send(self, worker_count=4):
		# type: (int) -> None
		"""Build the workflow and send it to the Gaia server to start running."""
		# TODO(jerry): Gaia will soon launch the workers, maybe with a given
		#  worker_count.
		self.launch_workers(worker_count)

		commands = self.get_commands()
		processes = self.get_processes()

		gaia = Gaia(GAIA_CONFIG)

		try:
			self.log_info('\nUploading {} tasks to Gaia'.format(len(commands)))
			gaia.command(commands)
			gaia.merge('sisyphus', processes)
			gaia.trigger('sisyphus')
		except ConnectionError as e:
			print('\n*** Did you set up port forwarding for gaia-base and'
				  ' zookeeper-prime? See runscripts/sisyphus/ssh-tunnel.sh ***\n')
			raise e
