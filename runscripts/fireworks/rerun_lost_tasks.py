#! /usr/bin/env python
"""
Script to rerun fireworks that get preempted when submitted to the owners node
on sherlock.  Preempted jobs do not have a way of updating the status in mongoDB
to FIZZLED and will be stuck on RUNNING.  This checks the name of RUNNING tasks
against those in slurm to see if any RUNNING tasks are not actually running.
A task is only rerun after it is missing from slurm in MISSES_BEFORE_RERUN
consecutive checks, since a single squeue call can miss queued jobs and
rerunning a task whose job is still pending leaves its children stuck.

For a RUNNING task missing from slurm whose launch COMPLETED is refreshed rather
than rerun: its job finished but the final workflow refresh timed out on the
workflow lock ("fw_id N locked. Can't refresh!" at the end of its log), which
is common for large workflows, where every task shares one lock and a refresh
can take minutes. Raising WFLOCK_EXPIRATION_SECS in FW_config.yaml makes
this rarer.

It also clears workflow locks left by a job killed while holding one:
whenever a rerun, refresh or release fails on a locked workflow (FireWorks
already waited WFLOCK_EXPIRATION_SECS for it), and when a lock is seen in two
consecutive checks with no firework job in slurm. Otherwise every task
completion, refresh and rerun in the workflow fails on the lock and the run
stalls.

It also releases simulation tasks left WAITING although every parent
COMPLETED. The workflow document keeps a stale parent state, often with a
stale RESERVED launch on the child.

Use the option -f to automatically set these tasks to rerun or release them,
otherwise only the fw_ids will be printed for manually fixing (releasing) them.

TODO:
- add a way of catching sim issues that cause tasks to be stuck on RUNNING
(MAX_RERUN is used to prevent the same task from being rerun many times if the
issue comes from the code)
- parse all options from command line for more flexibility
"""

import json
import os
import sys
import time
from typing import Any, Dict, Set

from fireworks import LaunchPad
from fireworks.core.launchpad import LockedWorkflowError, WFLock

from wholecell.utils import filepath


VERBOSE = False
MAX_RERUN = 5
MISSES_BEFORE_RERUN = 2
CACHE_WAIT = 60
SLEEP = 300
LAUNCHPAD_FILE = 'my_launchpad.yaml'
RUN_CMD = len(sys.argv) > 1 and sys.argv[1] == '-f'


def get_fws(criteria):
	result = json.loads(filepath.run_cmdline(f'lpad get_fws {criteria}', fallback='[]'))

	# Handle single fw case to always make a list of dicts
	if isinstance(result, dict):
		result = [result]
	return result

def find_stuck_waiting(lp):
	"""Return {fw_id: name} of simulation tasks that are WAITING although
	every parent COMPLETED."""
	# Load every task state and workflow once; each workflow document lists
	# all of its tasks, so querying it per task is slow for a large run
	fw_state = {d['fw_id']: d['state'] for d in lp.fireworks.find({}, {'fw_id': 1, 'state': 1})}
	parent_links = {}
	for wf in lp.workflows.find({}, {'nodes': 1, 'parent_links': 1}):
		for node in wf['nodes']:
			parent_links[node] = wf['parent_links'].get(str(node), [])

	stuck = {}
	for fw in lp.fireworks.find({'state': 'WAITING', 'name': {'$regex': '^Simulation'}},
			{'fw_id': 1, 'name': 1}):
		parents = parent_links.get(fw['fw_id'], [])
		if parents and all(fw_state.get(p) == 'COMPLETED' for p in parents):
			stuck[fw['fw_id']] = fw['name']
	return stuck

def launch_completed(lp, fw_id):
	"""Whether a current (not archived) launch of the task COMPLETED."""
	fw = lp.fireworks.find_one({'fw_id': fw_id}, {'launches': 1})
	return lp.launches.count_documents(
		{'launch_id': {'$in': fw.get('launches', [])}, 'state': 'COMPLETED'}) > 0

def refresh(lp, fw_id):
	"""Recompute the task and its children's states from their launches, as
	FireWorks does at the end of a launch, and return the task's new state."""
	with WFLock(lp, fw_id):
		wf_obj = lp.get_wf_by_fw_id_lzyfw(fw_id)
		lp._update_wf(wf_obj, wf_obj.refresh(fw_id))
	return lp.fireworks.find_one({'fw_id': fw_id}, {'state': 1})['state']

def release_waiting(lp, fw_id):
	"""Mark the parents COMPLETED in the workflow document, cancel stale
	reservations and refresh the task so it becomes READY."""
	wf = lp.workflows.find_one({'nodes': fw_id}, {'parent_links': 1, 'fw_states': 1})
	parents = wf['parent_links'].get(str(fw_id), [])
	with WFLock(lp, fw_id):
		for p in parents:
			if wf['fw_states'].get(str(p)) != 'COMPLETED':
				lp.workflows.update_one({'nodes': fw_id}, {'$set': {f'fw_states.{p}': 'COMPLETED'}})
	for launch in lp.launches.find({'fw_id': fw_id, 'state': 'RESERVED'}, {'launch_id': 1}):
		lp.cancel_reservation(launch['launch_id'])
	return refresh(lp, fw_id)

def clear_locks(lp, reason):
	"""
	Clear every workflow lock with `lpad admin unlock` when a rerun failed.
	Only called after FireWorks already waited WFLOCK_EXPIRATION_SECS for a
	lock, so the lock is almost certainly stale (left by a job killed while
	holding it).
	"""
	n = lp.workflows.update_many({'locked': {'$exists': True}}, {'$unset': {'locked': ''}}).modified_count
	print(f'Cleared {n} workflow lock(s): {reason}')

lp = LaunchPad.from_file(LAUNCHPAD_FILE)
rerun_counts = {}  # type: Dict[str, int]
miss_counts = {}  # type: Dict[str, int]
stuck_before = set()  # type: Set[int]
locked_before = set()  # type: Set[Any]
while True:
	# Get all fireworks with running status and all active job names
	running = get_fws('-s RESERVED') + get_fws('-s RUNNING')
	if VERBOSE:
		print(running)
	time.sleep(CACHE_WAIT)  # squeue is sometimes cached so wait until we know the new tasks will be returned
	if (squeue := filepath.run_cmdline(f'squeue -u {os.environ["USER"]} -o %j')) is None:
		continue
	else:
		tasks = set(squeue.split('\n'))
	if VERBOSE:
		print(tasks)

	# Clear workflow locks left by a job killed while holding one. Every task
	# completion, refresh and rerun in the workflow waits on the lock, so a
	# stale one stalls the whole run. A lock is only cleared when it was seen
	# in two consecutive checks with no firework job in slurm, since then
	# nothing can be holding it legitimately (a refresh takes minutes at most).
	fw_names = set(lp.fireworks.distinct('name'))
	no_fw_jobs = not (tasks & fw_names)
	locked = {wf['_id'] for wf in lp.workflows.find({'locked': {'$exists': True}}, {'_id': 1})}
	for wf_id in sorted(locked & locked_before) if no_fw_jobs else []:
		if RUN_CMD:
			lp.workflows.update_one({'_id': wf_id, 'locked': {'$exists': True}}, {'$unset': {'locked': ''}})
			print(f'Cleared stale lock on workflow {wf_id} (locked in 2 checks with no jobs in slurm)')
		else:
			print(f'Would clear stale lock on workflow {wf_id} if -f passed as arg')
	locked_before = locked if no_fw_jobs else set()

	# Count consecutive checks in which each running firework is missing from
	# the slurm queue
	missing = {str(fw['fw_id']) for fw in running if fw['name'] not in tasks}
	miss_counts = {fw_id: miss_counts.get(fw_id, 0) + 1 for fw_id in missing}
	fw_ids = ','.join(fw_id for fw_id, n in miss_counts.items() if n >= MISSES_BEFORE_RERUN)
	if VERBOSE:
		print(miss_counts)

	if fw_ids:
		# Get the current status of fireworks in case they finished between
		# running get_fws and squeue calls
		current = get_fws(f'-i {fw_ids}')
		if VERBOSE:
			print(current)

		# Rerun fireworks that still have a running status but are not in slurm
		fws = []
		for fw in current:
			# A launch that completed only missed its final workflow refresh
			# (the lock timed out), so refresh it instead of rerunning
			if fw['state'] == 'RUNNING' and launch_completed(lp, fw['fw_id']):
				if RUN_CMD:
					try:
						print(f'Refreshing fw_id {fw["fw_id"]} (launch completed, '
							f'refresh missed): now {refresh(lp, fw["fw_id"])}')
						miss_counts.pop(str(fw['fw_id']), None)
					except LockedWorkflowError:
						clear_locks(lp, f'workflow locked while refreshing fw_id {fw["fw_id"]}; will retry on the next check')
				else:
					print(f'Would refresh fw_id {fw["fw_id"]} (launch completed) if -f passed as arg')
				continue
			if fw['state'] == 'RUNNING' or fw['state'] == 'RESERVED':
				fw_id = str(fw['fw_id'])
				# Only count actual reruns so a dry run (no -f) would otherwise
				# drop a task from the report after MAX_RERUN checks:
				if RUN_CMD:
					rerun_counts[fw_id] = rerun_counts.get(fw_id, 0) + 1
				if rerun_counts.get(fw_id, 0) <= MAX_RERUN:
					fws.append(fw_id)
		fws_to_rerun = ','.join(fws)
		if fws_to_rerun:
			if RUN_CMD:
				print(f'Rerunning fw_ids {fws_to_rerun} (missing from squeue in {MISSES_BEFORE_RERUN} checks)')
				result = filepath.run_cmdline(f'lpad rerun_fws -i {fws_to_rerun}', input_='y', timeout=None)
				if result is None:
					clear_locks(lp, 'rerun failed (usually a locked workflow); will retry on the next check')
				else:
					for fw_id in fws:
						miss_counts.pop(fw_id, None)
			else:
				print(f'Would rerun fw_ids {fws_to_rerun} if -f passed as arg. Manually rerun with:'
					f'\n\tlpad rerun_fws -i {fws_to_rerun}')

	# Release tasks found stuck WAITING in two consecutive checks (a parent
	# that just completed may still be refreshing its children)
	stuck = find_stuck_waiting(lp)
	for fw_id in sorted(stuck.keys() & stuck_before):
		if stuck[fw_id] in tasks:
			continue
		if RUN_CMD:
			try:
				print(f'Released stuck WAITING fw_id {fw_id} ({stuck[fw_id]}): now {release_waiting(lp, fw_id)}')
			except LockedWorkflowError:
				clear_locks(lp, f'workflow locked while releasing fw_id {fw_id}; will retry on the next check')
		else:
			print(f'Would release stuck WAITING fw_id {fw_id} ({stuck[fw_id]}) if -f passed as arg')
	stuck_before = set(stuck)

	# Sleep for a bit
	print(f'{time.ctime()}: sleeping for {SLEEP} s...')
	time.sleep(SLEEP)
