#! /usr/bin/env python
"""
Script to repeatedly run the qlaunch rapidfire command, with an option to
specify the maximum number of jobs in the queue.
"""

import os
import argparse

VERBOSE=False

def run_qlaunch(max_jobs):
	while True:
		os.system(
			f"qlaunch -r rapidfire --nlaunches infinite --sleep 10 "
			f"--maxjobs_queue {max_jobs}")

		if VERBOSE:
			print("qlaunch aborted. Restarting...")

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"max_jobs",
		help="Number of maximum jobs to allow",
		type=int,
		default=64,
		)

	args = parser.parse_args().__dict__

	run_qlaunch(args["max_jobs"])
