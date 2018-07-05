"""
Memory leak debugging utility.
"""

from __future__ import absolute_import
from __future__ import division

from contextlib import contextmanager
import gc
import os


# See https://pymotw.com/2/gc/ for more info on using the gc interface.
# DEBUG_UNCOLLECTABLE causes the collector to report on objects it can't
# collect. You need to combine it with DEBUG_OBJECTS to print info about
# objects and DEBUG_INSTANCES to print info about instances of old-style
# classes (not derived from object).
TRACE_UNCOLLECTABLES = gc.DEBUG_UNCOLLECTABLE | gc.DEBUG_OBJECTS | gc.DEBUG_INSTANCES
TRACE_NONE = 0


@contextmanager
def detect_leaks(enabled=None):
	"""A context manager that optionally detects Python object leaks in the
	`with` statement body.

	Set `enabled` to True to enable leak detection, False to disable leak
	detection, or default to let the 'DEBUG_GC' environment variable
	enable/disable leak detection.

	Leak detection has some overhead including running a full collection and
	printing a list of uncollectable objects.

	Per https://docs.python.org/2/library/gc.html, "Objects that have
	__del__() methods and are part of a reference cycle cause the entire
	reference cycle to be uncollectable, including objects not necessarily in
	the cycle but reachable only from it. Python doesn't collect such cycles
	automatically because, in general, it isn't possible for Python to guess
	a safe order in which to run the __del__() methods."
	"""
	if enabled is None:
		enabled = os.environ.get('DEBUG_GC', False)

	saved_debug_flags = gc.get_debug()
	gc.set_debug(TRACE_UNCOLLECTABLES if enabled else TRACE_NONE)

	yield  # yield to the `with` statement body

	if enabled:
		gc.collect()  # prints lines like "gc: uncollectable <CleanupGraph 0x10045f810>"
		# Examine the gc.garbage list here?

	gc.set_debug(saved_debug_flags)
