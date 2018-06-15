from line_profiler import LineProfiler
from functools import wraps

def line_profile(func):
	@wraps(func)
	def wrapper(*args, **kwargs):
		lp = LineProfiler()
		try:
			return lp(func)(*args, **kwargs)
		finally:
			lp.print_stats()

	return wrapper
