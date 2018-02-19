import numpy as np
from wholecell.tests.utils.test_library_performance import time_it


def dot(x, y):
	for i in range(1000):
		np.dot(x, y)

def dot_convert(x, y):
	for i in range(1000):
		np.array(np.dot(x, y), np.int)

print 'int x int:'
int1 = np.array(np.random.rand(10000) * 100, np.int)
int2 = np.array(np.random.rand(10000) * 100, np.int)
time_it(lambda: dot(int1, int2))
print '\nfloat x float'
float1 = np.array(np.random.rand(10000) * 100, np.float)
float2 = np.array(np.random.rand(10000) * 100, np.float)
time_it(lambda: dot(float1, float2))
print '\nint x float'
ints = np.array(np.random.rand(10000) * 100, np.int)
floats = np.array(np.random.rand(10000) * 100, np.float)
time_it(lambda: dot(ints, floats))
print '\nfloat x float (converted)'
float1 = np.array(np.random.rand(10000) * 100, np.float)
float2 = np.array(np.random.rand(10000) * 100, np.float)
time_it(lambda: dot_convert(float1, float2))