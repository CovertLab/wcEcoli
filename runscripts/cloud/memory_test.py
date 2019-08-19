import time
import numpy as np

def grow():
	pile = []
	while True:
		time.sleep(1)
		pile.append(np.full((100000000,), 15))

if __name__ == '__main__':
	grow()
