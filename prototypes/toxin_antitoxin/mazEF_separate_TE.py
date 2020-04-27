
import numpy as np
import matplotlib.pyplot as plt

def run_model():

	# initial conditions
	DNA = 1
	DNA_bound = 0
	mRNA = 0
	mazE = 2
	mazF = 4
	mazEF = 0
	t = 0
	division_time = [0]

	# params
	mRNA_deg = np.log(2)/(2 * 60)
	mazE_deg = np.log(2)/(30*60)
	mazF_deg = np.log(2)/(4*60*60)
	c_EFassoc = np.log(2)/60
	c_EFdiss = np.log(2)/30
	c_DNAassoc = np.log(2)/(25) #  taken from hipAB rate
	c_DNAdiss = np.log(2)/(1)  # taken from hipAB rate
	c_trsc = 0.0005
	c_trl_max = 1
	mazF_inhib = 1 + 20/mazF
	c_trl_mazE = (0.101 * c_trl_max)/mazF_inhib # 1.01 from model
	c_trl_mazF = (0.039 * c_trl_max)/mazF_inhib # .39 from model


	# initialize
	steps = 100000
	out = np.ndarray((6,steps))
	t_vec = []
	a = np.ndarray((10))

	for step in range(0,steps):
		# dynamic rates
		if mazF > 0:
			mazF_inhib = 1 + mazF/30
		else:
			mazF_inhib = 1
				# reduce trl rate
		if step > 35000:
			c_trl_max = 0.1/mazF_inhib
		else:
			c_trl_max = 1/mazF_inhib

		c_trl_mazE = 0.101 * c_trl_max # 1.01 from model
		c_trl_mazF = 0.039 * c_trl_max # .39 from model

		# TODO change to a dictionary to avoid confusion
		a[0] = c_trsc * DNA
		a[1] = mRNA_deg * mRNA
		a[2] = c_trl_mazE * mRNA
		a[3] = c_trl_mazF * mRNA
		a[4] = mazE_deg * mazE
		a[5] = mazF_deg * mazF
		a[6] = c_EFassoc * mazE * mazF
		a[7] = c_EFdiss * mazEF
		a[8] = c_DNAassoc * mazEF * DNA
		a[9] = c_DNAdiss * DNA_bound

		atot = sum(a)

		r1 = np.random.random(1)
		r2 = np.random.random(1)

		tau = 1/atot * np.log(1/r1)
		pcut = atot * r2

		# choose reaction
		for n in range(len(a)):
			if np.cumsum(a)[n] > pcut:
				q = n
				break

		# execute reaction
		if q == 0:
			mRNA += 1
		elif q == 1:
			mRNA -= 1
		elif q == 2:
			mazE += 1
		elif q == 3:
			mazF += 1
		elif q == 4:
			mazE -= 1
		elif q == 5:
			mazF -= 1
		elif q == 6:
			if (mazE >= 2) & (mazF >= 4):
				mazEF += 1
				mazE -= 2
				mazF -= 4
		elif q == 7:
			mazEF -= 1
			mazE += 2
			mazF += 4
		elif q == 8:
			DNA_bound += 1
			DNA -= 1
			mazEF -= 1
		elif q == 9:
			DNA_bound -= 1
			DNA += 1
			mazEF += 1
		

		t += tau
		# cell division
		# if t > division_time[-1] + 1800:
		# 	division_time.append(t[0])
		# 	mazE = divide_molecule(mazE)
		# 	mazF = divide_molecule(mazF)
		# 	mazEF = divide_molecule(mazEF)
		# 	mRNA = divide_molecule(mRNA)

		# store outputs
		out[0, step] = mRNA
		out[1, step] = mazE
		out[2, step] = mazF
		out[3, step] = mazEF
		out[4, step] = DNA
		out[5, step] = DNA_bound
		t_vec.append(t[0])


	gap = 100
	plt.subplot(221)
	plt.plot(t_vec[0:-1:gap], out[0,0:-1:gap], label='mRNA')
	plt.legend()
	plt.subplot(222)
	plt.plot(t_vec[0:-1:gap], out[4,0:-1:gap], label='DNAfree')
	plt.legend()

	# plt.subplot(233)
	# plt.plot(t_vec[0:-1:gap], out[5,0:-1:gap], label='DNAbound')
	# plt.legend()





	plt.subplot(223)
	plt.plot(t_vec[0:-1:gap], out[1,0:-1:gap], label='mazE')
	plt.plot(t_vec[0:-1:gap], out[2,0:-1:gap], label='mazF')
	# plot_cell_division_lines(division_time)
	plt.axvline(t_vec[35000], c='r')
	plt.title('Protein')
	plt.legend()

	plt.subplot(224)
	plt.plot(t_vec[0:-1:gap], out[3,0:-1:gap])
	plt.title('mazEF Complex')
	plt.tight_layout()
	print("plot finished")
	plt.show()

# -------------------------------------------
def divide_molecule(prot):
	prot, r = np.divmod(prot, 2)
	prot += r * np.random.randint(2) # randomly assign remaining extra protein
	return prot

def plot_cell_division_lines(division_time):
	for l in division_time:
		plt.axvline(l, color='k')

# -------------------------------------------

if __name__ == '__main__':
	run_model()
