import numpy as np
import matplotlib.pyplot as plt
# eqs
# mRNA = c_trsc * [DNA] - mRNA_deg * [mRNA]
# mazE = c_trl * [mRNA] - mazE_deg * [mazE]
# mazF = c_trl * [mRNA] - mazF_deg * [mazF]

# params
mRNA_deg = np.log(2)/(2 * 60)
mazE_deg = np.log(2)/(30*60)
mazF_deg = np.log(2)/(4*60*60)
c_EFassoc = 0.002
c_EFdiss = 0.2
c_DNAassoc = np.log(2)/(25 * 60) #  taken from hipAB rate
c_DNAdiss = np.log(2)/(60)  # taken from hipAB rate
c_trsc = 0.0005
c_trl_mazE = 0.101
c_trl_mazF = 0.039
# initial conditions
DNA = 1
DNA_bound = 0
mRNA = 0
mazE = 0
mazF = 0
mazEF = 0
t = 0

# initialize
steps = 50000
out = np.ndarray((6,steps))
t_vec = []
a = np.ndarray((9))

for step in range(0,steps):
	a[0] = c_trsc * DNA
	a[1] = mRNA_deg * mRNA
	a[2] = c_trl * mRNA
	a[3] = mazE_deg * mazE
	a[4] = mazF_deg * mazF
	a[5] = c_EFassoc * mazE * mazF
	a[6] = c_EFdiss * mazEF
	a[7] = c_DNAassoc * mazEF * DNA
	a[8] = c_DNAdiss * DNA_bound

	atot = sum(a)

	r1 = np.random.random(1)
	r2 = np.random.random(1)

	tau = 1/atot * np.log(1/r1)
	pcut = atot * r2

	if a[0] > pcut:
		mRNA += 1
	elif sum(a[0:2]) > pcut:
		mRNA -= 1
	elif sum(a[0:3]) > pcut:
		mazE += 1
		mazF += 1
	elif sum(a[0:4]) > pcut:
		mazE -= 1
	elif sum(a[0:5]) > pcut:
		mazF -= 1
	elif sum(a[0:6]) > pcut:
		mazEF += 1
		mazE -= 2
		mazF -= 4
	elif sum(a[0:7]) > pcut:
		mazEF -= 1
		mazE += 2
		mazF += 4
	elif sum(a[0:8]) > pcut:
		DNA_bound += 1
		DNA -= 1
	elif sum(a[0:9]) > pcut:
		DNA_bound -= 1
		DNA += 1

	t += tau

	out[0, step] = mRNA
	out[1, step] = mazE
	out[2, step] = mazF
	out[3, step] = mazEF
	out[4, step] = DNA
	out[5, step] = DNA_bound
	t_vec.append(t[0])


gap = 50
plt.subplot(231)
plt.plot(t_vec[0:-1:gap], out[0,0:-1:gap], label='mRNA')
plt.legend()
plt.subplot(232)
plt.plot(t_vec[0:-1:gap], out[4,0:-1:gap], label='DNAfree')
plt.legend()

plt.subplot(233)
plt.plot(t_vec[0:-1:gap], out[5,0:-1:gap], label='DNAbound')
plt.legend()





plt.subplot(234)
plt.plot(t_vec[0:-1:gap], out[1,0:-1:gap], label='mazE')
plt.plot(t_vec[0:-1:gap], out[2,0:-1:gap], label='mazF')
plt.title('Protein')
plt.legend()

plt.subplot(235)
plt.plot(t_vec[0:-1:gap], out[3,0:-1:gap])
plt.title('mazEF Complex')
plt.tight_layout()
plt.show()

