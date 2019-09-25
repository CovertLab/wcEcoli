import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import copy

def model(X, t, params):
	A = X[0]
	T = X[1]

	kt = params['kt']
	n = params['n']
	kp1 = params['kp1']
	kp2 = params['kp2']
	kh = params['kh']
	p = params['p']
	sigma = params['sigma']
	alpha = params['alpha']
	umax = params['umax']
	lambda_A = params['lambda_A']
	lambda_T = params['lambda_T']

	X = 1 + (T**n)/(kt**n)
	Y = 1 + (A**2)/(kp1 ** 2) + ((2 * A**2 * T) / (kp2**2 * kh))**p + (A**2 * T**2)/(kp1**2 * kh**2)

	dAdt = (sigma * alpha) / (X * Y) - umax/X * A - lambda_A * A 
	dTdt = alpha / (X * Y) - umax/X * T - lambda_T * T

	return [dAdt, dTdt]
#---------------------------------------------------------------------

def hysteresis_graph(model, param_name, param_range, AT_inital, color='r'):
	p = {
			'A' : AT_inital[0],
			'T' : AT_inital[1],
			'lambda_A' : (np.log(2)/60),
			'lambda_T' : np.log(2)/(48*60),
			'sigma' : 10.0,
			'alpha' : 1.0,
			'kh' : 100.0,
			'kp1' : 1000.0,
			'kp2' : 10.0,
			'kt' : 10.0,
			'd' : 0, # noise scaling parameter
			'p' : 2.0,
			'n' : 2.5, 
			'umax' : (np.log(2)/30)} 

	t = np.linspace(0,10000.0,1000)
	# ---------------------- reproduce figure 3 from 2013 paper
	p_mod = copy.deepcopy(p)
	param_ratio = []
	T = []
	for m in np.linspace(param_range[0], param_range[1] , 20):
		p_mod[param_name] = m * state[param_name]
		x = odeint(model, AT_inital, t, args=(p_mod,))
		param_ratio.append(m)
		T.append(x[:,1])

	plt.plot(param_ratio, T, color=color)

#--------------------------------------------------------------------

p = {
		'lambda_A' : (np.log(2)/60),
		'lambda_T' : np.log(2)/(48*60),
		'sigma' : 10.0,
		'alpha' : 1.0,
		'kh' : 100.0,
		'kp1' : 1000.0,
		'kp2' : 10.0,
		'kt' : 10.0,
		'd' : 0, # noise scaling parameter
		'p' : 2.0,
		'n' : 2.5, 
		'umax' : 0.75*(np.log(2)/30),  # umax/umax0 = 0.75
		'dTdt' : 1,
		'dAdt' : 1} 

# loop through and vary A
A0 = np.linspace(0,200,25)

ode_out = []
for a0 in A0:
	x = odeint(model, [a0,200], np.linspace(0,5000.0,100), args=(p,))
	ode_out.append(x)


# plt.plot(np.linspace(0,2000,100), x[:,0])
# plt.savefig('plots/odeint_solution_A.png')
# plt.clf()

# plt.plot(np.linspace(0,2000,100), x[:,1])
# plt.savefig('plots/odeint_solution_T.png')
# plt.clf()

plt.figure()
axs = plt.gca()

points = np.array([x[:,0], x[:,1]]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap='rainbow')
lc.set_array(np.linspace(0,1,x[:,0].shape[0]))
lc.set_linewidth(3)
line = axs.add_collection(lc)

# plt.plot(x[:,0], x[:,1])
# plt.savefig('plots/odeint_solution_AvsT.png')
plt.xlim((0,80))
plt.ylim((0,80))
plt.xlabel('A')
plt.ylabel('T')

plt.savefig('plots/odeint_solution_AvsT.png')
plt.clf()