import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.backends.backend_pdf as mplpdf
import copy
import csv
import time

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

	dAdt = 1/(1 + T**n / kt**n) * (sigma * alpha)/(1 + (A**2/kp1**2) + \
		((2 * A**2 * T)/(kp2**2 * kh))**p + (A**2 * T**2)/(kp1**2 * kh**2)) \
		- umax/(1 + T**n/(kt**n)) * A - lambda_A * A

	dTdt = 1/(1 + T**n/kt**n) * alpha/( 1 + (A**2/kp1**2) + ((2 * A**2 * T)/(kp2**2 * kh))**p \
			+ (A**2 * T**2)/(kp1**2 * kh**2) ) - umax/(1 + T**n/kt**n) * T - lambda_T * T
	
	return [dAdt, dTdt]

def step_through_prange(param, p_mod, p, param_range):
	AT_initial = [1, 1]

	# forward from low state
	T1 = []
	param_ratio1 = []
	for m in param_range:
		p_mod[param] = m * p[param]
		x = odeint(model, AT_initial, t, args=(p_mod,))
		AT_initial =[x[-1,0], x[-1,1]]
		param_ratio1.append(m)
		T1.append(x[-1,1])

	# reverse from high state
	T2 = []
	param_ratio2 = []
	for m in param_range[::-1]:
		p_mod[param] = m * p[param]
		x = odeint(model, AT_initial, t, args=(p_mod,))
		AT_initial =[x[-1,0], x[-1,1]]
		param_ratio2.append(m)
		T2.append(x[-1,1])

	outdict = {'pr1' : param_ratio1, 'pr2' : param_ratio2, 'T1' : T1, 'T2' : T2}

	return 

def step_through_umax(p_mod, p):
	AT_initial = [1, 1]
	param_range = np.linspace(4, .001, 10)

	# forward from low state
	T1 = [None] * 10
	param_ratio1 = [None] * 10
	k = 0
	for m in param_range:
		p_mod['umax'] = m * p['umax']
		x = odeint(model, AT_initial, t, args=(p_mod,))
		AT_initial =[x[-1,0], x[-1,1]]
		param_ratio1[k] = m
		T1[k] = x[-1,1]
		k += 1

	# reverse from high state
	T2 = [None] * 10
	param_ratio2 = [None] * 10
	k = 0
	for m in param_range[::-1]:
		p_mod['umax'] = m * p['umax']
		x = odeint(model, AT_initial, t, args=(p_mod,))
		AT_initial =[x[-1,0], x[-1,1]]
		param_ratio2[k] = m
		T2[k] = x[-1,1]
		k += 1

	outdict = {'pr1' : param_ratio1, 'pr2' : param_ratio2, 'T1' : T1, 'T2' : T2}

	return outdict


#---------------------------------------------------------------------
outdir = '/Users/taryn/GoogleDrive/code/wcEcoli/prototypes/toxin_antitoxin/plots/multiD/'
# original parameters
p = {
		'A' : 1,
		'T' : 1,
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

# ranges to scan
N = 5 # number of slices in range

range_dict = {
				'A' : np.linspace(1, 200, N),
				'T' : np.linspace(1, 200, N),
				'lambda_A' : np.linspace(.01, 4, N),
				'lambda_T' : np.linspace(.01, 4, N),
				'sigma' : np.linspace(.1,3,N),
				'alpha' : np.linspace(1,5,N),
				'kh' : np.linspace(0.25, 4, N),
				'kp1' : np.linspace(0.1, 10, N),
				'kp2' : np.linspace(0.1, 3, N),
				'kt' : np.linspace(0.1, 5, N),
				'p' : np.linspace(1, 5, N),
				'n' : np.linspace(0.1, 1.5, N),
				'umax' : np.linspace(4, .001, 10)}

p_mod = copy.deepcopy(p)
t = np.linspace(0,50000,1000)
Tdat = step_through_prange('umax', p_mod, p, range_dict['umax'])

# print(Tdat['T1'])
# print(Tdat['T2'][::-1])
# print(sum(np.subtract(Tdat['T1'], Tdat['T2'][::-1])))

# set up csv for writing params
csv_header = ['set#','sigma', 'alpha', 'kh', 'lambdaA', 'lambdaT', 'kt', 'n']
with open( outdir + 'param_sets.csv', 'w') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	# writer.writerow(['set#'] + list(p_mod.keys()))
	writer.writerow(csv_header)

pdf = mplpdf.PdfPages(outdir + 'multiDscan.pdf')

setnum = 0 
runN = 0
timeleft = 'calculating...'

t0 = time.time()
toploopcount = 0

for n_sigma in range_dict['sigma']:
	p_mod['sigma'] = n_sigma * p['sigma']
	toploopcount += 1
	print('approx time left:' + timeleft + 'm')

	for n_alpha in range_dict['alpha']:
		p_mod['alpha'] = n_sigma * p['alpha']

		for n_kh in range_dict['kh']:
			p_mod['kh'] = n_kh * p['kh']

			for n_lambdaA in range_dict['lambda_A']:
				p_mod['lambda_A'] = n_lambdaA * p['lambda_A']

				for n_lambdaT in range_dict['lambda_T']:
					p_mod['lambda_T'] = n_lambdaT * p['lambda_T']

					for n_kt in range_dict['kt']:
						p_mod['kt'] = n_kt * p['kt']

						for n_n in range_dict['n']:
							p_mod['n'] = n_n * p['n']

							runN +=1
							# print('run # ' + str(runN))
							Tdat = step_through_umax(p_mod, p)
							sim_difference = np.sum( np.abs((np.subtract(Tdat['T1'], Tdat['T2'][::-1]))) > 10)

							if sim_difference > (0.2 * len(Tdat['T1'])):
								setnum += 1

								multiplier_list = [setnum, n_sigma, n_alpha, n_kh,  n_lambdaA, n_lambdaT, n_kt, n_n]


								with open(outdir + 'param_sets.csv', 'a') as csvfile:
									writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
									# writer.writerow([setnum] + list(p_mod.values()))
									writer.writerow(multiplier_list)

								# print params to console
								print('set# ' + str(setnum))	
								print(sim_difference)

								for key in p_mod:
									print(key + ': ' + "{:.4f}".format(p_mod[key]))
								print('---------------------')

								# make and save figure
								fig = plt.figure(figsize=(8,5))
								plt.subplot(1,2,1)
								plt.plot(Tdat['pr1'], Tdat['T1'], color='b')
								plt.plot(Tdat['pr2'], Tdat['T2'], color='r')
								plt.xlabel('umax/umax0')
								plt.ylabel('[T]ss')
								plt.title('param set #: ' + str(setnum))

								plt.subplot(1,2,2)
								pdf.savefig(fig)
								# plt.savefig(outdir + str(setnum))
								plt.close()

							t1 = time.time() - t0
							timeleft = (t1*N - t1*toploopcount)/60

pdf.close()




					