import datatable as dt
import matplotlib.pyplot as plt

ptable = dt.fread('/Users/taryn/GoogleDrive/code/wcEcoli/prototypes/toxin_antitoxin/plots/multiD/param_sets.csv')
features = ['sigma', 'alpha', 'kh', 'lambda_A', 'lambda_T', 'kt', 'n']

ptable = ptable.to_numpy()[:,1:]
N=5
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



plt.boxplot(ptable,0, '')
x = 1
for f in features:
	plt.scatter(x, min(range_dict[f]), color='r')
	plt.scatter(x, max(range_dict[f]), color='r')
	x += 1

plt.xticks(range(1,len(features)+1), labels=features)