import datatable as dt
from sklearn.decomposition import PCA

ptable = dt.fread('/Users/taryn/GoogleDrive/code/wcEcoli/prototypes/toxin_antitoxin/plots/multiD/param_sets.csv')
features = list(ptable.keys())[1:]

ptable = ptable.to_numpy()[:,1:]

ptable_reg = np.subtract(ptable, np.tile(np.mean(ptable, axis=0), (ptable.shape[0],1)))/np.tile(np.std(ptable, axis=0), (ptable.shape[0],1))

n_comp = 4
pca = PCA(n_components=n_comp)

W = pca.fit_transform(ptable_reg)

plt.figure(figsize=(19,3))
plt.subplot(151)
plt.bar(np.arange(1,n_comp+1), pca.explained_variance_ratio_)
plt.xticks(np.arange(1, n_comp+1))
plt.xlabel('components')
plt.ylabel('variance explained')

for s in range(0,n_comp):
	plt.subplot(1,5,s+2)
	plt.bar(np.arange(0,7), pca.components_[s, :])
	plt.xticks(np.arange(0,7),labels=features, rotation=45, fontsize=6)
	plt.axhline(0, color='k')
	plt.title('comp ' + str(s+1))

plt.tight_layout()