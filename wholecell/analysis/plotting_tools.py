"""
Reusable plotting functions and tools
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import numpy as np

DEFAULT_MATPLOTLIB_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

COLORS_LARGE = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
		"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
		"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
		"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
		"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
		"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
		"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
		"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

		"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
		"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
		"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
		"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
		"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
		"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
		"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
		"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]

COLORS_SMALL = ["#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#000000", "#007FFF",
		"#236B8E", "#70DB93", "#B5A642", "#5F9F9F", "#B87333", "#2F4F2F", "#9932CD", "#871F78", "#855E42",
		"#545454", "#8E2323", "#238E23", "#CD7F32", "#527F76",
		"#9F9F5F", "#8E236B", "#2F2F4F", "#CFB53B", "#FF7F00", "#DB70DB",
		"#5959AB", "#8C1717", "#238E68", "#6B4226", "#8E6B23", "#00FF7F",
		"#38B0DE", "#DB9370", "#5C4033", "#4F2F4F", "#CC3299", "#99CC32"]

CMAP_COLORS_255 = [
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

with plt.style.context('seaborn-v0_8-colorblind'):
	COLORS_COLORBLIND = plt.rcParams['axes.prop_cycle'].by_key()['color']

def plotSplom(arrayOfdataArrays, nameArray="", stdArrays=None, labels=None, fig=None, plotCorrCoef=True, formatString='o'):
	"""
	Plot a scatterplot matrix (Splom) of data contained in arrayOfdataArrays,
	with labels in the same order held within nameArray.
	"""

	if len(arrayOfdataArrays) != len(nameArray):
		raise IndexError("Your array of data arrays and the array of names must be the same length.")

	if stdArrays is None:
		stdArrays = [None]*len(arrayOfdataArrays)

	if len(stdArrays) != len(arrayOfdataArrays):
		raise IndexError("If you provide an array of standard deviations, there must be one entry per input data array. Entries can be None.")

	if fig is None:
		fig = plt.figure()

	num_entries = len(arrayOfdataArrays)
	plottingIndex = 1
	for rowNum in range(1,num_entries+1):
		for colNum in range(1,num_entries+1):
			if colNum < plottingIndex:
				continue
			plt.subplot(num_entries,num_entries,num_entries*(rowNum-1)+(colNum))

			plt.errorbar(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1], xerr=stdArrays[colNum-1], yerr=stdArrays[rowNum-1], fmt=formatString)

			if nameArray != "":
				plt.xlabel(nameArray[colNum-1])
				plt.ylabel(nameArray[rowNum-1])

			if plotCorrCoef:
				corr_coef, pValue = stats.pearsonr(arrayOfdataArrays[colNum-1], arrayOfdataArrays[rowNum-1])
				plt.title("R = %.4f" % (corr_coef,))
		plottingIndex += 1

	return fig

def labeled_indexable_hist(obj, ax, data, gen_data, gen_start, gen_end, colors, xlabel, bin_width=1., xlim=None, sf=1, font_size=9):
	"""
	Creates a histogram of (subset of) data, with label for mean and standard
	deviation of data for each variant

	Args:
		obj: specify the Plot object
		ax: Axes object
		data: data to plot
		gen_data: generation index corresponding to each data point
		gen_start: index of generation to start from (inclusive)
		gen_end: index of generation to end at (exclusive)
		colors: list of colors to use for each variant
		xlabel: x-axis label for plot
		bin_width: used in conjunction with xlim to determine number of bins
		xlim: specify x-axis plotting region
		sf: scale factor
		font_size: font size for labeling axes

	Returns:
		histogram of data, colored by variant, for data corresponding to
			generation indexes in [gen_start:gen_end]

	"""

	if xlim:
		bins = np.histogram(range(xlim[0],xlim[1]+1),
							bins=int(np.ceil((xlim[1]-xlim[0])/bin_width)))[1]

	for variant, variant_data in data.items():
		variant_gen_data = gen_data[variant]
		variant_data = variant_data[(variant_gen_data >= gen_start) &
									(variant_gen_data < gen_end)]

		if not variant_data.any():
			continue

		color = colors[variant % len(colors)]
		if not xlim:
			bins = max(1, int(np.ceil((variant_data.max() - variant_data.min())
									  / bin_width)))
		mean = variant_data.mean()
		std = variant_data.std()
		ax.hist(variant_data, bins, color=color, alpha=0.5,
			label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
		ax.axvline(mean, color=color, linestyle='--', linewidth=1)

	if xlim:
		ax.set_xlim(xlim)
	obj.remove_border(ax)
	ax.set_xlabel(xlabel, fontsize=font_size)
	ax.tick_params(labelsize=font_size)
	ax.legend()

def labeled_indexable_scatter(obj, ax, xdata, ydata, gen_data, gen_start,
							  gen_end, colors, xlabel, ylabel, xlim=None,
							  ylim=None, sf=1, font_size=9):
	"""
	Creates a scatterplot of (subset of) data, with label for mean and
	standard deviation of data for each variant

	Args:
		obj: specify the Plot object
		ax: Axes object
		xdata: data to plot on x axes
		ydata: data to plot on y axes
		gen_data: generation index corresponding to each data point
		gen_start: index of generation to start from (inclusive)
		gen_end: index of generation to end at (exclusive)
		colors: list of colors to use for each variant
		xlabel: x-axis label for plot
		ylabel: y-axis label for plot
		xlim: specify x-axis plotting region
		ylim: specify y-axis plotting region
		sf: scale factor
		font_size: font size for labeling axes

	Returns:
		scatterplot of data, colored by variant, for data corresponding
	 		to generation indexes in [gen_start:gen_end]

	"""

	for variant, variant_data in xdata.items():
		variant_gen_data = gen_data[variant]
		variant_xdata = variant_data[(variant_gen_data >= gen_start) &
									 (variant_gen_data < gen_end)]
		variant_ydata = ydata[variant][(variant_gen_data >= gen_start) &
									   (variant_gen_data < gen_end)]

		if not (variant_xdata.any() or variant_ydata.any()):
			continue

		color = colors[variant % len(colors)]
		mean = variant_ydata.mean()
		std = variant_ydata.std()
		mean_new_gene = variant_xdata.mean()
		ax.scatter(variant_xdata, variant_ydata, color=color, alpha=0.5,
			label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
		ax.scatter(mean_new_gene,mean,color=color,alpha=0.5,marker='x')

	if xlim:
		ax.set_xlim(xlim)
	if ylim:
		ax.set_ylim(ylim)
	obj.remove_border(ax)
	ax.set_xlabel(xlabel, fontsize=font_size)
	ax.set_ylabel(ylabel, fontsize=font_size)
	ax.tick_params(labelsize=font_size)
	ax.legend()

def  heatmap(obj, ax, mask, data, completion_data, xticklabels, yticklabels,
			xlabel="", ylabel="", title="", box_text_size ="medium",
			ax_font_size=9, title_font_size=9, percent_completion_threshold = 0.88):
	"""
	Args:
		obj: specify the Plot object
		ax: Axes object
		mask: Only plot values where mask is true, must match dimensions of
		data
		data: 2-dimensional numpy array of data to plot
		completion_data: Percent of seeds that successfully completed all
		generations that contributed to this value, must match dimensions of
		data
		xticklabels: tick values for x-axis
		yticklabels: tick values for y-axis
		xlabel: x-axis label for plot
		ylabel: y-axis label for plot
		title: plot title
		box_text_size: size of text value to be printed in box
		ax_font_size: font size for labeling axes
		title_font_size: font size for title
		percent_completion_threshold: If the percent completion for this
		parameter combination is lower than the threshold, the number in the
		box will be colored red. If the threshold is 0, no numbers will be
		colored red.

	Returns:
		heatmap of data, where squares are colored by value and numbers
		are colored by the percent of simulations that successfuly completed,
		for data corresponding to different parameter values in 2 dimensions

	"""

	assert(mask.shape == completion_data.shape == data.shape)

	grid_colors = [(255 / 255, 255 / 255, 255 / 255),
				   (22 / 255, 110 / 255, 164 / 255)]
	cmap_name = 'blue_cmap'
	blue_cmap = LinearSegmentedColormap.from_list(cmap_name, grid_colors, N=100)

	im = ax.imshow(data, cmap=blue_cmap)
	ax.set_xticks(np.arange(len(xticklabels)))
	ax.set_xticklabels(xticklabels)
	ax.set_yticks(np.arange(len(yticklabels)))
	ax.set_yticklabels(yticklabels)
	plt.setp(
		ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
	for i in range(len(yticklabels)):
		for j in range(len(xticklabels)):
			if mask[i,j]:
				col = "k"
				if completion_data[i,j] == 0 and data[i,j] == -1:
					continue
				if completion_data[i,j] < percent_completion_threshold:
					col = "r"
				text = ax.text(
					j, i, data[i, j], ha="center", va="center", color=col,
					fontsize=box_text_size)
	ax.set_xlabel(xlabel, fontsize=ax_font_size)
	ax.set_ylabel(ylabel, fontsize=ax_font_size)
	ax.set_title(title, fontsize=title_font_size)
