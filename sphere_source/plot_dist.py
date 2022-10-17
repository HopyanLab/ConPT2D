#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import gamma
from scipy.optimize import curve_fit
from pathlib import Path
import argparse

def k_gamma (x, k):
	return k**k*x**(k-1)*np.exp(-k*x)/gamma(k)

def k_theta_gamma (x, k, theta):
	return x**(k-1)*np.exp(-x/theta)/gamma(k)/theta**k

################################################################################

def plot_distribution(sim_dir, plot_file, return_ax = False):
	fig, ax = plt.subplots(1)
	x = np.linspace(0.0001,5.5,1000)
	k = 2.1
	theta = 1.0
	n_params = np.array([24, 32, 48, 64, 128])
	p0_fits = np.array([3.7451,3.7711,3.7898,3.8016,3.8098])
	markers = np.array(['o', 'v', 's', 'P', 'X'])
	colors = np.array(['purple', 'blue', 'green',
						'orange', 'red'])
	num_bins = 200
	full_bins = np.zeros((len(n_params),num_bins))
	full_data = np.zeros((len(n_params),num_bins))
	for csvfile in sim_dir.glob('n_[0-9]*_p0_[0-9]*.[0-9]*.csv'):
		n_value = int(csvfile.stem.split('_')[1])
		p_value = float(csvfile.stem.split('_')[-1])
		indices = (n_value == n_params)
		if not np.any(indices):
			continue
		shift = 3.814 - p0_fits[indices]
		if p_value < 3.681 - shift:# and \
		#	 int(csvfile.stem.split('_')[1]) == 92:
		#	print(str(csvfile))

			if p_value > p0_fits[n_value == n_params]-0.1:
				continue
			color = colors[indices][0]
			marker = markers[indices][0]
			data = np.genfromtxt(csvfile)
			avg = np.mean(data)
			binned_data, bins = np.histogram(data/avg,
									bins = num_bins, range = (0,6))
			bin_size = bins[1]-bins[0]
			bins = (bins[:-1]+bins[1:])/2
			binned_data = binned_data / (bin_size * np.sum(binned_data))
			full_bins[indices] = bins
			full_data[indices] += binned_data
	for index in range(len(n_params)):
		if np.any(full_data[index] != 0.):
			bin_size = full_bins[index,1]-full_bins[index,0]
			full_data[index] = full_data[index] / \
								(bin_size * np.sum(full_data[index]))
			ax.plot(full_bins[index], full_data[index],
						color = colors[index],
						marker = markers[index],
						markersize = 3.0,
						label = 'n = {0:d}'.format(n_params[index]),
						linestyle = '')
#	fit_params, pcov = curve_fit(k_gamma, full_bins.flatten(),
#									full_data.flatten(), p0=[k, theta])
	fit_params, pcov = curve_fit(k_gamma, full_bins.flatten(),
									full_data.flatten(), p0=[k])
	print('fitk {0:f}'.format(fit_params[0]))
	print('vark {0:f}'.format(np.sqrt(pcov[0,0])))
#	print('fitθ {0:f}'.format(fit_params[1]))
#	print('varθ {0:f}'.format(np.sqrt(pcov[1,1])))
	ax.plot(x, k_gamma(x,*fit_params),
					color = 'black',
					linewidth = 3.0,
					linestyle = '-')
	ax.set_yscale('log')
	ax.grid(True)
#	ax.set_title('Distribution of Energy Barriers')
	ax.set_xlabel(
		r'$\Delta \varepsilon / \overline{ \Delta \varepsilon }$')
	ax.set_ylabel(
		r'$\rho ( \Delta \varepsilon / \overline{ \Delta \varepsilon } )$')
	ax.set_xlim([-0.1,5.5])
	ax.set_ylim([1.5e-3,1.2e0])
	ax.legend(loc = 'best', fancybox = True, framealpha = 1.)
	plt.savefig(plot_file)
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(plot_file.with_suffix('.pgf'))
	if return_ax:
		return ax
	#plt.show()

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Plot results of simulation.')
	parser.add_argument('-o', '--outfile',
						default = '../plots/sphere_gamma_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('sim_dir',
						nargs = '?',
						default = '../sphere_output',
						type = str,
						help = 'directory with simulation data')
	args = parser.parse_args()
	plot_distribution(sim_dir = Path(args.sim_dir),
					plot_file = Path(args.outfile) )

################################################################################
# EOF
