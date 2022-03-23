#!/usr/bin/env /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import curve_fit
from pathlib import Path
import argparse
import tikzplotlib

def k_gamma (x, k):
	return k**k*x**(k-1)*np.exp(-k*x)/gamma(k)

def k_theta_gamma (x, k, theta):
	return x**(k-1)*np.exp(-x/theta)/gamma(k)/theta**k

def plot_dist (out_dir = '../flat_output/',
			   plot_file = '../plots/flat_gamma_plot.svg',
			   return_ax = False):
	fig, ax = plt.subplots(1)
	x = np.linspace(0.001,4.8,1000)
	k = 2.1
	theta = 1.0
	r_params = np.array([8.0, 2.0, 1.0, 0.5, 0.1])
	markers = np.array(['*', 'p', 's', 'v', 'o'])
	colors = np.array(['purple', 'blue', 'green', 'red', 'black'])
	num_bins = 200
	full_bins = np.zeros((len(r_params),num_bins))
	full_data = np.zeros((len(r_params),num_bins))
	for csvfile in out_dir.glob('r_[0-9]*.[0-9]*_p0_[0-9]*.[0-9]*.csv'):
		if float(csvfile.stem.split('_')[-1]) < 3.701:
			file_r = float(csvfile.stem.split('_')[-3])
			color = colors[file_r == r_params][0]
			marker = markers[file_r == r_params][0]
			data = np.genfromtxt(csvfile)
			avg = np.mean(data)
			binned_data, bins = np.histogram(data/avg,
									bins = num_bins, range = (0,6))
			bin_size = bins[1]-bins[0]
			bins = (bins[:-1]+bins[1:])/2
			binned_data = binned_data / (bin_size * np.sum(binned_data))
			full_bins[file_r == r_params] = bins
			full_data[file_r == r_params] += binned_data

	for index in range(len(r_params)):
		bin_size = full_bins[index,1]-full_bins[index,0]
		full_data[index] = full_data[index] / \
								(bin_size * np.sum(full_data[index]))
		ax.plot(full_bins[index], full_data[index],
				color = colors[index],
				marker = markers[index],
				markersize = 2.0,
				label = 'r = {0:1.1f}'.format(r_params[index]),
				linestyle = '')

#	fit_params, pcov = curve_fit(k_theta_gamma, full_bins.flatten(),
#									full_data.flatten(), p0=[k,theta])
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
	ax.set_ylim([7e-4,1.2e0])
	ax.legend(loc = 'best', fancybox = True, framealpha = 1.)
	plt.savefig(plot_file)
#	tikzplotlib.Flavors.latex.preamble()
#	tikzplotlib.save(plot_file.with_suffix('.tex'))
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(plot_file.with_suffix('.pgf'))
	if return_ax:
		return ax
#	plt.show()

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Plot results of simulation.')
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = ['../plots/flat_gamma_plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('data_path',
						nargs = '?',
						default = '../flat_output/',
						type = str,
						help = 'path to simulation result data')
	args = parser.parse_args()
	plot_dist(out_dir = Path(args.data_path),
			  plot_file = Path(args.outfile[0]) )

################################################################################
# EOF
