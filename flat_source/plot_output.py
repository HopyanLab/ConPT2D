#!/usr/bin/env /usr/bin/python3
import numpy as np
import argparse
from matplotlib import pyplot as plt
from pathlib import Path
from scipy import optimize

################################################################################
#===============================================================================
# plot_output.py
#===============================================================================
################################################################################

def smooth_runaway (x):
	return np.where(x>0, np.exp(-1/x), 0)

def smooth_step (x):
	return smooth_runaway(x)/(smooth_runaway(x)+smooth_runaway(1-x))

def smooth_k_step (x,k): # k should be >=1
	return np.where(x>0, np.where(x<1,
		0.5*(np.tanh(k*(2*x-1)/2/np.sqrt(x*(1-x)))+1), 1), 0)

################################################################################

def F (x, x0, alpha, m, b, c):
#	return smooth_step((x0-x)/alpha) * G(x, m, b, c)
	return ((np.tanh((x0-x)/alpha) + 1)/2)**2 * G(x, m, b, c)

def G (x, m, b, c):
	return -m*x+b + np.sqrt( (-m*x+b)**2 + c**2 )

################################################################################

def plot_output (result_file = Path('../flat_output/results.csv'),
				 output_file = Path('../plots/flat_output_plot.svg'),
				 return_ax = False):
	fig, ax = plt.subplots(1)
	x_values = np.linspace(2.8,4.2,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,1.7,11)
	ax.plot(np.ones_like(y_values)*3.813, y_values,
			color = 'black', linestyle = '--', marker = '')
		#	label = '3.813')
	x = np.linspace(3,4,1000)
	energy_data = np.loadtxt(result_file, comments = '#', delimiter = ',')
	r_param = energy_data[:,0]
	p0_param = energy_data[:,1]
	energy = energy_data[:,2]
	error = energy_data[:,3]
	r_params = np.unique(r_param)[::-1]
	markers = ['*', 'p', 's', 'v', 'o']
	colors = ['purple', 'blue', 'green', 'red', 'black']
	p0 = [3.8, 0.08, 1.0, 3.6, 0.1]
	test = np.array([])
	for index, current_r in enumerate(r_params):
		ax.plot(p0_param[r_param == current_r],
				energy[r_param == current_r],
				color = colors[index], linestyle = '',
				marker = markers[index], markersize = '4.',
				label = f'r = {current_r:1.1f}')
#		ax.errorbar(p0_param[r_param == current_r],
#					energy[r_param == current_r],
#					error[r_param == current_r],
#					color = colors[index], linestyle = '',
#					marker = markers[index], markersize = '4.',
#					label = f'r = {current_r:1.1f}')
		params, covar = optimize.curve_fit(F, p0_param[r_param == current_r],
											  energy[r_param == current_r],
											  p0 = p0,
											  bounds = (0,24))
	#	if current_r == 0.5:
	#		print(params)
		p0 = params
		test = np.append(test,p0[0])
		y = F (x,*params)
		ax.plot(x,y, color = colors[index], linestyle = '-')
	print(np.mean(test))
	print(np.std(test))
	ax.grid(True)
	ax.set_xlim((2.95,4.05))
	ax.set_ylim((-0.05,1.65))
	ax.set_xticks(np.linspace(3.0,4.0,11))
#	ax.set_title('Energy Barrier as a Function of Shape Index')
	ax.set_xlabel(r'$p_0$')
	ax.set_ylabel(r'$\Delta \varepsilon$')
	ax.legend()
	plt.savefig(output_file)
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(output_file.with_suffix('.pgf'))
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
						default = ['../plots/flat_output_plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('result_file',
						nargs = '?',
						default = '../flat_output/results.csv',
						type = str,
						help = 'file with simulation result data')
	args = parser.parse_args()
	plot_output( result_file = Path(args.result_file),
				 output_file = Path(args.outfile[0]) )

################################################################################
# EOF
