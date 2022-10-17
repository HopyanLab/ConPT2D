#!/usr/bin/env /usr/bin/python3
import numpy as np
import argparse
from matplotlib import pyplot as plt
from pathlib import Path
from scipy import optimize
import tikzplotlib as tpl

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

def interpolate (x, x0, alpha):
	return ((np.tanh((x0-x)/alpha) + 1)/2)**2

################################################################################

def F (x, x0, alpha, m, b, c):
#	return smooth_step((x0-x)/alpha) * G(x, m, b, c)
	return ((np.tanh((x0-x)/alpha) + 1)/2)**2 * G(x, m, b, c)

def G (x, m, b, c):
	return -m*x+b + np.sqrt( (-m*x+b)**2 + c**2 )

def F_flat (x, x0):
	return F(x+x0, 3.80573442, 0.0796802, 0.95420221, 3.66025486, 0.01169476)
#	return F(x+x0, 3.92888633, 0.04715518, 1.20032791, 4.53721188, 0.01816722)

def F_bound (x, x0):
#	return F(x+x0, 3.91715964,0.06754252,1.09835944,4.20105076,0.02856569)
	return F(x+x0, 3.91653039, 0.06823231, 1.08406733, 4.1480409, 0.02766134)

################################################################################

def plot_output (result_file = Path('../hyperbolic_output/results.csv'),
				 bounded_result_file = Path('../bounded_output/results.csv'),
				 output_file = Path('../plots/hyperbolic_output_plot.svg'),
				 return_ax = False):
	fig, ax = plt.subplots(1)
	x_values = np.linspace(2.8,4.2,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,1.7,11)
	ax.plot(np.ones_like(y_values)*3.813, y_values,
			color = 'black', linestyle = '--', marker = '',
			zorder = 1)
		#	label = '3.813')
	x = np.linspace(3,4,1000)
	energy_data = np.loadtxt(result_file, comments = '#', delimiter = ',')
	bounded_data = np.loadtxt(bounded_result_file, comments = '#',
													delimiter = ',')
	bounded_r_param = bounded_data[:,0]
	bounded_p0_param = bounded_data[:,1]
	bounded_energy = bounded_data[:,2]
	bounded_p0_param = bounded_p0_param[bounded_r_param == 0.5]
	bounded_energy = bounded_energy[bounded_r_param == 0.5]
	area_param = energy_data[:,0]
	p0_param = energy_data[:,1]
	energy = energy_data[:,2]
	error = energy_data[:,3]
	area_params = np.unique(area_param)[::-1]
	#markers = ['*', 'P', 'p', 's', 'v', 'o']
	#colors = ['purple', 'blue', 'green', 'orange', 'red', 'black']
	markers = np.array(['h', 'H', 'p', 'd', '^'])
	colors = np.array(['tab:purple', 'tab:blue', 'tab:green',
							'tab:orange', 'tab:red'])
	p0 = [3.8, 0.08, 1.0, 1.0, 0.1]
	test = np.array([])
	curve_handles = np.empty(len(area_params)+1, dtype=object)
	black_handles = np.empty(len(area_params)+1, dtype=object)
	color_handles = np.empty(len(area_params)+1, dtype=object)
	for index, current_area in enumerate(area_params):
		number = int(np.round((current_area-0.1)/np.pi))
		if number == 0:
			number = int(np.round(np.pi/current_area))
			string = '$\pi$/'+str(number)
		elif number == 1:
			string = '$\pi$'
		else:
			string = str(number)+'$\pi$'
		black_handles[index], = ax.plot(p0_param[area_param == current_area],
									energy[area_param == current_area],
									color = 'black',
									linestyle = '',
									marker = markers[index],
									markersize = '6.',
									zorder = 9-index)
		color_handles[index], = ax.plot(p0_param[area_param == current_area],
									energy[area_param == current_area],
									color = colors[index],
									linestyle = '',
									marker = markers[index],
									markersize = '3.',
									zorder = 10-index,
									label = f'A = {string:s}')
#		ax.errorbar(p0_param[r_param == current_r],
#					energy[r_param == current_r],
#					error[r_param == current_r],
#					color = colors[index], linestyle = '',
#					marker = markers[index], markersize = '4.',
#					label = f'r = {current_r:1.1f}')
		params, covar = optimize.curve_fit(F,
									p0_param[area_param == current_area],
									energy[area_param == current_area],
									p0 = p0,
									bounds = (0,24))
	#	if current_r == 0.5:
	#		print(params)
		p0 = params
		test = np.append(test,p0[0])
		y = F (x,*params)
		curve_handles[index], = ax.plot(x, y,
									color = colors[index],
									linestyle = '-',
									zorder = 8-index)
		shift,_ = optimize.curve_fit(F_bound,
									p0_param[area_param == current_area],
									energy[area_param == current_area],
									p0 = [0.0],
									bounds = (-1.,1.))
		print(shift[0])
	curve_handles[-1], = ax.plot(x, F_bound(x,0),
			color = 'gray', linestyle = '-', marker = '', zorder = 3)
	black_handles[-1], = ax.plot(bounded_p0_param, bounded_energy,
			color = 'black', linestyle = '', marker = '.',
			markersize = 6.,
			zorder = 3)
	color_handles[-1], = ax.plot(bounded_p0_param, bounded_energy,
			color = 'gray', linestyle = '', marker = '.',
			markersize = 3.,
			label = r'A $\rightarrow 0$ (flat)',
			zorder = 3)
	ax.plot(x, F_flat(x,0),
			color = 'black', linestyle = '-', marker = '',
			label = r'N $\rightarrow \infty$ (flat)',
			zorder = 2)
	ax.grid(True)
	ax.set_xticks(np.arange(3.7,4.0,0.05))
	ax.set_xlim((3.69,3.96))
	ax.set_ylim((-0.02,0.41))
#	ax.set_title('Energy Barrier as a Function of Shape Index')
	ax.set_xlabel(r'$p_0$')
	ax.set_ylabel(r'$\Delta \varepsilon$')
	handles, labels = ax.get_legend_handles_labels()
#	order = np.array([0,1,2,3,4,5,10,9,8,7,6])
#	handles = np.array(handles,dtype=object)[order]
#	labels = np.array(labels,dtype=object)[order]
	combined_handles = handles
	combined_handles[:len(area_params)+1] = zip(curve_handles,
												black_handles,
												color_handles)
	ax.legend(combined_handles,labels,
		  loc = 'best', fancybox = True, framealpha = 1.)
#	plt.tight_layout()
	plt.savefig(output_file)
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(output_file.with_suffix('.pgf'))
	tpl.save(output_file.with_suffix('.tex'))
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
						default = ['../plots/hyperbolic_output_plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('-b', '--boundfile',
						nargs = 1,
						default = ['../bounded_output/results.csv'],
						type = str,
						required = False,
						help = 'file with bounded flat sim data')
	parser.add_argument('result_file',
						nargs = '?',
						default = '../hyperbolic_output/results.csv',
						type = str,
						help = 'file with simulation result data')
	args = parser.parse_args()
	plot_output( result_file = Path(args.result_file),
				 bounded_result_file = Path(args.boundfile[0]),
				 output_file = Path(args.outfile[0]) )

################################################################################
# EOF
