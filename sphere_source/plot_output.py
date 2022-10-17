#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import argparse
from pathlib import Path
from scipy import optimize
from scipy.special import erf

################################################################################
#===============================================================================
# plot_spherical_output.py
#===============================================================================
################################################################################

#3.80573442,0.0796802,0.95420221,3.66025486,0.01169476

################################################################################

def smooth_runaway (x):
	return np.where(x>0, np.exp(-1/x), 0)

def smooth_step (x):
	return smooth_runaway(x)/(smooth_runaway(x)+smooth_runaway(1-x))

################################################################################

def F (x,x0, alpha, m, b, c):
	return smooth_step((x0-x)/alpha) * G(x, m, b, c)

def G (x, m, b, c):
	return -m*x+b + np.sqrt( (-m*x+b)**2 + c**2 )

def F_flat (x,x0, alpha, m, b, c):
	return ((np.tanh((x0-x)/alpha) + 1)/2)**2 * G(x, m, b, c)

def F_shift (x, x0):
	return F_flat(x+x0,3.80573442,0.0796802,0.95420221,3.66025486,0.01169476)

################################################################################

#def CDF (x,mu,sigma):
#	return 0.5*(1+erf((x-mu)/sigma/np.sqrt(2)))

#def F (x, m, b, mu, sigma):
#	(-m*x-b)*(1-CDF(x,mu,sigma))

################################################################################

n_params = np.array([24, 32, 48, 64, 128])
markers = np.array(['o', 'v', 's', 'P', 'X'])
colors = np.array(['purple', 'blue', 'green',
					'darkorange', 'red'])
fit_function = F
p0 = [3.8, 0.1, 1., 3.6, 0.1]
#fit_function = F_shift
#p0 = [0.04]

def plot_output (result_file = Path('../sphere_output/results.csv'),
				 output_file = Path('../plots/sphere_output_plot.svg'),
				 fine_scale = False, draw_box = False, return_ax = False):
	fig, ax = plt.subplots(1)
	output_file.parent.mkdir(parents = True, exist_ok = True)
	energy_data = np.loadtxt(result_file, comments = '#', delimiter = ',')
	n_param = energy_data[:,0]
	p0_param = energy_data[:,1]
	median = energy_data[:,2]
	mean = energy_data[:,3]
	error = energy_data[:,4]
	x_values = np.linspace(3.,4.,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,0.5,11)
	ax.plot(np.ones_like(y_values)*3.813, y_values,
			color = 'black', linestyle = '--', marker = '')
#	p0_values = p0_fits = np.array([3.7451,3.7711,3.7898,3.8016,3.8098])
#	for index, current_n in enumerate(n_params):
#		ax.plot(np.ones_like(y_values)*p0_values[index], y_values,
#				color = colors[index], linestyle = '--', marker = '')
	x = np.linspace(3.58,3.92,200)
	curve_handles = np.empty(len(n_params), dtype=object)
	black_handles = np.empty(len(n_params), dtype=object)
	color_handles = np.empty(len(n_params), dtype=object)
	for index, current_n in enumerate(n_params):
		fit_mask = (n_param == current_n)
		plot_mask = (n_param == current_n) & ((p0_param*1000)%10 == 0)
		params, covar = optimize.curve_fit( fit_function,
										p0_param[fit_mask],
										mean[fit_mask],
										p0 = p0,
										sigma = error[fit_mask],
										bounds = (0,24) )
	#	print(f'{params[0]:1.4f} +/-{covar[0,0]:1.4f}')
		y = fit_function(x,*params)
		curve_handles[index], = ax.plot(x,y,
									color = colors[index],
									linestyle = '-')
		black_handles[index], = ax.plot(p0_param[plot_mask],
									mean[plot_mask],
									color = 'black', linestyle = '',
									marker = markers[index], markersize = 6.)
		color_handles[index], = ax.plot(p0_param[plot_mask],
									mean[plot_mask],
									color = colors[index], linestyle = '',
									marker = markers[index], markersize = 3.,
									label = 'N = {0:d}'.format(current_n))
	ax.plot(x, F_shift(x,0),
			 color = 'black', linestyle = '-', marker = '',
			 label = r'N $\rightarrow \infty$ (flat)')
	ax.grid(True)
	if(not fine_scale):
		ax.set_ylim((-0.01,0.355))
		ax.set_xlim((3.645,3.905))
		if(draw_box):
			xx = np.linspace(3.715, 3.825, 200)
			yy = np.linspace(-0.005, 0.14, 200)
			ax.plot(xx, np.ones_like(xx)*yy[0],
					 color = 'gray', linestyle = 'dotted', marker = '')
			ax.plot(xx, np.ones_like(xx)*yy[-1],
					 color = 'gray', linestyle = 'dotted', marker = '')
			ax.plot(np.ones_like(yy)*xx[0], yy,
					 color = 'gray', linestyle = 'dotted', marker = '')
			ax.plot(np.ones_like(yy)*xx[-1], yy,
					 color = 'gray', linestyle = 'dotted', marker = '')
	else:
		ax.set_ylim((-0.005,0.14))
		ax.set_xlim((3.715,3.825))
#	ax.set_xticks((p0_param[n_param == 24])[::5])
#	ax.set_title('Energy Barrier versus Shape Index')
	ax.set_xlabel(r'$p_0$')
	ax.set_ylabel(r'$\Delta \varepsilon$')
	handles, labels = ax.get_legend_handles_labels()
	combined_handles = handles
	combined_handles[:len(n_params)] = zip(curve_handles,
											black_handles,
											color_handles)
	ax.legend(combined_handles,labels,
		  loc = 'best', fancybox = True, framealpha = 1.)
#	plt.tight_layout()
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
						required = False,
						default = '../plots/sphere_output_plot.svg',
						type = str,
						help = 'file to put plot in')
	parser.add_argument('-b', '--box', dest='draw_box',
						action='store_const',
						const=True, default=False,
						help = 'Put a dotted box around finer scale.')
	parser.add_argument('-f', '--fine', dest='fine_scale',
						action='store_const',
						const=True, default=False,
						help = 'Plot finer scale.')
	parser.add_argument('result_file',
						nargs = '?',
						default = '../sphere_output/results.csv',
						type = str,
						help = 'file with simulation result data')
	args = parser.parse_args()
	plot_output(result_file = Path(args.result_file),
				output_file = Path(args.outfile),
				fine_scale = args.fine_scale,
				draw_box = args.draw_box )

################################################################################
# EOF
