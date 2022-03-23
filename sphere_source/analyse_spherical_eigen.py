#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import argparse
import pickle
from scipy.optimize import curve_fit

################################################################################
#===============================================================================
# analyse_spherical_eigen.py
#===============================================================================
################################################################################

def F(x, x0, a, b):
	return np.where( (x-x0) <= b/a, 0, a-b/(x-x0) )

def plot_data (data_dir = Path('../sphere_eigenvalues'),
			   output_file = Path('../plots/flat_eigen_plot.svg'),
			   return_ax = False):
	fig, ax = plt.subplots(1)
	n_params = np.array([24, 32, 48, 64, 128])
	markers = np.array(['o', 'v', 's', 'P', 'X'])
	colors = np.array(['purple', 'blue', 'green',
						'orange', 'red'])
	dof = np.array([87, 119, 183, 247, 503])
	run_values = np.genfromtxt(data_dir/'run_values.csv',
								delimiter=',', dtype=int)
	n_values = np.genfromtxt(data_dir/'n_values.csv',
								delimiter=',', dtype=int)
	p_values = np.genfromtxt(data_dir/'p_values.csv',
								delimiter=',', dtype = float)
	if (data_dir / 'data.pkl').exists():
		with open(data_dir / 'data.pkl', 'rb') as filestream:
			data = pickle.load(filestream)
	else:
		data = np.array([], dtype = object)
		for n_index,n_param in enumerate(n_values):
			new_data = np.zeros((len(run_values), len(p_values),
								 dof[n_index]), dtype = float)
			data = np.append(data,[0])
			data [-1] = new_data
		for run_index, run_number in enumerate(run_values):
			run_dir = (data_dir/f'run_{run_number:d}')
			if not run_dir.is_dir():
				continue
			for n_index,n_param in enumerate(n_values):
				for p_index,p_param in enumerate(p_values):
					csvfile = run_dir / ('eigenvalues' + \
										 '_n_{0:d}'.format(n_param) + \
										 '_p0_{0:1.3f}.csv'.format(p_param))
					(data[n_index])[run_index,p_index,:] = np.genfromtxt(
																csvfile,
																delimiter=',')
		with open(data_dir / 'data.pkl', 'wb') as filestream:
			pickle.dump(data, filestream)
	for n_index,n_data in enumerate(data):
		data[n_index] = np.count_nonzero(n_data<=0,axis=2)
		data[n_index] = np.median(data[n_index], axis=0)
		data[n_index] = data[n_index]-2
		data[n_index][data[n_index]<0] = 0
#	factor = correction_factors(1,6,n_values)
	factor = np.ones_like(n_values)
	ax.grid(True)
	x_values = np.linspace(3.,4.,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,0.3,11)
	ax.plot(np.ones_like(y_values)*3.813, y_values,
			color = 'black', linestyle = '--', marker = '')
	x = np.linspace(3.6,3.9,300)
	for n_index,n_param in enumerate(n_values):
		params, covar = curve_fit(F,
								  p_values/factor[n_index],
								  data[n_index]/(dof[n_index]-2),
								  p0 = [-0.2,0.25,1])
		ax.plot(x, F(x,*params),
				color = colors[n_index],
				linestyle = '-',
				marker = '')
		print(params[2]/params[1]+params[0])
	#	print(np.amin(p_values[data[n_index]>0])/factor[n_index])
		ax.plot(p_values/factor[n_index], data[n_index]/(dof[n_index]-2),
				color = colors[n_index], linestyle = '',
				marker = markers[n_index],
				markersize = 4,
				label = 'n = {0:d}'.format(n_param))
	ax.plot(x, F(x,3.79501129,0.3143082,0.00594078),
			color = 'black',
			linestyle = '-',
			marker = '',
			label = r'n $\rightarrow \infty$ (flat)')
	ax.set_ylim((-0.01,0.26))
	ax.set_xlim((3.715,3.885))
	ax.set_xlabel('$p_0$')
	ax.set_ylabel('$N(\lambda=0)$')
	ax.legend(loc = 'best', fancybox = True, framealpha = 1.)
#	if output_file.suffix == '.pgf':
#		plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(output_file)
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(output_file.with_suffix('.pgf'))
	if return_ax:
		return ax

################################################################################

def correction_factors (A, sides, n_values):
	R = np.sqrt(n_values/4/np.pi)
	l_inf = np.sqrt(4*np.tan(np.pi/sides)*A/sides)
	def length (R, n, A):
		return 2*R * np.arccos(np.cos(np.pi/n)/np.cos(A/2/R**2/n-np.pi/n))
	factor = length(R,sides,A)/l_inf
	return factor

################################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
							description = '')
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = '../plots/sphere_eigen_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('datadir',
						nargs = '?',
						default = '../sphere_eigenvalues',
						type = str,
						help = 'directory of eigenvalues data')
	args = parser.parse_args()
	plot_data(data_dir = Path(args.datadir),
			  output_file = Path(args.outfile) )

################################################################################
# EOF
