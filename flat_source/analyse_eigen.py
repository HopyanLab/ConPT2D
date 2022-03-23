#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import argparse
import pickle
from scipy.optimize import curve_fit

################################################################################
#===============================================================================
# analyse_eigen.py
#===============================================================================
################################################################################

def F(x, x0, a, b):
	return np.where( (x-x0) <= b/a, 0, a-b/(x-x0) )

def plot_data (data_dir = Path('../flat_eigenvalues'),
			   output_file = Path('../plots/flat_eigen_plot.svg'),
			   return_ax = False):
	fig, ax = plt.subplots(1)
	run_values = np.genfromtxt(data_dir/'run_values.csv',
								delimiter=',',dtype = int)
	r_values = np.genfromtxt(data_dir/'r_values.csv',
								delimiter=',', dtype = float)
	p_values = np.genfromtxt(data_dir/'p_values.csv',
								delimiter=',', dtype = float)
	if (data_dir / 'data.pkl').exists():
		with open(data_dir / 'data.pkl', 'rb') as filestream:
			data = pickle.load(filestream)
	else:
		data = np.zeros((len(run_values), len(r_values),
						 len(p_values), 256), dtype = float)
		for run_index, run_number in enumerate(run_values):
			run_dir = (data_dir/f'run_{run_number:d}')
			if not run_dir.is_dir():
				continue
			for r_index,r_param in enumerate(r_values):
				for p_index,p_param in enumerate(p_values):
					csvfile = run_dir / ('eigenvalues' + \
										 '_r_{0:1.1f}'.format(r_param) + \
										 '_p0_{0:1.3f}.csv'.format(p_param))
					data[run_index,r_index,p_index,:] = np.genfromtxt(csvfile,
																delimiter=',')
		with open(data_dir / 'data.pkl', 'wb') as filestream:
			pickle.dump(data, filestream)
	floppy_modes = np.count_nonzero(data<=0,axis=3)
	floppy_modes = floppy_modes - 2
	floppy_modes[floppy_modes<0] = 0
	floppy_modes = floppy_modes/254
	floppy_modes = np.median(floppy_modes, axis=0)
	markers = ['*', 'p', 's', 'v', 'o']
	colors = ['purple', 'blue', 'green', 'red', 'black']
#	for r_index,r_param in enumerate(r_values):
#		ax.plot(p_values,floppy_modes[r_index,:],
#				color = colors[r_index],
#				marker = markers[r_index])
	floppy_modes = np.mean(floppy_modes, axis=0)
	params, covar = curve_fit(F,
							  p_values,
							  floppy_modes,
							  p0 = [-0.2,0.25,1])
#	print(params)
	print(params[2]/params[1]+params[0])
	print(np.sqrt((params[2]/params[1]**2)**2*covar[1,1] + \
				  (1/params[1])**2*covar[2,2] + covar[0,0]))
	print(params[1])
	ax.grid(True)
	x_values = np.linspace(3.,4.,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,0.3,11)
	ax.plot(np.ones_like(y_values)*3.813, y_values,
			color = 'black', linestyle = '--', marker = '')
	x = np.linspace(3.6,3.9,300)
	ax.plot(x, F(x,*params),
			color = 'black',
			linestyle = '-',
			marker = '')
	ax.plot(p_values,floppy_modes,
			color = 'black',
			linestyle = '',
			marker = 'o')
	ax.set_ylim((-0.01,0.26))
	ax.set_xlim((3.715,3.885))
	ax.set_xlabel('$p_0$')
	ax.set_ylabel('$N(\lambda=0)$')
#	ax.legend()
	plt.savefig(output_file)
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(output_file.with_suffix('.pgf'))
	if return_ax:
		return ax

################################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
							description = '')
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = '../plots/flat_eigen_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('datadir',
						nargs = '?',
						default = '../flat_eigenvalues.paper',
						type = str,
						help = 'directory of eigenvalues data')
	args = parser.parse_args()
	plot_data(data_dir = Path(args.datadir),
			  output_file = Path(args.outfile) )

################################################################################
# EOF
