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

def phi (x, l):
	if l == 0:
		return np.exp(-x)
	else:
		return np.where(l*x<1, (1-l*x)**(1/l), 0)

def f (x, a, b):
	return phi(phi(x,b),a)

def g (x, x0, A, a, b):
	return A*f((x-x0)/a/b,a,b)

def plot_eigen (data_dir = Path('../bounded_eigenvalues'),
				output_file = Path('../plots/bounded_eigen_plot.svg'),
				return_ax = False):
	fig, ax = plt.subplots(1)
	run_values = np.genfromtxt(data_dir/'run_values.csv',
								delimiter=',',dtype = int)
	r_values = np.genfromtxt(data_dir/'r_values.csv',
								delimiter=',', dtype = float)
	r_values = np.atleast_1d(r_values)
	p_values = np.genfromtxt(data_dir/'p_values.csv',
								delimiter=',', dtype = float)
	bad_runs = np.empty(0, dtype = int)
	if (data_dir / 'data.pkl').exists():
		with open(data_dir / 'data.pkl', 'rb') as filestream:
			data = pickle.load(filestream)
	else:
		data = np.zeros((len(run_values), len(r_values),
						 len(p_values), 256), dtype = float)
		for run_index, run_number in enumerate(run_values):
			run_bad = False
			run_dir = (data_dir/f'run_{run_number:02d}')
			if not run_dir.is_dir():
				bad_runs = np.append(bad_runs, run_index)
				continue
			for r_index,r_param in enumerate(r_values):
				if run_bad:
					continue
				for p_index,p_param in enumerate(p_values):
					if run_bad:
						continue
					try:
						csvfile = run_dir / ('eigenvalues' + \
										 '_r_{0:1.1f}'.format(r_param) + \
										 '_p0_{0:1.3f}.csv'.format(p_param))
						temp_data = np.genfromtxt(csvfile, delimiter=',')
						if len(temp_data) > data.shape[-1]:
							data = np.pad(data,((0,0),(0,0),(0,0),
										(0,len(temp_data)-data.shape[-1])),
										mode = 'constant',
										constant_values = 1000.)
						if len(temp_data) < data.shape[-1]:
							temp_data = np.pad(temp_data,
									((0,data.shape[-1]-len(temp_data))),
											mode = 'constant',
										constant_values = 1000.)
						data[run_index,r_index,p_index,:] = temp_data
					except:
						run_bad = True
						bad_runs = np.append(bad_runs, run_index)
		data = np.delete(data, bad_runs, axis=0)
		with open(data_dir / 'data.pkl', 'wb') as filestream:
			pickle.dump(data, filestream)
#	floppy_modes = np.zeros_like(data[:,:,:,0])
#	for run_index in np.arange(data.shape[0]):
#		for r_index in np.arange(data.shape[1]):
#			for p_index in np.arange(data.shape[2]):
#				floppy_modes[run_index, r_index, p_index] = np.count_nonzero(
#					data[run_index, r_index, p_index,:] <= \
#						1e-8 * np.mean(data[run_index, r_index, p_index,:]))
	total_modes = data.shape[3]
	floppy_modes = np.count_nonzero(data<=1e-8 * \
					np.mean(data,axis=3)[:,:,:,np.newaxis],axis=3)
	floppy_modes[floppy_modes<0] = 0
	floppy_modes = floppy_modes/total_modes
	floppy_modes = np.median(floppy_modes, axis=0)
	markers = ['*', 'p', 's', 'v', 'o']
	colors = ['purple', 'blue', 'green', 'red', 'black']
#	for r_index,r_param in enumerate(r_values):
#		ax.plot(p_values,floppy_modes[r_index,:],
#				color = colors[r_index],
#				marker = markers[r_index])
	floppy_modes = np.mean(floppy_modes, axis=0)
	params, covar = curve_fit(g,
							  p_values,
							  floppy_modes,
							  p0 = [3.85, 0.2,0.1,0.1])
	print(params)
#	np.savetxt('bounded_eigen.csv',np.vstack([p_values,floppy_modes]).T,
#				delimiter = ',')
	ax.grid(True)
	x_values = np.linspace(3.,5.,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '')
	y_values = np.linspace(0.,0.5,11)
	x = np.linspace(3.80,4.1,1000)
	ax.plot(x, g(x,*params),
			color = 'black',
			linestyle = '-',
			marker = '')
	ax.plot(p_values,floppy_modes,
			color = 'black',
			linestyle = '',
			marker = 'o')
	ax.set_ylim((-0.01,0.26))
	ax.set_xlim((3.795,4.105))
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
						default = '../plots/bounded_eigen_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('datadir',
						nargs = '?',
						default = '../bounded_eigenvalues',
						type = str,
						help = 'directory of eigenvalues data')
	args = parser.parse_args()
	plot_eigen(data_dir = Path(args.datadir),
			   output_file = Path(args.outfile) )

################################################################################
# EOF
