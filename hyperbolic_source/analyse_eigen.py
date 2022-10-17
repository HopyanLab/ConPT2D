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

def F_flat (x):
	return F(x, 3.80055715, 2.78255922e-1, 3.17909568e-3)

def phi (x, l):
	if l == 0:
		return np.exp(-x)
	else:
		return np.where(l*x<1, (1-l*x)**(1/l), 0)

def f (x, a, b):
	return phi(phi(x,b),a)

def g (x, x0, A, a, b):
	return A*f((x-x0)/a/b,a,b)

def plot_eigen (data_dir = Path('../hyperbolic_eigenvalues'),
				output_file = Path('../plots/hyperbolic_eigen_plot.svg'),
				return_ax = False):
	fig, ax = plt.subplots(1)
	ax.grid(True)
	x_values = np.linspace(3.,4.2,11)
	ax.plot(x_values, np.zeros_like(x_values),
			color = 'darkgrey', linestyle = '-', marker = '', zorder = 3)
	y_values = np.linspace(0.,0.3,11)
	run_values = np.genfromtxt(data_dir/'run_values.csv',
								delimiter=',',dtype = int)
	area_values = np.genfromtxt(data_dir/'area_values.csv',
								delimiter=',', dtype = float)
	p_values = np.genfromtxt(data_dir/'p_values.csv',
								delimiter=',', dtype = float)
	if (data_dir / 'data.pkl').exists():
		with open(data_dir / 'data.pkl', 'rb') as filestream:
			data = pickle.load(filestream)
	else:
		data = np.zeros((len(run_values), len(area_values),
						 len(p_values), 0), dtype = float)
		for run_index, run_number in enumerate(run_values):
			run_dir = (data_dir/f'run_{run_number:02d}')
			if not run_dir.is_dir():
				continue
			for area_index,area_param in enumerate(area_values):
				for p_index,p_param in enumerate(p_values):
					try:
						csvfile = run_dir / ('eigenvalues' + \
									 '_area_{0:1.4f}'.format(area_param) + \
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
						data[run_index,area_index,p_index,:] = temp_data
					except:
						pass
		with open(data_dir / 'data.pkl', 'wb') as filestream:
			pickle.dump(data, filestream)
	total_modes = data.shape[3]
	floppy_modes = np.count_nonzero(data<1e-8 * \
					np.mean(data,axis=3)[:,:,:,np.newaxis],axis=3)
	floppy_modes = np.median(floppy_modes, axis=0)
	floppy_modes = floppy_modes/total_modes
	markers = np.array(['h', 'H', 'p', 'd', '^', '.'])
	colors = np.array(['tab:purple', 'tab:blue', 'tab:green',
						'tab:orange', 'tab:red', 'grey'])
	area_labels = np.array(['A = 2$\pi$','A = $\pi$','A = $\pi/2$',
							'A = $\pi/4$','A = $\pi/8$',
							r'A $\rightarrow 0$ (flat)'])
	x = np.linspace(3.6,4.4,1000)
#	ax.plot(x, F(x, 3.80541353, 0.2986759, 0.02161341),
#			color = 'grey',
#			linestyle = '-',
#			marker = '')
	curve_handles = np.empty(len(area_values), dtype=object)
	black_handles = np.empty(len(area_values), dtype=object)
	color_handles = np.empty(len(area_values), dtype=object)
	p0 = np.array([[3.96,0.198, 0.4, 0.4],
				   [3.94,0.200, 0.4, 0.4],
				   [3.91,0.202, 0.2, 0.2],
				   [3.90,0.210, 0.1, 0.1],
				   [3.90,0.215, 0.1, 0.1]])
	for area_index,area_param in enumerate(area_values):
		black_handles[area_index], = ax.plot(p_values,
				floppy_modes[area_index,:],
				color = 'black', #colors[area_index],
				marker = markers[area_index],
				markersize = 5.,
				linestyle = '',
				zorder = 10+area_index)
		color_handles[area_index], = ax.plot(p_values,
				floppy_modes[area_index,:],
				color = colors[area_index],
				marker = markers[area_index],
				markersize = 3.,
				linestyle = '',
				label = area_labels[area_index],
				zorder = 12+area_index)
#		params, covar = curve_fit(F,
#								  p_values,
#								  floppy_modes[area_index,:],
#								  p0 = [3.8, 0.25, 0.01])
		params, covar = curve_fit(g,
								  p_values,
								  floppy_modes[area_index,:],
								  p0 = p0[area_index])
		curve_handles[area_index], = ax.plot(x, g(x,*params),
				color = colors[area_index],
				linestyle = '-',
				marker = '',
				zorder = 6+area_index)
		print(params)
	if Path('bounded_eigen.csv').exists():
		bounded_data = np.genfromtxt(Path('bounded_eigen.csv'), delimiter=',')
		params, covar = curve_fit(g,
								  bounded_data[:,0],
								  bounded_data[:,1],
								  p0 = [3.88, 0.22, 0.58, 0.09])
		print(params)
		black_handles = np.append(black_handles, None)
		color_handles = np.append(color_handles, None)
		curve_handles = np.append(curve_handles, None)
		black_handles[-1], = ax.plot(
								bounded_data[:,0],
								bounded_data[:,1],
								color = 'black',
								marker = markers[-1],
								markersize = 6.,
								linestyle = '',
								zorder = 9)
		color_handles[-1], = ax.plot(
								bounded_data[:,0],
								bounded_data[:,1],
								color = colors[-1],
								marker = markers[-1],
								markersize = 3.,
								linestyle = '',
								label = area_labels[-1],
								zorder = 11)
		curve_handles[-1], = ax.plot(x, g(x,*params),
								color = colors[-1],
								linestyle = '-',
								marker = '',
								zorder = 5)
	ax.plot(np.ones_like(y_values)*3.812, y_values,
			color = 'black', linestyle = '--', marker = '',
			zorder = 4)
	ax.plot(x, F_flat(x),
			color = 'black', linestyle = '-', marker = '',
			label = r'N $\rightarrow \infty$ (flat)',
			zorder = 4)
	ax.set_ylim((-0.005,0.255))
#	ax.set_xlim((3.815,4.105))
	ax.set_xlim((3.795,4.105))
	ax.set_xlabel('$p_0$')
	ax.set_ylabel('$N(\lambda=0)$')
	handles, labels = ax.get_legend_handles_labels()
	combined_handles = handles
	combined_handles[:len(handles)-1] = zip(curve_handles,
										  black_handles,
										  color_handles)
	ax.legend(combined_handles,labels,
			  loc = 'upper left', fancybox = True,
			  framealpha = 1.)
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
						default = '../plots/hyperbolic_eigen_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('datadir',
						nargs = '?',
						default = '../hyperbolic_eigenvalues',
						type = str,
						help = 'directory of eigenvalues data')
	args = parser.parse_args()
	plot_eigen(data_dir = Path(args.datadir),
			   output_file = Path(args.outfile) )

################################################################################
# EOF
