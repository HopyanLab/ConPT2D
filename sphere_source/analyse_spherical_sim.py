#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess
import glob
from pathlib import Path
import argparse

################################################################################
#===============================================================================
# analyse_spherical_sim.py
#===============================================================================
################################################################################

def combine_runs (out_dir = Path('../sphere_output') ):
	for run_dir in out_dir.glob('run_[0-9]*'):
		if not run_dir.is_dir():
			continue
	#	print(str(run_dir))
		for csvfile in run_dir.glob('n_[0-9]*_p0_[0-9]*.[0-9]*_results.csv'):
	#		print(csvfile.name)
			n_param = int(str(csvfile.stem).split('_')[-4])
			p0_param = float(str(csvfile.stem).split('_')[-2])
			data = np.genfromtxt(csvfile)
			mean = np.mean(data)
			median = np.median(data)
			width = np.std(data)/np.sqrt(len(data))
			with open(out_dir/'n_{0:d}_p0_{1:1.3f}.csv'.format(
								n_param, p0_param), 'a') as data_file:
					np.savetxt(data_file,data)
			with open(out_dir/'stats_n_{0:d}_p0_{1:1.3f}.csv'.format(
								n_param, p0_param), 'a') as stats_file:
					stats_file.write(
						f'{median:1.9f},{mean:1.9f},{width:1.9f}\n')

################################################################################

def compile_stats (out_dir = Path('../sphere_output') ):
	data = np.zeros((0,5), dtype = float)
	for csvfile in out_dir.glob('n_[0-9]*_p0_[0-9]*.[0-9]*.csv'):
	#	print(csvfile.name)
		n_param = float(str(csvfile.stem).split('_')[-3])
		p0_param = float(str(csvfile.stem).split('_')[-1])
		file_data = np.genfromtxt(csvfile, delimiter = ',', dtype = float)
		stat_data = np.genfromtxt(csvfile.with_name('stats_' + csvfile.name),
									delimiter = ',', dtype = float)
		entry = np.array([ [n_param, p0_param,
							np.median(file_data),
							np.mean(stat_data[:,1]),
							np.sqrt(np.sum(stat_data[:,2]**2)) / \
									np.sqrt(len(stat_data[:,1]))] ])
		data = np.append(data, entry, axis=0)
	data.view('int,float,float,float,float').sort(
								order = ['f0', 'f1'], axis = 0)
	form = '%d', '%1.3f', '%1.6f', '%1.6f', '%1.6f'
	np.savetxt(out_dir / 'results.csv', data, delimiter = ',', fmt = form)

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Combine simulation results.')
	parser.add_argument('data_path',
						nargs = '?',
						default = '../sphere_output',
						type = str,
						help = 'path to simulation result data')
	args = parser.parse_args()
	out_dir = Path(args.data_path)
	for csvfile in out_dir.glob('n_[0-9]*_p0_[0-9]*.[0-9]*.csv'):
		csvfile.unlink()
	for csvfile in out_dir.glob('stats_n_[0-9]*_p0_[0-9]*.[0-9]*.csv'):
		csvfile.unlink()
	(out_dir/'results.csv').unlink
	combine_runs(out_dir)
	compile_stats(out_dir)
	subprocess.Popen(['python3','./plot_spherical_output.py',
						str(out_dir/'results.csv')])
	subprocess.Popen(['python3','./plot_spherical_dist.py',str(out_dir)])

################################################################################
# EOF
