#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess
import glob
from pathlib import Path

################################################################################
#===============================================================================
# analyse_sim.py
#===============================================================================
################################################################################

def combine_runs (out_dir = Path('../bounded_output') ):
#	for run_dir in out_dir.glob('quick_test'):
	for run_dir in out_dir.glob('run_[0-9]*'):
		if not run_dir.is_dir():
			continue
		print(str(run_dir))
		for csvfile in \
				run_dir.glob('r_[0-9]*.[0-9]*_p0_[0-9]*.[0-9]*_results.csv'):
			print(csvfile.name)
			r_param = float(str(csvfile.stem).split('_')[-4])
			p0_param = float(str(csvfile.stem).split('_')[-2])
			with open(csvfile, 'r') as instream, \
				 open(out_dir/'r_{0:2.1f}_p0_{1:1.3f}.csv'.format(
								r_param, p0_param), 'a') as outstream:
				for line in instream:
					outstream.write(line)

################################################################################

def compile_stats (out_dir = Path('../bounded_output') ):
	data = np.zeros((0,4), dtype = float)
	for csvfile in out_dir.glob('r_[0-9]*.[0-9]*_p0_[0-9]*.[0-9]*.csv'):
		r_param = float(str(csvfile.stem).split('_')[-3])
		p0_param = float(str(csvfile.stem).split('_')[-1])
		file_data = np.genfromtxt(csvfile, delimiter = ',', dtype = float)
		entry = np.array([ [r_param, p0_param,
							np.mean(file_data),
							np.std(file_data)/np.sqrt(len(file_data))] ])
		data = np.append(data, entry, axis=0)
	data.view('float,float,float,float').sort(order = ['f0', 'f1'], axis = 0)
	form = '%1.1f', '%1.3f', '%1.6f', '%1.6f'
	np.savetxt(out_dir / 'results.csv', data, delimiter = ',', fmt = form)

################################################################################

if __name__ == '__main__':
	combine_runs()
	compile_stats()

################################################################################
# EOF
