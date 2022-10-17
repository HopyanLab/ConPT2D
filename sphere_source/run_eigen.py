#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess as sp
import multiprocessing as mp
from pathlib import Path
from timer import timer
from make_initial import make_initial

################################################################################
#===============================================================================
# run_spherical_eigen.py
#===============================================================================
################################################################################

base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'sphere_eigenvalues'
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 16
r_param = 0.5
n_values = np.array([24, 32, 48, 64, 128])
#run_type = 'coarse'
run_type = 'fine'
# Coarse
if run_type == 'coarse':
	p_values = np.arange(3.70, 3.91, 0.02)
# Fine
elif run_type == 'fine':
	p_values = np.arange(3.720, 3.881, 0.002)
else:
	p_values = np.arange(3.6, 4.01, 0.1)
out_dir.mkdir(exist_ok = True)
os.chdir(str(out_dir))
np.savetxt('./n_values.csv', n_values, delimiter=',', fmt='%d')
np.savetxt('./p_values.csv', p_values, delimiter=',', fmt='%1.3f')
np.savetxt('./run_values.csv', np.arange(1,number_sims+1),
			delimiter=',', fmt='%d')

def run_eigen (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Run through the parameter values.
	for n_param in n_values:
		# Generate an initial state.
		if (sim_dir/('initial_state_{0:d}.fe'.format(n_param))).exists():
			pass
		else:
			make_initial( N = n_param,
						  suevfile = sim_dir / ('initial_state_' + \
												'{0:d}.fe'.format(n_param)),
						  shape_index = 3.0,
						  perimeter_modulus = r_param )
		for p_param in p_values:
			if (sim_dir/('eigenvalues_n_{0:d}'.format(n_param) + \
						'_p0_{0:1.3f}.csv'.format(p_param))).exists():
				continue
			else:
				with open(sim_dir/'eigen.fe','w') as eigen_script:
					eigen_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p_param))
					eigen_script.write('relax_system(1000);\n')
					eigen_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
					eigen_script.write('relax_system(10000);\n')
					eigen_script.write('ritz(-1,2*vertex_count)>>"temp.txt"\n')
					eigen_script.write('quit 1\n')
				# Relax system and output eigenvalues.
				eigen = sp.Popen(['evolver','-feigen.fe',
								 '-x','initial_state_{0:d}.fe'.format(n_param)])
				eigen.wait()
				(sim_dir / 'eigen.fe').unlink()
				eigenvalues = np.genfromtxt('./temp.txt', usecols = (1),
											skip_header = 2, skip_footer = 1)
				np.savetxt(sim_dir / \
							('eigenvalues_n_{0:d}'.format(n_param) + \
							 '_p0_{0:1.3f}.csv'.format(p_param)),
							eigenvalues,
							delimiter = ',')
				(sim_dir / 'temp.txt').unlink()

if __name__ == '__main__':
	code_timer = timer()
	code_timer.start()
	Path.mkdir(out_dir, exist_ok = True)
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_eigen, range(1,number_sims+1))
	code_timer.stop()

################################################################################
# EOF
