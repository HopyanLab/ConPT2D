#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess as sp
import multiprocessing as mp
from pathlib import Path
from timer import timer
from make_voronoi import *

################################################################################
#===============================================================================
# run_eigen.py
#===============================================================================
################################################################################

base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'flat_eigenvalues'
out_dir.mkdir(exist_ok=True)
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 16
number_cells = 64
#run_type = 'coarse'
run_type = 'fine'
r_values = np.array([0.5, 1., 2., 8.])

# Coarse
if run_type == 'coarse':
	p_values = np.arange(3.5, 3.91, 0.02)
# Fine
elif run_type == 'fine':
	p_values = np.arange(3.720, 3.881, 0.002)
else:
	p_values = np.arange(3.4, 4.01, 0.1)

def run_eigen (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Generate an initial state.
	if (sim_dir/'initial_state.fe').exists():
		pass
	else:
		make_voronoi( N = number_cells,
					  vertfile = sim_dir / 'vertices.txt',
					  edgefile = sim_dir / 'edges.txt',
					  facefile = sim_dir / 'faces.txt',
					  suevfile = sim_dir / 'initial_state.fe',
					  shape_index = 3.0,
					  perimeter_modulus = 0.5 )
	
	# Run through the parameter values.
	for r_param in r_values:
		for p_param in p_values:
			if (sim_dir/('eigenvalues_r_{0:1.1f}'.format(r_param) + \
						'_p0_{0:1.3f}.csv'.format(p_param))).exists():
				continue
			else:
				with open(sim_dir/'eigen.fe','w') as eigen_script:
					eigen_script.write('r_peri_modulus := {0:1.3f};\n'.format(
																	r_param))
					eigen_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p_param))
					eigen_script.write('relax_system(10000);\n')
					eigen_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
					eigen_script.write('relax_system(10000);\n')
					eigen_script.write('ritz(-1,2*vertex_count)>>"temp.txt"\n')
					eigen_script.write('quit 1\n')
				# Relax system and output eigenvalues.
				eigen = sp.Popen(['evolver','-feigen.fe',
								 '-x','initial_state.fe'])
				eigen.wait()
				(sim_dir / 'eigen.fe').unlink()
				eigenvalues = np.genfromtxt('./temp.txt', usecols = (1),
											skip_header = 2, skip_footer = 1)
				np.savetxt(sim_dir / \
							('eigenvalues_r_{0:1.1f}'.format(r_param) + \
							 '_p0_{0:1.3f}.csv'.format(p_param)),
							eigenvalues,
							delimiter = ',')
				(sim_dir / 'temp.txt').unlink()

if __name__ == '__main__':
	os.chdir(str(out_dir))
	np.savetxt('./r_values.csv', r_values, delimiter=',', fmt='%1.3f')
	np.savetxt('./p_values.csv', p_values, delimiter=',', fmt='%1.3f')
	np.savetxt('./run_values.csv', np.arange(1,number_sims+1),
				delimiter=',', fmt='%d')
	code_timer = timer()
	code_timer.start()
	Path.mkdir(out_dir, exist_ok = True)
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_eigen, range(1,number_sims+1))
	code_timer.stop()

################################################################################
# EOF
