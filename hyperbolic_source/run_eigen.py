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
# run_eigen.py
#===============================================================================
################################################################################

base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'hyperbolic_eigenvalues'
out_dir.mkdir(exist_ok=True)
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 16
r_param = 0.5
number_cells = 128
polygon_number = 12
area_values = np.array([2.0,1.0,0.5,0.25,0.125])*np.pi
#run_type = 'coarse'
run_type = 'medium'
#run_type = 'fine'

# Coarse
if run_type == 'coarse':
	p_values = np.arange(3.800, 4.201, 0.020)
# Fine
elif run_type == 'fine':
	p_values = np.arange(3.800, 4.201, 0.002)
else:
	p_values = np.arange(3.800, 4.201, 0.010)

def run_eigen (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Run through the parameter values.
	for area_param in area_values:
		# Generate an initial state.
		if (sim_dir/('initial_state_{0:1.4f}.fe'.format(area_param))).exists():
			pass
		else:
			make_initial(N = number_cells,
						 polygon_sides = polygon_number,
						 outfile = sim_dir / ('initial_state_' + \
									'{0:1.4f}.fe'.format(area_param)),
						 polygon_area = area_param,
						 shape_index = 3.6,
						 perimeter_modulus = 0.5)
		for p_param in p_values:
			if (sim_dir/('eigenvalues_area_{0:1.4f}'.format(area_param) + \
						'_p0_{0:1.3f}.csv'.format(p_param))).exists():
				continue
			else:
				with open(sim_dir/'eigen.fe','w') as eigen_script:
					eigen_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p_param))
					eigen_script.write('relax_system(1000);\n')
					eigen_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
					eigen_script.write('relax_system(10000);\n')
					eigen_script.write('ritz(-1000,2*vertex_count)')
					eigen_script.write('>>"temp.txt"\n')
					eigen_script.write('quit 1\n')
				# Relax system and output eigenvalues.
				eigen = sp.Popen(['evolver', '-feigen.fe', '-x',
							'initial_state_{0:1.4f}.fe'.format(area_param)])
				eigen.wait()
				(sim_dir / 'eigen.fe').unlink()
				print(run_number)
				eigenvalues = np.genfromtxt('./temp.txt', usecols = (1),
											skip_header = 2, skip_footer = 1)
				np.savetxt(sim_dir / \
							('eigenvalues_area_{0:1.4f}'.format(area_param) + \
							 '_p0_{0:1.3f}.csv'.format(p_param)),
							eigenvalues,
							delimiter = ',')
				(sim_dir / 'temp.txt').unlink()

if __name__ == '__main__':
	os.chdir(str(out_dir))
	np.savetxt('./area_values.csv', area_values, delimiter=',', fmt='%1.4f')
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
