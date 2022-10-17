#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess as sp
import multiprocessing as mp
from pathlib import Path
from timer import timer
from make_initial import *

################################################################################
#===============================================================================
# run_sim.py
#===============================================================================
################################################################################

base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'hyperbolic_output'
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 4
r_param = 0.5
number_cells = 128
polygon_number = 12
area_values = np.array([2.0,1.0,0.5,0.25,0.125])*np.pi
#run_type = 'coarse'
run_type = 'medium'
#run_type = 'fine'

# Coarse
if run_type == 'coarse':
	p0_values = np.arange(3.8, 4.01, 0.1)
# Medium
elif run_type == 'medium':
	p0_values = np.arange(3.80, 4.01, 0.02)
# Fine
elif run_type == 'fine':
	p0_values = np.arange(3.800, 4.001, 0.002)
else:
	p0_values = np.arange(3.0, 4.01, 0.1)

p0_values = np.array([3.2,3.4,3.6])

def run_sim (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Run through the parameter values.
	for area_param in area_values:
		# Generate an initial state.
		if (sim_dir/('initial_state_{0:1.2f}.fe'.format(area_param))).exists():
			pass
		else:
			make_initial(N = number_cells,
						 polygon_sides = polygon_number,
						 outfile = sim_dir / ('initial_state_' + \
									'{0:1.2f}.fe'.format(area_param)),
						 polygon_area = area_param,
						 shape_index = 3.6,
						 perimeter_modulus = 0.5)
		for p0_param in p0_values:
			if (sim_dir/('relaxed_state_area_{0:1.2f}'.format(area_param) + \
						'_p0_{0:1.3f}.fe'.format(p0_param))).exists():
				continue
			else:
				# Make a simulation script.
				with open(sim_dir/'sim.fe','w') as sim_script:
					sim_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p0_param))
					sim_script.write('relax_system(1000);\n')
					sim_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
					sim_script.write('relax_system(10000);\n')
					sim_script.write('run_sim();\n')
					sim_script.write('dump "relaxed_state')
					sim_script.write('_area_{0:1.2f}'.format(area_param))
					sim_script.write('_p0_{0:1.3f}'.format(p0_param))
					sim_script.write('.fe";\n')
					sim_script.write('quit 1\n')
				# Run simulation.
				sim = sp.Popen(['evolver', '-fsim.fe', '-x',
							'initial_state_{0:1.2f}.fe'.format(area_param)])
				sim.wait()
				(sim_dir / 'sim.fe').unlink()

if __name__ == '__main__':
	code_timer = timer()
	code_timer.start()
	Path.mkdir(out_dir, exist_ok = True)
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_sim, range(1,number_sims+1))
	code_timer.stop()

################################################################################
# EOF
