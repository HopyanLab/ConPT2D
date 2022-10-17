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
# run_sim.py
#===============================================================================
################################################################################

base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'sphere_output'
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 4
r_param = 0.5
n_values = np.array([24, 32, 48, 64, 128])
#n_values = np.array([24, 32, 48])
run_type = 'coarse'
#run_type = 'fine'

# Coarse
if run_type == 'coarse':
	p0_values = np.arange(3.60, 3.91, 0.02)
# Fine
elif run_type == 'fine':
	p0_values = np.arange(3.65, 3.821, 0.002)
#	p0_values = np.arange(3.720, 3.821, 0.002)
else:
	p0_values = np.arange(3.6, 4.01, 0.1)

def run_sim (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Run through the parameter values.
	for n_param in n_values:
		if (sim_dir/('initial_state_{0:d}.fe'.format(n_param))).exists():
			pass
		else:
			# Generate an initial state.
			make_initial( N = n_param,
						  suevfile = sim_dir / ('initial_state_' + \
												'{0:d}.fe'.format(n_param)),
						  shape_index = 3.6,
						  perimeter_modulus = r_param )
		for p0_param in p0_values:
			if (sim_dir/('relaxed_state_n_{0:d}'.format(n_param) + \
						'_p0_{0:1.3f}.fe'.format(p0_param))).exists():
				continue
			else:
				# Make a simulation script.
				with open(sim_dir/'sim.fe','w') as sim_script:
					sim_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p0_param))
					sim_script.write('relax_system(100);\n')
					sim_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
					sim_script.write('relax_system(10000);\n')
					sim_script.write('run_sim();\n')
					sim_script.write('dump "relaxed_state')
					sim_script.write('_n_{0:d}'.format(n_param))
					sim_script.write('_p0_{0:1.3f}'.format(p0_param))
					sim_script.write('.fe";\n')
					sim_script.write('quit 1\n')
				# Run simulation.
				sim = sp.Popen(['evolver','-fsim.fe','-x',
								'initial_state_{0:d}.fe'.format(n_param)])
				sim.wait()
				(sim_dir / 'sim.fe').unlink()

################################################################################

if __name__ == '__main__':
	code_timer = timer()
	code_timer.start()
	Path.mkdir(out_dir, exist_ok = True)
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_sim, range(1,number_sims+1))
	code_timer.stop()

################################################################################
# EOF
