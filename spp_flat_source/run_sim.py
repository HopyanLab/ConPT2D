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
out_dir = base_dir.parent/'spp_flat_output'
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 8
number_cells = 64

p0_values = np.arange(3.6, 4.01, 0.025)

def run_sim (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Generate an initial state.
	if (sim_dir/'initial_state.fe').exists():
		pass
	else:
		make_initial( N = number_cells,
					  suevfile = sim_dir / 'initial_state.fe',
					  shape_index = 3.6,
					  perimeter_modulus = 1.0,
					  D_r = 1.0,
					  v_0 = 0.1,
					  dt = 0.02,
					  n_t = 1e5,
					  track_positions = True)
	# Run through the parameter values.
	for p0_param in p0_values:
		p0_dir = sim_dir/('p0_{0:1.3f}'.format(p0_param))
		Path.mkdir(p0_dir, exist_ok = True)
		os.chdir(str(p0_dir))
		if (p0_dir/('p0_{0:1.3f}.fe'.format(p0_param))).exists():
			continue
		else:
			# Make a simulation script.
			with open(p0_dir/'sim.fe','w') as sim_script:
				sim_script.write('p0_shape_index := {0:1.3f};\n'.format(
																p0_param))
				sim_script.write('relax_system(1000);\n')
				sim_script.write('run_sim(n_t);\n')
				sim_script.write('dump "p0_{0:1.3f}.fe";\n'.format(p0_param))
				sim_script.write('quit 1\n')
			# Run simulation.
			with open(os.devnull, 'w') as nowhere:
				sim = sp.Popen(['evolver','-fsim.fe',
								'-x','-y','../initial_state.fe'],
										stdout=nowhere)
				sim.wait()
			(p0_dir / 'sim.fe').unlink()
		os.chdir(str(sim_dir))

if __name__ == '__main__':
	code_timer = timer()
	code_timer.start()
	Path.mkdir(out_dir, exist_ok = True)
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_sim, range(1,number_sims+1))
	code_timer.stop()

################################################################################
# EOF
