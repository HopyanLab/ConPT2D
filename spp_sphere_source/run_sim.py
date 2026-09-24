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
out_dir = base_dir.parent/'spp_sphere_output'
cores_to_use = mp.cpu_count() - 4 # = 12
number_sims = cores_to_use * 8
r_param = 0.5
#n_values = np.array([32, 64])
n_values = np.array([24, 32, 48, 64, 128])

#p0_values = np.arange(3.6, 4.01, 0.025)
#p0_values = np.array([3.6, 3.7, 3.80, 3.9, 4.0])
#p0_values = np.array([3.775, 3.8, 3.825])
p0_values = np.array([3.8])

def run_sim (run_number):
	sim_dir = out_dir/'run_{0:02d}'.format(run_number)
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Run through the parameter values.
	for n_param in n_values:
		n_dir = sim_dir/('n_{0:d}'.format(n_param))
		Path.mkdir(n_dir, exist_ok = True)
		if (n_dir/('initial_state_{0:d}.fe'.format(n_param))).exists():
			pass
		else:
			# Generate an initial state.
			make_initial( N = n_param,
						  suevfile = n_dir / ('initial_state_' + \
												'{0:d}.fe'.format(n_param)),
						  shape_index = 3.6,
						  perimeter_modulus = 1.0,
						  D_r = 4.0,
						  v_0 = 0.1,
						  dt = 0.02,
						  n_t = 1e5,
						  track_positions = True)
		for p0_param in p0_values:
			p0_dir = n_dir/('p0_{0:1.3f}'.format(p0_param))
			Path.mkdir(p0_dir, exist_ok = True)
			os.chdir(str(p0_dir))
			if (p0_dir/('n_{0:d}'.format(n_param) + \
						'_p0_{0:1.3f}.fe'.format(p0_param))).exists():
				continue
			else:
				# Make a simulation script.
				with open(p0_dir/'sim.fe','w') as sim_script:
					sim_script.write('p0_shape_index := {0:1.3f};\n'.format(
																	p0_param))
					sim_script.write('relax_system(1000);\n')
					sim_script.write('run_sim(n_t);\n')
					sim_script.write('dump ')
					sim_script.write('"n_{0:d}'.format(n_param))
					sim_script.write('_p0_{0:1.3f}'.format(p0_param))
					sim_script.write('.fe";\n')
					sim_script.write('quit 1\n')
				# Run simulation.
				with open(os.devnull, 'w') as nowhere:
					sim = sp.Popen(['evolver','-fsim.fe','-x',
						str(n_dir)+'/initial_state_{0:d}.fe'.format(n_param)],
											stdout=nowhere)
					sim.wait()
				(p0_dir / 'sim.fe').unlink()
			os.chdir(str(sim_dir))

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
