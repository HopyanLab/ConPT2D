#!/usr/bin/env /usr/bin/python3
import numpy as np
import os
import subprocess as sp
import multiprocessing as mp
import shutil
from pathlib import Path
from make_voronoi import *

################################################################################
#===============================================================================
# run_test.py
#===============================================================================
################################################################################

cores_to_use = 11 # mp.cpu_count() - 2
number_cells = 64
r_value = 0.5
p0_values = np.arange(3.7, 3.91, 0.02)
base_dir = Path(__file__).resolve().parent
out_dir = base_dir.parent/'flat_output'

def run_test (p0_value):
	# Make a simulation script.
	script_name = 'sim_{0:1.2f}.fe'.format(p0_value)
	with open(sim_dir/script_name,
					'w') as sim_script:
		sim_script.write('r_peri_modulus := {0:f};\n'.format(r_value))
		sim_script.write('p0_shape_index := {0:f};\n'.format(p0_value))
		sim_script.write('relax_system(10000);\n')
		sim_script.write('J;\n0.01\nrelax_system(100);\nJ;\n')
		sim_script.write('relax_system(1000);\n')
		sim_script.write('run_sim();\n')
		sim_script.write('quit 1\n')
	# Give'er
	sim = sp.Popen(['evolver', '-f' + script_name, '-x', 'initial_state.fe'])
	sim.wait()
	(sim_dir / 'sim_{0:1.2f}.fe'.format(p0_value)).unlink()

################################################################################

if __name__ == '__main__':
	Path.mkdir(out_dir, exist_ok = True)
	sim_dir = out_dir/'quick_test'
	Path.mkdir(sim_dir, exist_ok = True)
	os.chdir(str(sim_dir))
	# Generate an initial state.
	make_voronoi( N = number_cells,
				  vertfile = sim_dir / 'vertices.txt',
				  edgefile = sim_dir / 'edges.txt',
				  facefile = sim_dir / 'faces.txt',
				  suevfile = sim_dir / 'initial_state.fe',
				  shape_index = 3.0,
				  perimeter_modulus = 0.5 )
	with mp.Pool(processes = cores_to_use) as pool:
		pool.map(run_test, p0_values)

################################################################################
# EOF
