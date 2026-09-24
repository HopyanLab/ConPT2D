#!/usr/bin/env /usr/bin/python3
import numpy as np
import argparse
from pathlib import Path
from numpy.random import Generator, PCG64
from make_voronoi import make_voronoi

################################################################################
#===============================================================================
# make_initial.py
#===============================================================================
################################################################################

def make_initial (N = 64,
				  suevfile = Path('../spp_flat_output/initial_state.fe'),
				  shape_index = 3.8,
				  perimeter_modulus = 0.5,
				  length_threshold = 1.e-2,
				  energy_threshold = 1.e-8,
				  D_r = 1.0,
				  v_0 = 0.1,
				  dt = 0.02,
				  n_t = 1e5,
				  track_positions = False,
				  testing = False):
	vertices, edges, wrappings, faces = make_voronoi(N)
	# Write surface evolver script
	length_tolerance = 1e-9
	save_nums = np.concatenate(
				( np.linspace(1.,9.,num=9,endpoint=True), 
				  np.round(np.logspace(1,5,num=41,endpoint=True)) )
								).astype(int)
	with open(suevfile, 'w+') as outfile:
		outfile.write('// Periodic Voronoi Tesselation\n')
		outfile.write('\nSTRING')
		outfile.write('\nSURFACE_DIMENSION 1')
		outfile.write('\nSPACE_DIMENSION   2')
		outfile.write('\nTORUS\n')
		outfile.write('\nDEFINE save_nums integer[{0:d}]'.format(
															len(save_nums)))
		outfile.write(' = {')
		for index in np.arange(0,len(save_nums)):
			if index > 0:
				outfile.write(',')
			outfile.write('{0:d}'.format(save_nums[index]))
		outfile.write('}\n')
		outfile.write('\nDEFINE edge ATTRIBUTE flipped integer')
		outfile.write('\nDEFINE facet ATTRIBUTE n_x real')
		outfile.write('\nDEFINE facet ATTRIBUTE n_y real\n')
		outfile.write('\nPARAMETER p0_shape_index = {0:1.3f}'.format(
															shape_index))
		outfile.write('\nPARAMETER r_peri_modulus = {0:1.3f}'.format(
															perimeter_modulus))
		outfile.write('\nPARAMETER length_threshold = {0:1.12f}'.format(
															length_threshold))
		outfile.write('\nPARAMETER energy_threshold = {0:1.12f}'.format(
															energy_threshold))
		outfile.write('\nPARAMETER length_tolerance = {0:1.12f}'.format(
															length_tolerance))
		outfile.write('\nPARAMETER D_r = {0:1.12f}'.format(D_r))
		outfile.write('\nPARAMETER v_0 = {0:1.12f}'.format(v_0))
		outfile.write('\nPARAMETER dt = {0:1.12f}'.format(dt))
		outfile.write('\nPARAMETER n_t = {0:d}'.format(
													np.floor(n_t).astype(int)))
		outfile.write('\nPARAMETER track_positions = {0:b}'.format(
													track_positions))
		outfile.write('\nQUANTITY fixed_edge INFO_ONLY METHOD edge_length\n')
		for face_index, face in enumerate(faces):
			outfile.write('\nMETHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_peri METHOD edge_length\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos METHOD edge_torus_area\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg METHOD edge_torus_area\n')
			outfile.write('PARAMETER cell_{0:d}_corr '.format(face_index+1))
			outfile.write('= 0\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_energy ENERGY FUNCTION\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos.value - ')
			outfile.write('cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg.value + ')
			outfile.write('cell_{0:d}'.format(face_index+1))
			outfile.write('_corr - 1)^2 +\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_peri.value - p0_shape_index)^2')
			outfile.write('/r_peri_modulus\n')
		outfile.write('\nperiods\n')
		outfile.write('{0:f}\t0.0\n'.format(np.sqrt(N)))
		outfile.write('0.0\t{0:f}\n'.format(np.sqrt(N)))
		outfile.write('\nvertices\n')
		for vertex_index, vertex in enumerate(vertices):
			outfile.write('{0:d}\t{1:1.10f}\t{2:1.10f}\n'.format(
														vertex_index+1,
														vertex[0]*np.sqrt(N),
														vertex[1]*np.sqrt(N)))
		outfile.write('\nedges\n')
		for edge_index, edge in enumerate(edges):
			edge_wrapping = ['', '']
			for i in [0,1]:
				if wrappings[edge_index, i] == 0:
					edge_wrapping[i] = '*'
				elif wrappings[edge_index, i] == 1:
					edge_wrapping[i] = '+'
				else:
					edge_wrapping[i] = '-'
			# Write with Surface Evolver notation for wrapping
			outfile.write('{0:d}\t{1:d}\t{2:d}\t{3:s}\t{4:s}'.format(
										edge_index+1, edge[0], edge[1],
										edge_wrapping[0], edge_wrapping[1]))
			outfile.write('\ttension\t0')
			for face_index, face in enumerate(faces):
				for edge_check in face:
					if edge_check == edge_index+1:
						outfile.write('\tcell_{0:d}_peri'.format(face_index+1))
						outfile.write('\tcell_{0:d}'.format(face_index+1))
						outfile.write('_area_pos')
					elif -edge_check == edge_index+1:
						outfile.write('\tcell_{0:d}_peri'.format(face_index+1))
						outfile.write('\tcell_{0:d}'.format(face_index+1))
						outfile.write('_area_neg')
			outfile.write('\n')
		outfile.write('\nfaces\n')
		for face_index, face in enumerate(faces):
			outfile.write('{0:d}'.format(face_index+1))
			for index in face:
				outfile.write('\t{0:d}'.format(index))
			outfile.write('\n')
		outfile.write('\nbodies\n')
		for face_index, face in enumerate(faces):
			outfile.write('{0:d}\t{0:d}\n'.format(face_index+1))
		outfile.write('\nread\n')
		outfile.write('\nclipped')
		outfile.write('\nconj_grad on') # Much faster convergence!
		outfile.write('\nautorecalc on\n')
		rg = Generator(PCG64())
		outfile.write('\nrandom_seed := {0:d}\n'.format(
										rg.integers(0,1000)))
		# Procedures in seperate file.
		with open(Path(__file__).resolve().parent / \
					'procedures.fe','r') as commandfile:
			for line in commandfile.read():
				outfile.write(line)
		if testing:
			# Graphics command part.
			outfile.write('\nshow\nq\n')
			# Try to find energy minimum.
			outfile.write('\nrelax_system(10000);')
			outfile.write('\nrun_sim(n_t);')
		outfile.write('\n')

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Generate initial condition')
	# N = int(sys.argv[1])
	parser.add_argument('-n', '--number',
						nargs = 1,
						default = [64],
						type = int,
						required = False,
						help = 'number of points')
	parser.add_argument('-s', '--suevfile',
						nargs = 1,
						default = ['../spp_flat_output/initial_state.fe'],
						type = str,
						required = False,
						help = 'file to put surface evolver script in')
	parser.add_argument('-p', '--p0_param',
						nargs = 1,
						default = [3.8],
						type = float,
						required = False,
						help = 'shape index parameter')
	parser.add_argument('-r', '--r_param',
						nargs = 1,
						default = [0.5],
						type = float,
						required = False,
						help = 'inverse perimeter modulus')
	parser.add_argument('-t', '--track', dest='track_positions',
						action='store_const',
						const=True, default=False,
						help = 'flag to output positions during simulation')
	parser.add_argument('-x', '--test', dest='testing',
						action='store_const',
						const=True, default=False,
						help = 'flag to output positions during simulation')
	args = parser.parse_args()
	suevfile = Path(args.suevfile[0])
	suevfile.parent.mkdir(exist_ok = True)
	make_initial(N = args.number[0],
				 suevfile = suevfile,
				 shape_index = args.p0_param[0],
				 perimeter_modulus = args.r_param[0],
				 track_positions = args.track_positions,
				 testing = args.testing)

################################################################################
# EOF
