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
				  suevfile = Path('../spp_sphere_output/initial_state.fe'),
				  shape_index = 3.8,
				  perimeter_modulus = 0.5,
				  length_threshold = 1.e-2,
				  energy_threshold = 1.e-8,
				  D_r = 4.0,
				  v_0 = 0.1,
				  dt = 0.02,
				  n_t = 1e5,
				  track_positions = False,
				  testing = False):
	# Radius of sphere is chosen so that surface area is number of points.
	R = np.sqrt(N/4/np.pi)
	# Make Voronoi diagram on sphere.
	vertices, edges, faces, regions = make_voronoi(N)
	# Write the Surface Evolver initial state
	length_tolerance = 1.e-9
	save_nums = np.unique(np.concatenate(
				( np.linspace(1.,9.,num=9,endpoint=True), 
				  np.round(np.logspace(1,5,num=41,endpoint=True)) )
								).astype(int))
	constraint_tolerance = 1.e-12
	with open(suevfile, 'w+') as outfile:
		outfile.write('// Spherical Voronoi Tesselation\n')
		outfile.write('\nSTRING')
		outfile.write('\nSURFACE_DIMENSION 1')
		outfile.write('\nSPACE_DIMENSION   3')
		outfile.write('\nDEFINE save_nums integer[{0:d}]'.format(
															len(save_nums)))
		outfile.write(' = {')
		for index in np.arange(0,len(save_nums)):
			if index > 0:
				outfile.write(',')
			outfile.write('{0:d}'.format(save_nums[index]))
		outfile.write('}\n')
		outfile.write('\nDEFINE edge ATTRIBUTE flipped integer')
		outfile.write('\nDEFINE facet ATTRIBUTE use_north integer')
		outfile.write('\nDEFINE facet ATTRIBUTE last_normal_x real')
		outfile.write('\nDEFINE facet ATTRIBUTE last_normal_y real')
		outfile.write('\nDEFINE facet ATTRIBUTE last_normal_z real')
		outfile.write('\nDEFINE facet ATTRIBUTE normal_x real')
		outfile.write('\nDEFINE facet ATTRIBUTE normal_y real')
		outfile.write('\nDEFINE facet ATTRIBUTE normal_z real')
		outfile.write('\nDEFINE facet ATTRIBUTE direction_x real')
		outfile.write('\nDEFINE facet ATTRIBUTE direction_y real')
		outfile.write('\nDEFINE facet ATTRIBUTE direction_z real')
		outfile.write('\nDEFINE vertex ATTRIBUTE motion_x real')
		outfile.write('\nDEFINE vertex ATTRIBUTE motion_y real')
		outfile.write('\nDEFINE vertex ATTRIBUTE motion_z real\n')
		#outfile.write('\nSCALE 0.001 FIXED')
		outfile.write('\nLENGTH_METHOD_NAME "spherical_arc_length"')
		#outfile.write('\nAREA_METHOD_NAME "spherical_arc_area_n"\n')
		outfile.write('\nPARAMETER p0_shape_index   = {0:1.3f}'.format(
															shape_index))
		outfile.write('\nPARAMETER r_peri_modulus   = {0:1.3f}'.format(
															perimeter_modulus))
		outfile.write('\nPARAMETER length_threshold = {0:1.12f}'.format(
															length_threshold))
		outfile.write('\nPARAMETER energy_threshold = {0:1.12f}'.format(
															energy_threshold))
		outfile.write('\nPARAMETER length_tolerance = {0:1.12f}'.format(
															length_tolerance))
		outfile.write('\nPARAMETER radius           = {0:f}\n'.format(R))
		outfile.write('\nPARAMETER D_r = {0:1.12f}'.format(D_r))
		outfile.write('\nPARAMETER v_0 = {0:1.12f}'.format(v_0))
		outfile.write('\nPARAMETER dt = {0:1.12f}'.format(dt))
		outfile.write('\nPARAMETER n_t = {0:d}'.format(
													np.floor(n_t).astype(int)))
		outfile.write('\nPARAMETER track_positions = {0:b}'.format(
													track_positions))
		outfile.write('\nCONSTRAINT_TOLERANCE {0:1.12f}'.format(
														constraint_tolerance))
		outfile.write('\nCONSTRAINT 1')
		outfile.write('\nformula: x**2 + y**2 + z**2 = radius**2\n')
		outfile.write('\nQUANTITY fixed_edge INFO_ONLY ')
		outfile.write('METHOD spherical_arc_length\n')
		for face_index, face in enumerate(faces):
			outfile.write('\nMETHOD_INSTANCE ')
			outfile.write('cell_{0:d}_area_pos_n '.format(face_index+1))
			outfile.write('METHOD spherical_arc_area_n\n')
			outfile.write('METHOD_INSTANCE ')
			outfile.write('cell_{0:d}_area_neg_n '.format(face_index+1))
			outfile.write('METHOD spherical_arc_area_n\n')
			outfile.write('METHOD_INSTANCE ')
			outfile.write('cell_{0:d}_area_pos_s '.format(face_index+1))
			outfile.write('METHOD spherical_arc_area_s\n')
			outfile.write('METHOD_INSTANCE ')
			outfile.write('cell_{0:d}_area_neg_s '.format(face_index+1))
			outfile.write('METHOD spherical_arc_area_s\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_peri METHOD spherical_arc_length\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_area INFO_ONLY FUNCTION\n')
			outfile.write('\tcell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos_n.value + ')
			outfile.write('cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos_s.value - \n')
			outfile.write('\tcell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg_n.value - ')
			outfile.write('cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg_s.value\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_energy ENERGY FUNCTION\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_area.value - 1)^2 +\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_peri.value - p0_shape_index)^2')
#			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
#			outfile.write('_energy ENERGY FUNCTION\n')
#			outfile.write('\t(cell_{0:d}'.format(face_index+1))
#			outfile.write('_area_pos_n.value + ')
#			outfile.write('cell_{0:d}'.format(face_index+1))
#			outfile.write('_area_pos_s.value - \n')
#			outfile.write('\tcell_{0:d}'.format(face_index+1))
#			outfile.write('_area_neg_n.value - ')
#			outfile.write('cell_{0:d}'.format(face_index+1))
#			outfile.write('_area_neg_s.value - 1)^2 +\n')
#			outfile.write('\t(cell_{0:d}'.format(face_index+1))
#			outfile.write('_peri.value - p0_shape_index)^2')
			outfile.write('/r_peri_modulus\n')
		outfile.write('\nvertices\n')
		for vertex_index, vertex in enumerate(vertices):
			outfile.write('{0:d}\t{1:1.10f}\t{2:1.10f}\t{3:1.10f}\t'.format(
														vertex_index+1,
														vertex[0],
														vertex[1],
														vertex[2]) + \
							'constraint 1\n')
		outfile.write('\nedges\n')
		for edge_index, edge in enumerate(edges):
			outfile.write('{0:d}\t{1:d}\t{2:d}'.format(
										edge_index+1, edge[0]+1, edge[1]+1))
			outfile.write('\ttension\t0')
			for face_index, face in enumerate(faces):
				if np.sum(vertices[regions[face_index]], axis=0)[2] > 0:
					for edge_check in face:
						if edge_check == edge_index+1:
							outfile.write('\tcell_{0:d}_peri'.format(
																face_index+1))
							outfile.write('\tcell_{0:d}'.format(face_index+1))
							outfile.write('_area_pos_n')
						elif -edge_check == edge_index+1:
							outfile.write('\tcell_{0:d}_peri'.format(
																face_index+1))
							outfile.write('\tcell_{0:d}'.format(face_index+1))
							outfile.write('_area_neg_n')
				else:
					for edge_check in face:
						if edge_check == edge_index+1:
							outfile.write('\tcell_{0:d}_peri'.format(
																face_index+1))
							outfile.write('\tcell_{0:d}'.format(face_index+1))
							outfile.write('_area_pos_s')
						elif -edge_check == edge_index+1:
							outfile.write('\tcell_{0:d}_peri'.format(
																face_index+1))
							outfile.write('\tcell_{0:d}'.format(face_index+1))
							outfile.write('_area_neg_s')
			outfile.write('\n')
		outfile.write('\nfaces\n')
		for face_index, face in enumerate(faces):
			outfile.write('{0:d}'.format(face_index+1))
			for entry in face:
				outfile.write('\t{0:d}'.format(entry))
			if np.sum(vertices[regions[face_index]], axis=0)[2] > 0:
				outfile.write('\tuse_north\t1')
			else:
				outfile.write('\tuse_north\t0')
			outfile.write('\n')
		outfile.write('\nbodies\n')
		for face_index, face in enumerate(faces):
			outfile.write('{0:d}\t{0:d}\n'.format(face_index+1))
		outfile.write('\nread\n')
		outfile.write('\nconj_grad on')
		outfile.write('\nautorecalc on\n')
		rg = Generator(PCG64())
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
						default = ['../spp_sphere_output/initial_state.fe'],
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
