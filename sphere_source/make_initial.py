#!/usr/bin/env /usr/bin/python3
import numpy as np
from pathlib import Path
import argparse
from make_voronoi import make_voronoi

################################################################################
#===============================================================================
# make_initial.py
#===============================================================================
################################################################################

def make_initial (N = 64,
				  suevfile = Path('../sphere_output/initial_state.fe'),
				  shape_index = 3.8,
				  perimeter_modulus = 0.5,
				  length_threshold = 1.e-1,
				  energy_threshold = 1.e-8):
	# Radius of sphere is chosen so that surface area is number of points.
	R = np.sqrt(N/4/np.pi)
	# Make Voronoi diagram on sphere.
	vertices, edges, faces, regions = make_voronoi(N)
	# Write the Surface Evolver initial state
	length_tolerance = 1e-9
	constraint_tolerance = 1e-12
	with open(suevfile, 'w+') as outfile:
		outfile.write('// Spherical Voronoi Tesselation\n')
		outfile.write('\nSTRING')
		outfile.write('\nSURFACE_DIMENSION 1')
		outfile.write('\nSPACE_DIMENSION   3')
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
		outfile.write('\nCONSTRAINT_TOLERANCE {0:1.12f}'.format(
														constraint_tolerance))
		outfile.write('\nCONSTRAINT 1')
		outfile.write('\nformula: x**2 + y**2 + z**2 = radius**2\n')
		outfile.write('\nQUANTITY fixed_edge INFO_ONLY ')
		outfile.write('METHOD spherical_arc_length\n')
		for face_index, face in enumerate(faces):
			if np.sum(vertices[regions[face_index]], axis=0)[2] > 0:
				outfile.write('\nMETHOD_INSTANCE ')
				outfile.write('cell_{0:d}_area_pos '.format(face_index+1))
				outfile.write('METHOD spherical_arc_area_n\n')
				outfile.write('METHOD_INSTANCE ')
				outfile.write('cell_{0:d}_area_neg '.format(face_index+1))
				outfile.write('METHOD spherical_arc_area_n\n')
			else:
				outfile.write('\nMETHOD_INSTANCE ')
				outfile.write('cell_{0:d}_area_pos '.format(face_index+1))
				outfile.write('METHOD spherical_arc_area_s\n')
				outfile.write('METHOD_INSTANCE ')
				outfile.write('cell_{0:d}_area_neg '.format(face_index+1))
				outfile.write('METHOD spherical_arc_area_s\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_peri METHOD spherical_arc_length\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_energy ENERGY FUNCTION\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos.value - ')
			outfile.write('cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg.value - 1)^2 +\n')
			outfile.write('\t(cell_{0:d}'.format(face_index+1))
			outfile.write('_peri.value - p0_shape_index)^2')
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
			for entry in face:
				outfile.write('\t{0:d}'.format(entry))
			outfile.write('\n')
		outfile.write('\nbodies\n')
		for face_index, face in enumerate(faces):
			outfile.write('{0:d}\t{0:d}\n'.format(face_index+1))
		outfile.write('\nread\n')
		outfile.write('\nconj_grad on')
		outfile.write('\nautorecalc on\n')
		# Procedures in seperate file.
		with open(Path(__file__).resolve().parent / \
					'procedures.fe','r') as commandfile:
			for line in commandfile.read():
				outfile.write(line)
		# Graphics command part. Uncomment for testing.
	#	outfile.write('\nshow\nq\n')
		# Try to find energy minimum. Uncomment for testing.
	#	outfile.write('\nrelax_system(1000);')
	#	outfile.write('\nJ;\n0.01\nrelax_system(100);\nJ;')
	#	outfile.write('\nrelax_system(10000);\n')
	#	outfile.write('ritz(-100, 2*vertex_count);\n')

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
						default = ['../sphere_output/initial_state.fe'],
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
	args = parser.parse_args()
	suevfile = Path(args.suevfile[0])
	suevfile.parent.mkdir(exist_ok = True)
	make_initial(N = args.number[0],
				 suevfile = suevfile,
				 shape_index = args.p0_param[0],
				 perimeter_modulus = args.r_param[0])

################################################################################
# EOF
