#!/usr/bin/env /usr/bin/python3
import numpy as np
from scipy.spatial import Voronoi
import argparse
from pathlib import Path
from make_bounded_voronoi import make_bounded_voronoi
from plot_system import plot_system

################################################################################
#===============================================================================
# make_initial.py
#===============================================================================
################################################################################

def make_initial (N = 64,
				  polygon_sides = 6,
				  outfile = Path('../bounded_output/initial_state.fe'),
				  scale_factor = None,
				  shape_index = 3.8,
				  perimeter_modulus = 0.5,
				  length_threshold = 1.e-1,
				  energy_threshold = 1.e-8,
				  make_plot = False):
	if scale_factor is None:
		scale_factor = np.sqrt(N * 2 / polygon_sides / \
							np.sin(2*np.pi/polygon_sides))
	# generate Voronoi diagram in polygon
	vertices, edges, faces, boundary_vertices, boundary_edges, \
			outer_vertices, outer_edges, outer_faces = \
					make_bounded_voronoi(N, polygon_sides, scale_factor,
												1.e-6)
	# make plot of initial geometry
	if make_plot:
		from matplotlib import pyplot as plt
		fig = plt.figure()
		ax = plt.axes()
		ax = plot_system(ax, vertices, edges, faces,
							boundary_edges, vertices[outer_vertices])
		plt.savefig(outfile.parent / 'voronoi.svg')
	# boundaries are lines forming outside of polygon.
	angle = 2*np.pi/polygon_sides
	distance = np.cos(angle/2)*scale_factor
	angles = (np.arange(polygon_sides)+1/2) * angle
	x_coeffs = np.cos(angles)
	y_coeffs = np.sin(angles)
	length_tolerance = 1e-9
	# write surface evolver script
	with open(outfile, 'w+') as outfile:
		outfile.write('// Bounded Voronoi Tesselation\n')
		outfile.write('\nSTRING')
		outfile.write('\nSURFACE_DIMENSION 1')
		outfile.write('\nSPACE_DIMENSION   2')
		outfile.write('\nDEFINE facet ATTRIBUTE type integer\n')
		outfile.write('\nPARAMETER p0_shape_index = {0:1.3f}'.format(
															shape_index))
		outfile.write('\nPARAMETER r_peri_modulus = {0:1.3f}'.format(
															perimeter_modulus))
		outfile.write('\nPARAMETER length_threshold = {0:1.12f}'.format(
															length_threshold))
		outfile.write('\nPARAMETER energy_threshold = {0:1.12f}\n'.format(
															energy_threshold))
		outfile.write('\nPARAMETER length_tolerance = {0:1.12f}'.format(
															length_tolerance))
		for boundary_index in range(polygon_sides):
			outfile.write('\nCONSTRAINT {0:d}'.format(boundary_index+1))
			outfile.write('\n\tformula: ')
			outfile.write('{0:1.9f}*x + '.format(x_coeffs[boundary_index]))
			outfile.write('{0:1.9f}*y - '.format(y_coeffs[boundary_index]))
			outfile.write('{0:1.9f} = 0\n'.format(distance))
		outfile.write('\nQUANTITY fixed_edge INFO_ONLY METHOD edge_length\n')
		for face_index, face in enumerate(faces):
			outfile.write('\nMETHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_peri METHOD edge_length\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos METHOD edge_area\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg METHOD edge_area\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_energy ENERGY FUNCTION\n')
			if outer_faces[face_index]:
				outfile.write('\t0\n\n')
			else:
				outfile.write('\t(cell_{0:d}'.format(face_index+1))
				outfile.write('_area_pos.value - ')
				outfile.write('cell_{0:d}'.format(face_index+1))
				outfile.write('_area_neg.value - 1)^2 +\n')
				outfile.write('\t(cell_{0:d}'.format(face_index+1))
				outfile.write('_peri.value - p0_shape_index)^2')
				outfile.write('/r_peri_modulus\n\n')
		outfile.write('\nvertices\n')
		for vertex_index, vertex in enumerate(vertices):
			outfile.write('{0:d}\t{1:1.10f}\t{2:1.10f}'.format(
														vertex_index+1,
														vertex[0],
														vertex[1]))
			if outer_vertices[vertex_index]:
				outfile.write('\tfixed')
			for boundary_index in range(boundary_vertices.shape[1]):
				if boundary_vertices[vertex_index,boundary_index]:
					outfile.write('\tconstraint {0:d}'.format(
														boundary_index+1))
			outfile.write('\n')
		outfile.write('\nedges\n')
		for edge_index, edge in enumerate(edges):
			outfile.write('{0:d}\t{1:d}\t{2:d}'.format(
										edge_index+1, edge[0], edge[1]))
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
		outfile.write('\nconj_grad on') # Much faster convergence!
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
#		outfile.write('ritz(-1, 2*vertex_count);\n')

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Generate initial condition')
	# N = int(sys.argv[1])
	parser.add_argument('-n', '--number',
						nargs = 1,
						default = [96],
						type = int,
						required = False,
						help = 'number of points')
	parser.add_argument('-g', '--polygon',
						nargs = 1,
						default = [12],
						type = int,
						required = False,
						help = 'number of sides of the polygon area')
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = ['../bounded_output/initial_state.fe'],
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
	parser.add_argument('-m', '--make_plot',
						action='store_const',
						const=True, default=False,
						help = 'Output plot of initial Voronoi diagram')
	args = parser.parse_args()
	outfile = Path(args.outfile[0])
	outfile.parent.mkdir(exist_ok = True)
	scale_factor = np.sqrt(args.number[0] * 2 / args.polygon[0] / \
							np.sin(2*np.pi/args.polygon[0]))
	make_initial(N = args.number[0],
				 polygon_sides = args.polygon[0],
				 outfile = outfile,
				 scale_factor = scale_factor,
				 shape_index = args.p0_param[0],
				 perimeter_modulus = args.r_param[0],
				 make_plot = args.make_plot)

################################################################################
# EOF
