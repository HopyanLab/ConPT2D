#!/usr/bin/env /usr/bin/python3
import numpy as np
from scipy.spatial import Voronoi
from scipy.integrate import quad, dblquad
from scipy.optimize import fsolve
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
				  polygon_sides = 12,
				  outfile = Path('../hyperbolic_output/initial_state.fe'),
				  polygon_area = 2*np.pi,
				  shape_index = 3.8,
				  perimeter_modulus = 1.0,
				  length_threshold = 1.e-1,
				  energy_threshold = 1.e-6,
				  make_plot = False):
	def integrand(y,x):
		return x/(x**2+y**2)*(1/np.sqrt(1-x**2-y**2)-1)
	def area(polygon_sides = 4, scale = 1):
		x_0 = scale*np.cos(np.pi/polygon_sides)
		y_0 = scale*np.sin(np.pi/polygon_sides)
		return 2*polygon_sides*quad(integrand, 0, y_0, args=(x_0))[0]
	def equation(scale, polygon_sides = 4, polygon_area = 4*np.pi):
		return area(polygon_sides, scale) - polygon_area
	scale_factor = fsolve(equation, 0.99, args=(polygon_sides, polygon_area))[0]
	#print('Max Radius: {0:1.5f}'.format(scale_factor))
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
	#	outfile.write('\nKLEIN_METRIC')
	#	outfile.write('\nMETRIC')
	#	outfile.write('\n(1-y^2)/(1-x^2-y^2)^2')
	#	outfile.write('\tx*y/(1-x^2-y^2)^2')
	#	outfile.write('\nx*y/(1-x^2-y^2)^2')
	#	outfile.write('\t(1-x^2)/(1-x^2-y^2)^2\n')
		outfile.write('\nSURFACE_DIMENSION 1')
		outfile.write('\nSPACE_DIMENSION   2')
		outfile.write('\nDEFINE facet ATTRIBUTE type integer\n')
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
		outfile.write('\nPARAMETER correction = {0:1.12f}'.format(
													np.sqrt(N/polygon_area)))
		outfile.write('\nPARAMETER polygon_area = {0:1.12f}'.format(
													polygon_area))
		for boundary_index in range(polygon_sides):
			outfile.write('\nCONSTRAINT {0:d}'.format(boundary_index+1))
			outfile.write('\n\tformula: ')
			outfile.write('{0:1.9f}*x + '.format(x_coeffs[boundary_index]))
			outfile.write('{0:1.9f}*y - '.format(y_coeffs[boundary_index]))
			outfile.write('{0:1.9f} = 0\n'.format(distance))
		outfile.write('\nQUANTITY fixed_edge INFO_ONLY METHOD')
		outfile.write(' edge_length\n')
		outfile.write('\nMETHOD_INSTANCE line_tension METHOD')
		outfile.write(' edge_general_integral\n')
		outfile.write('\tscalar_integrand:\n')
		outfile.write('\t\tsqrt(x3^2 + x4^2 - (x1*x4 - x2*x3)^2) *')
		outfile.write('\t\t(1/(1-x1^2-x2^2))\n')
		outfile.write('QUANTITY push_appart ENERGY FUNCTION\n')
		outfile.write('\t-1*line_tension.value*correction\n')
		for face_index, face in enumerate(faces):
			outfile.write('\nMETHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_peri METHOD edge_general_integral\n')
			outfile.write('\tscalar_integrand:\n')
			outfile.write('\t\tsqrt(x3^2 + x4^2 - (x1*x4 - x2*x3)^2) *\n')
			outfile.write('\t\t(1/(1-x1^2-x2^2)) * correction\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_pos METHOD edge_general_integral\n')
			outfile.write('\tscalar_integrand:\n')
			outfile.write('\t\t(x1*x4 - x2*x3) * (1/(x1^2+x2^2)) *\n')
			outfile.write('\t\t(1/sqrt(1-x1^2-x2^2) - 1) * correction^2\n')
			outfile.write('METHOD_INSTANCE cell_{0:d}'.format(face_index+1))
			outfile.write('_area_neg METHOD edge_general_integral\n')
			outfile.write('\tscalar_integrand:\n')
			outfile.write('\t\t(x1*x4 - x2*x3) * (1/(x1^2+x2^2)) *\n')
			outfile.write('\t\t(1/sqrt(1-x1^2-x2^2) - 1) * correction^2\n')
			outfile.write('QUANTITY cell_{0:d}'.format(face_index+1))
			outfile.write('_energy ENERGY FUNCTION\n')
			if outer_faces[face_index]:
				outfile.write('\t0\n')
			else:
				outfile.write('\t(cell_{0:d}'.format(face_index+1))
				outfile.write('_area_pos.value - ')
				outfile.write('cell_{0:d}'.format(face_index+1))
				outfile.write('_area_neg.value - 1)^2 +\n')
				outfile.write('\t(cell_{0:d}'.format(face_index+1))
				outfile.write('_peri.value - p0_shape_index)^2')
				outfile.write('/r_peri_modulus\n')
		for edge_index, edge in enumerate(edges):
			outfile.write('\nQUANTITY edge_{0:d}'.format(edge_index+1))
			outfile.write('_length INFO_ONLY METHOD')
			outfile.write(' edge_general_integral\n')
			outfile.write('\tscalar_integrand: ')
			outfile.write('sqrt(x3^2 + x4^2 - (x1*x4 - x2*x3)^2) *\n')
			outfile.write('\t\t(1/(1-x1^2-x2^2)) * correction\n')
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
			outfile.write('\tedge_{0:d}_length'.format(edge_index+1))
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
		outfile.write('\nset scale 1e-6')
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
	#	outfile.write('ritz(-1000, 2*vertex_count);\n')

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Generate initial condition')
	# N = int(sys.argv[1])
	parser.add_argument('-n', '--number',
						nargs = 1,
						default = [128],
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
						default = ['../hyperbolic_output/initial_state.fe'],
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
	parser.add_argument('-a', '--area',
						nargs = 1,
						default = [2*np.pi],
						type = float,
						required = False,
						help = 'area of polygon in klein disk')
	parser.add_argument('-m', '--make_plot',
						action='store_const',
						const=True, default=False,
						help = 'Output plot of initial Voronoi diagram')
	args = parser.parse_args()
	outfile = Path(args.outfile[0])
	outfile.parent.mkdir(exist_ok = True)
	make_initial(N = args.number[0],
				 polygon_sides = args.polygon[0],
				 outfile = outfile,
				 polygon_area = args.area[0],
				 shape_index = args.p0_param[0],
				 perimeter_modulus = args.r_param[0],
				 make_plot = args.make_plot)

################################################################################
# EOF
