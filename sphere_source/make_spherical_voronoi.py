#!/usr/bin/env /usr/bin/python3
import numpy as np
from copy import deepcopy
from scipy.spatial import SphericalVoronoi
from pathlib import Path
import argparse

################################################################################
#===============================================================================
# make_spherical_voronoi.py
#===============================================================================
################################################################################

def random_spherical_points (number = 1, radius = 1.0):
	phi = np.random.uniform(low = 0., high = 2*np.pi, size = number)
	theta = np.arccos(np.random.uniform(low = -1., high = 1., size = number))
	return np.vstack([radius * np.sin(theta) * np.cos(phi),
					  radius * np.sin(theta) * np.sin(phi),
					  radius * np.cos(theta)]).T

################################################################################

def make_spherical_voronoi (N = 64,
						vertfile = Path('../sphere_output/vertices.txt'),
						edgefile = Path('../sphere_output/edges.txt'),
						facefile = Path('../sphere_output/faces.txt'),
						suevfile = Path('../sphere_output/initial_state.fe'),
						shape_index = 3.0,
						perimeter_modulus = 1.0,
						length_threshold = 1.e-1,
						energy_threshold = 1.e-8):
	# Radius of sphere is chosen so that surface area is number of points.
	R = np.sqrt(N/4/np.pi)
	# Generate random points (Poisson Process)
	np.random.seed()
	points = random_spherical_points(number = N, radius = R)
	# Generate a Spherical Voronoi diagram based on those points.
	vor = SphericalVoronoi(points, radius = R, center = [0,0,0])
	# Sort the vertices in clockwise or counterclockwise order.
	vor.sort_vertices_of_regions()
	vertices = vor.vertices
	regions = vor.regions
	# Reverse the clockwise regions.
	for i,region in enumerate(regions):
		edge1 = vertices[region[1]] - vertices[region[0]]
		edge2 = vertices[region[2]] - vertices[region[1]]
		cross_product = np.cross(edge1, edge2)
		if np.dot(cross_product, vertices[region[1]]) < 0:
			regions[i] = list(np.flip(region))
	# Write vertices
	fmt = '%d', '%f', '%f', '%f'
	np.savetxt(vertfile, np.insert(vertices, 0 ,
							np.arange(1, len(vertices)+1, 1),
								axis = 1), fmt=fmt, delimiter='\t')
	# Build arrays for edges and faces
	edges = np.zeros((0,2), dtype = int)
	faces = deepcopy(regions)
	for region_index, region in enumerate(regions):
		for index in np.arange(len(region)):
			new_edge = np.sort(np.array([region[index],
							region[(index+1)%len(region)]]))
			edge_direction = 1 if new_edge[0] == region[index] else -1
			edge_check = np.equal([new_edge], edges).all(axis=1)
			if edge_check.any():
				edge_index = np.where(edge_check)[0][0]
			else:
				edge_index = len(edges)
				edges = np.append(edges, [new_edge], axis=0)
			faces[region_index][index] = edge_direction * (edge_index + 1)
	# Write edges
	fmt = '%d', '%d', '%d'
	np.savetxt(edgefile, np.insert(edges, 0 ,
					np.arange(1, len(edges)+1, 1), axis = 1),
				delimiter='\t', fmt=fmt)
	# Write faces
	with open(facefile, 'w') as outstream:
		for face_index, face in enumerate(faces):
			outstream.write('{0:d}'.format(face_index+1))
			for entry in face:
				outstream.write('\t{0:d}'.format(entry))
			outstream.write('\n')
	# Write the Surface Evolver initial state
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
		outfile.write('\nPARAMETER length_threshold = {0:1.9f}'.format(
															length_threshold))
		outfile.write('\nPARAMETER energy_threshold = {0:1.9f}'.format(
															energy_threshold))
		outfile.write('\nPARAMETER radius           = {0:f}\n'.format(R))
		outfile.write('\nCONSTRAINT_TOLERANCE 1.e-9')
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
					'spherical_procedures.fe','r') as commandfile:
			for line in commandfile.read():
				outfile.write(line)
		# Graphics command part. Uncomment for testing.
	#	outfile.write('\nshow\nq\n')
		# Try to find energy minimum. Uncomment for testing.
	#	outfile.write('\nrelax_system(10000);')
	#	outfile.write('\nJ;\n0.01\nrelax_system(100);\nJ;')
	#	outfile.write('\nrelax_system(1000);\n')
		# Graphics command part. Uncomment for testing.
	#	outfile.write('\nshow\nq\n')

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
	parser.add_argument('-v', '--vertfile',
						nargs = 1,
						default = ['../sphere_output/vertices.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	parser.add_argument('-e', '--edgefile',
						nargs = 1,
						default = ['../sphere_output/edges.txt'],
						type = str,
						required = False,
						help = 'file to put edge data in')
	parser.add_argument('-f', '--facefile',
						nargs = 1,
						default = ['../sphere_output/faces.txt'],
						type = str,
						required = False,
						help = 'file to put face data in')
	parser.add_argument('-s', '--suevfile',
						nargs = 1,
						default = ['../sphere_output/initial_state.fe'],
						type = str,
						required = False,
						help = 'file to put surface evolver script in')
	parser.add_argument('-p', '--p0_param',
						nargs = 1,
						default = [3.6],
						type = float,
						required = False,
						help = 'shape index parameter')
	parser.add_argument('-r', '--r_param',
						nargs = 1,
						default = [1.0],
						type = float,
						required = False,
						help = 'inverse perimeter modulus')
	args = parser.parse_args()
	vertfile = Path(args.vertfile[0])
	vertfile.parent.mkdir(exist_ok = True)
	edgefile = Path(args.edgefile[0])
	edgefile.parent.mkdir(exist_ok = True)
	facefile = Path(args.facefile[0])
	facefile.parent.mkdir(exist_ok = True)
	suevfile = Path(args.suevfile[0])
	suevfile.parent.mkdir(exist_ok = True)
	make_spherical_voronoi(N = args.number[0],
						vertfile = vertfile,
						edgefile = edgefile,
						facefile = facefile,
						suevfile = suevfile,
						shape_index = args.p0_param[0],
						perimeter_modulus = args.r_param[0])

################################################################################
# EOF
