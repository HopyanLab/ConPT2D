#!/usr/bin/env /usr/bin/python3
import numpy as np
from mayavi import mlab

################################################################################
#===============================================================================
# plot_spherical_system.py
#===============================================================================
################################################################################

def subdivide_triangular_surface (vertices, faces):
	new_vertices = vertices.copy()
	new_faces = np.zeros((4*len(faces), 3), dtype = int)
	new_index_map = np.zeros((0,3), dtype = int)
	#
	def split_edge(vertex_1_index, vertex_2_index):
		nonlocal new_vertices
		nonlocal new_faces
		nonlocal new_index_map
		# 'new_index_map' is stuctured as
		#  [[vertex_1_index, vertex_2_index, new_index], ...]
		#  with vertex_1_index < vertex_2_index
		if vertex_2_index < vertex_1_index:
			vertex_1_index, vertex_2_index = vertex_2_index, vertex_1_index
		try:
			new_index_entry = new_index_map[
				 np.logical_and(new_index_map[:,0] == vertex_1_index,
								new_index_map[:,1] == vertex_2_index)]
		except:
			new_index_entry = np.zeros((0,3), dtype = int)
		if len(new_index_entry) == 0:
			# Need to make new vertex, add to list, then return index.
			new_vertex_index = len(new_vertices)
			new_vertex_position = (vertices[vertex_1_index,:] + \
								   vertices[vertex_2_index,:]) / 2.
			# Need to put new index on the sphere too.
			new_vertex_position /= np.linalg.norm(new_vertex_position)
			new_vertices = np.append(new_vertices,
									[new_vertex_position], axis=0)
			new_index_map = np.append(new_index_map,
									np.array([vertex_1_index, vertex_2_index,
									new_vertex_index]))
			return new_vertex_index
		else:
			# Don't need to make new vertex just return index.
			return new_index_entry[0,2]
	#
	for face in faces:
		v1_index, v2_index, v3_index = face
		v4_index = split_edge(v1_index, v2_index)
		v5_index = split_edge(v2_index, v3_index)
		v6_index = split_edge(v3_index, v1_index)
		new_face = np.array([[v1_index, v4_index, v6_index],
							 [v4_index, v2_index, v5_index],
							 [v6_index, v5_index, v3_index],
							 [v4_index, v5_index, v6_index]])
		new_faces = np.append(new_faces, new_face, axis=0)
	return new_vertices, new_faces

################################################################################

def spherical_plot (vertices, faces):
	mlab.figure(bgcolor=(1.,1.,1.), size=(1000,1000))
	for face in faces:
		center_point = np.sum(vertices[face], axis=0)
		center_point /= np.linalg.norm(center_point)
		points = np.append(vertices[face], [center_point], axis=0)
		triangles = np.zeros((len(points)-1,3), dtype = int)
		triangles[-1,:] = np.array([len(points)-2, 0, len(points)-1])
		for i in np.arange(0, len(points)-2):
			triangles[i,:] = np.array([i, i+1, len(points)-1])
		points, triangles = subdivide_triangular_surface(points, triangles)
		mlab.triangular_mesh(points[:,0], points[:,1], points[:,2],
							 triangles, color = (0,1,0),
							 representation = 'surface')
	for face in faces:
		points = vertices[face]
		number = 2 * len(face)
		points = np.zeros((number,3), dtype=float)
		points[:number:2,:] = vertices[face]
		points[1:number:2,:] = vertices[face] - (vertices[face] - \
					np.roll(vertices[face], -1, axis = 0)) / 2
		points[1:number:2,:] = (points[1:number:2,:].T / \
					np.linalg.norm(points[1:number:2,:], axis = 1)).T
		for i in range(0, number):
			mlab.plot3d([points[i, 0],points[(i+1)%number, 0]],
						[points[i, 1],points[(i+1)%number, 1]],
						[points[i, 2],points[(i+1)%number, 2]],
						representation = 'surface',
						color = (0,0,1), tube_radius = 0.01, opacity = 1.0)
		mlab.points3d(points[:,0], points[:,1], points[:,2],
						color = (0,0,1), scale_factor = 0.02)
	mlab.points3d(vertices[:,0], vertices[:,1], vertices[:,2],
					color = (1,0,0), scale_factor = 0.05)
	mlab.show()

################################################################################

def plot_system ():
	#TODO

################################################################################

if __name__ == '__main__':
	# Filenames can be passed as arguements so
	#  we use argparse module to parse them.
	parser = argparse.ArgumentParser(
							description = 'Plot spherical Voronoi diagram')
	# outfile = sys.argv[1]
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = ['./output_spherical/plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('-v', '--vertfile',
						nargs = 1,
						default = ['./output_spherical/vertices.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	parser.add_argument('-e', '--edgefile',
						nargs = 1,
						default = ['./output_spherical/edges.txt'],
						type = str,
						required = False,
						help = 'file to put edge data in')
	parser.add_argument('-f', '--facefile',
						nargs = 1,
						default = ['./output_spherical/faces.txt'],
						type = str,
						required = False,
						help = 'file to put face data in')
	args = parser.parse_args()
	plot_system(vertfile = args.vertfile[0],
				edgefile = args.edgefile[0],
				facefile = args.facefile[0],
				savefile = args.outfile[0]
				)

################################################################################
# EOF
