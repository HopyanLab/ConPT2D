#!/usr/bin/env /usr/bin/python3
import numpy as np
from scipy.spatial import Voronoi
from random_polygon import rotate_point
from make_voronoi import make_voronoi
from plot_system import plot_system

################################################################################
#===============================================================================
# make_bounded_voronoi.py
#===============================================================================
################################################################################

def make_bounded_voronoi (N = 64,
						  polygon_sides = 6,
						  scale_factor = np.sqrt(1/2),
						  outer_length = 1.e-12):
	# generate Voronoi diagram in polygon
	vertices, edges, faces, boundary_vertices, boundary_edges, points = \
		make_voronoi(N, polygon_sides, scale_factor)
	# add new vertices, edges, and faces outside the polygon.
	new_vertices = np.zeros((polygon_sides,2), dtype = float)
	new_edges = np.zeros((2*polygon_sides,2), dtype = int)
	new_faces = np.empty(polygon_sides, dtype = object)
	for index in range(polygon_sides): # revisit length
		# new vertices are just projected out from corners of polygon
		new_vertices[index] = (scale_factor + outer_length) * \
						np.array([np.cos(2*np.pi*index/polygon_sides),
								  np.sin(2*np.pi*index/polygon_sides)])
		# new edges join the corners radially outward to new vertices
		#  and also around the outside joining those new vertices together
		if index > 0:
			radial = np.array([
						np.where(np.logical_and(
								boundary_vertices[:,index-1],
								boundary_vertices[:,index]
										))[0][0] + 1 + polygon_sides,
							index+1], dtype = int)
		else:
			radial = np.array([
						np.where(np.logical_and(
								boundary_vertices[:,polygon_sides-1],
								boundary_vertices[:,0]
										))[0][0] + 1 + polygon_sides,
							1], dtype = int)
		if index < polygon_sides - 1:
			outer = np.array([index+1, index+2], dtype = int)
		else:
			outer = np.array([index+1, 1], dtype = int)
		new_edges[2*index] = radial
		new_edges[2*index+1] = outer
		# new faces are outside each edge of polygon
		edges_on = np.nonzero(boundary_edges[:,index])[0] # np.where works same
		new_face = np.zeros(len(edges_on)+3, dtype = int)
		bound_x, bound_y = rotate_point(vertices[edges[edges_on]-1,0],
									vertices[edges[edges_on]-1,1],
									-2*np.pi*index/polygon_sides)
		sorted_x = np.sort(bound_x, axis=1)
		edge_order = np.argsort(sorted_x[:,0])
		swaped = (bound_x[:,0] != sorted_x[:,0])
		new_face[:-3] = (edges_on[edge_order] + 1 + 2*polygon_sides) * \
							(1-2*swaped[edge_order])
		new_face[-3] = 2*index+1
		new_face[-2] = 2*index+2
		if index < polygon_sides-1:
			new_face[-1] = -2*index-3
		else:
			new_face[-1] = -1
		new_faces[index] = new_face
	# next we renumber the vertices and edges of old arrays
	#  since the new ones will be inserted at the beginning
	edges += polygon_sides
	for face_index in range(len(faces)):
		for edge_index in range(len(faces[face_index])):
			if faces[face_index][edge_index] > 0:
				faces[face_index][edge_index] += 2*polygon_sides
			else:
				faces[face_index][edge_index] -= 2*polygon_sides
	vertices = np.concatenate((new_vertices, vertices), axis=0,
								dtype = float)
	edges = np.concatenate((new_edges, edges), axis=0,
								dtype = int)
	faces = np.concatenate((new_faces, faces), dtype = object)
	boundary_vertices = np.vstack((np.zeros((polygon_sides, polygon_sides),
										dtype = bool),  boundary_vertices))
	boundary_edges = np.vstack((np.zeros((2*polygon_sides, polygon_sides),
										dtype = bool), boundary_edges))
	outer_vertices = np.zeros(vertices.shape[0], dtype = bool)
	outer_vertices[:polygon_sides] = True
	outer_edges = np.zeros(edges.shape[0], dtype = bool)
	outer_edges[:2*polygon_sides] = True
	outer_faces = np.zeros(faces.shape[0], dtype = bool)
	outer_faces[:polygon_sides] = True
	return vertices, edges, faces, \
			boundary_vertices, boundary_edges, \
			outer_vertices, outer_edges, outer_faces


################################################################################

if __name__ == '__main__':
	vertices, edges, faces, boundary_vertices, boundary_edges, \
					outer_points, outer_edges, outer_faces = \
											make_bounded_voronoi(64,6)
	from matplotlib import pyplot as plt
	fig = plt.figure()
	ax = plt.axes()
	ax = plot_system(ax, vertices, edges, faces,
					 boundary_edges, vertices[outer_points])
	plt.show()

################################################################################
# EOF
