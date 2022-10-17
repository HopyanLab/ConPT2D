#!/usr/bin/env /usr/bin/python3
import numpy as np
from scipy.spatial import Voronoi
from pathlib import Path
from functools import reduce
from operator import concat
from random_polygon import random_polygon
from random_polygon import rotate_point
from plot_system import plot_system

################################################################################
#===============================================================================
# make_voronoi.py
#===============================================================================
################################################################################

def make_voronoi (N = 64, polygon_sides = 8,):
	# Generate random points in our polygon.
	points, mirror_points = random_polygon(N, polygon_sides)
	# Generate voronoi diagram for that those points.
	vor = Voronoi(np.vstack([points, mirror_points]))
	# We only want the regions inside the boundary
	region_indices = vor.point_region[:N]
	regions = np.array(vor.regions, dtype=object)[region_indices]
	# Easy to make a backward map. We need to reverse it.
	used_vertices = np.unique(reduce(concat, regions))
	vertices = vor.vertices[used_vertices]
	vertex_map = np.ones(len(vor.vertices), dtype=int) * -1
	vertex_map[used_vertices] = np.arange(len(used_vertices), dtype=int)
	# Use this map to renumber regions.
	remapped_regions = np.empty(len(regions), dtype=object)
	for index,region in enumerate(regions):
		remapped_regions[index] = vertex_map[region]
	regions = remapped_regions
	# make sure they are all counterclockwise
	for region_index, region in enumerate(regions):
		if is_clockwise(vertices, region):
			regions[region_index] = region[::-1]
	# Build edge and face lists in Surface Evolver notation.
	edges = np.empty((0,2), dtype = int)
	boundary_edges = np.empty(0, dtype = bool)
	faces = np.empty(len(regions), dtype = object)
	for face_index, region in enumerate(regions):
		face = np.zeros(len(region), dtype = int)
		for i in range(len(region)):
			edge = np.array([region[i-1]+1, region[i]+1])
			orientation = 1
			if region[i] < region[i-1]:
				edge = edge[::-1]
				orientation = -1
			edge_search = np.where((edges[:,0] == edge[0]) & \
								   (edges[:,1] == edge[1]))[0]
			if len(edge_search) == 0:
				edge_index = len(edges)
				edges = np.append(edges, edge[np.newaxis,:], axis=0)
				boundary_edges = np.append(boundary_edges, True)
			else:
				edge_index = edge_search[0]
				boundary_edges[edge_index] = not(boundary_edges[edge_index])
			face[i] = orientation * (edge_index + 1)
	#	faces = np.append(faces, [None])
		faces[face_index] = face
	# make lists of vertices on each boundary
	boundary_vertices = np.zeros((vertices.shape[0], polygon_sides),
									dtype = bool)
	slope = -np.sin(2*np.pi/polygon_sides)/(1-np.cos(2*np.pi/polygon_sides))
	for boundary_edge in edges[boundary_edges]:
		for index in range(polygon_sides):
			x1, y1 = rotate_point(vertices[boundary_edge[0]-1,0],
								  vertices[boundary_edge[0]-1,1],
									-2*np.pi*(index+1/2)/polygon_sides)
			x2, y2 = rotate_point(vertices[boundary_edge[1]-1,0],
								  vertices[boundary_edge[1]-1,1],
									-2*np.pi*(index+1/2)/polygon_sides)
			if x1 > 0 and np.abs(x1-x2) < 1e-8:
				boundary_vertices[boundary_edge[0]-1,index] = True
				boundary_vertices[boundary_edge[1]-1,index] = True
				break
	boundary_edges = np.zeros((len(edges),polygon_sides), dtype = bool)
	for index in range(polygon_sides):
		boundary_edges[:,index] = np.logical_and(
				boundary_vertices[edges[:,0]-1,index],
				boundary_vertices[edges[:,1]-1,index])
	return vertices, edges, faces, \
		   boundary_vertices, boundary_edges, \
		   points
	# edges and faces are in surface evolver notation!

################################################################################

def is_clockwise (vertices, region):
	points = vertices[region]
	edge_sum = 0
	for index, point in enumerate(points):
		edge_sum += (point[0] - points[index-1,0]) * \
					(point[1] + points[index-1,1])
	if edge_sum > 0:
		return True
	else:
		return False

################################################################################

if __name__ == '__main__':
	vertices, edges, faces, boundary_vertices, boundary_edges, points = \
		make_voronoi(64,8)
	from matplotlib import pyplot as plt
#	plt.rcParams["axes.prop_cycle"] = \
#						plt.cycler("color", plt.cm.hsv(
#											np.linspace(0,1,len(faces))))
	fig = plt.figure()
	ax = plt.axes()
	ax = plot_system(ax, vertices, edges, faces, boundary_edges, points)
#	for boundary_index in range(boundary_vertices.shape[1]):
#		for point_index in range(boundary_vertices.shape[0]):
#			if boundary_vertices[point_index,boundary_index]:
#				ax.text(vertices[point_index,0],vertices[point_index,1],
#							str(boundary_index))
	plt.show()

################################################################################
# EOF
