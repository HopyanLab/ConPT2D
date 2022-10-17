#!/usr/bin/env /usr/bin/python3
import numpy as np
from copy import deepcopy
from scipy.spatial import SphericalVoronoi
from pathlib import Path
import argparse

################################################################################
#===============================================================================
# make_voronoi.py
#===============================================================================
################################################################################

def random_points (number = 1, radius = 1.0):
	np.random.seed()
	phi = np.random.uniform(low = 0., high = 2*np.pi, size = number)
	theta = np.arccos(np.random.uniform(low = -1., high = 1., size = number))
	return np.vstack([radius * np.sin(theta) * np.cos(phi),
					  radius * np.sin(theta) * np.sin(phi),
					  radius * np.cos(theta)]).T

################################################################################

def make_voronoi (N = 64):
	# Radius of sphere is chosen so that surface area is number of points.
	R = np.sqrt(N/4/np.pi)
	# Generate random points
	points = random_points(number = N, radius = R)
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
	return vertices, edges, faces, regions
