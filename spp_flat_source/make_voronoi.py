#!/usr/bin/env /usr/bin/python3
import numpy as np
from scipy.spatial import Voronoi
import argparse
from pathlib import Path

################################################################################
#===============================================================================
# make_voronoi.py
#===============================================================================
################################################################################

def make_voronoi (N = 64, box_size = 1.):
	# Generate random points (Poisson Process)
	np.random.seed()
	points = np.random.uniform(0, box_size, size=(N,2))
	# Tile points around the region we care about
	point_tile = tile_points(points)
	# Generate voronoi diagram for that tiled space
	vor = Voronoi(point_tile)
	# Get just the vertices we care about (undo tiling) and
	#   create an index map from old vertex index to new vertex index
	tile_vertices = vor.vertices
	vertices, index_map = get_vertices(tile_vertices)
	# Get edges
	tile_edges = vor.ridge_vertices
	edges, wrappings = get_edges(tile_vertices, tile_edges, index_map)
	# Get faces
	# New index map finds the points that are just outside the box
	#   but will correspond to the indices of faces built in tile
	index_map = get_new_index_map(tile_vertices, vertices, index_map)
	# Faces in tiled vertices
	regions = vor.regions
	faces = get_faces(regions, edges, wrappings, index_map, vertices)
	edges += 1
	return vertices, edges, wrappings, faces

################################################################################

# Tile points around original region
def tile_points(points, box_size = 1.):
	#
	N = len(points[:,0])
	#
	# 9 tiles surrounding individual coordinates
	point_tile = np.zeros((9 * N, 2))
	#
	# original coordinates
	point_tile[:N] = points
	#
	# upper left
	point_tile[N:2*N, 0] = points[:,0] - box_size
	point_tile[N:2*N, 1] = points[:,1] + box_size
	#
	# upper centre
	point_tile[2*N:3*N, 0] = points[:,0]
	point_tile[2*N:3*N, 1] = points[:,1] + box_size
	#
	# upper right
	point_tile[3*N:4*N, 0] = points[:,0] + box_size
	point_tile[3*N:4*N, 1] = points[:,1] + box_size
	#
	# centre right
	point_tile[4*N:5*N, 0] = points[:,0] + box_size
	point_tile[4*N:5*N, 1] = points[:,1]
	#
	# lower right
	point_tile[5*N:6*N, 0] = points[:,0] + box_size
	point_tile[5*N:6*N, 1] = points[:,1] - box_size
	#
	# lower centre
	point_tile[6*N:7*N, 0] = points[:,0]  
	point_tile[6*N:7*N, 1] = points[:,1] - box_size
	#
	# lower left
	point_tile[7*N:8*N,0] = points[:,0] - box_size
	point_tile[7*N:8*N,1] = points[:,1] - box_size
	#
	# centre left
	point_tile[8*N:,0] = points[:,0] - box_size
	point_tile[8*N:,1] = points[:,1]
	#
	return point_tile

################################################################################

# Get vertices inside original bounds
def get_vertices(tile_vertices, box_size = 1.):
	vertices = []
	index_map = {}
	for i,(x,y) in enumerate(tile_vertices):
		if x >= 0 and x <= box_size:
			if y >= 0 and y <= box_size:
				current_index = len(vertices)
				index_map[i] = current_index
				vertices.append((x,y))
	return np.array(vertices, dtype = float), index_map

################################################################################

# Get edges with respect to periodic boundaries
#   vertices and edges are from voronoi with tiled coordinates
def get_edges(vertices, edges, index_map, box_size = 1.):
	reduced_edges = []
	wrappings = []
	for edge in edges:
		point_1 = edge[0]
		point_2 = edge[1]
		# Case 1: point_1 and point_2 in index map
		#   Append indices for vertex list wrt periodic bounds
		if point_1 in index_map and point_2 in index_map:
			reduced_edges.append((index_map[point_1],index_map[point_2]))
			wrappings.append((0,0))
			continue
		# Case 2: neither in index map
		#   Do nothing (because this is an edge totally outside the region)
		if point_1 not in index_map and point_2 not in index_map:
			continue
		# Case 3: point_1 or point_2 in index map, but not both.
		#   This edge wraps the boundary and so need to be associated
		#   with point in the region
		if point_1 in index_map and point_2 not in index_map:
			# find the vertex that it wraps around to
			# vin = vertex in plane
			vin = np.array(vertices[point_1])
			# vout = vertex out of plane
			vout = np.array(vertices[point_2])
			wrap_x = 0
			wrap_y = 0
			if vout[0] < 0:
				vout[0] = box_size + vout[0]
				wrap_x = -1
			if vout[0] > box_size:
				vout[0] = vout[0] - box_size
				wrap_x = 1
			if vout[1] < 0:
				vout[1] = box_size + vout[1]
				wrap_y = -1
			if vout[1] > box_size:
				vout[1] = vout[1] - box_size
				wrap_y = 1
			# # # find index of this vertex 
			for key in index_map:
				if abs(vertices[key][0] - vout[0]) < 10**-9 and \
				   abs(vertices[key][1] - vout[1]) < 10**-9 and \
				   (index_map[key], index_map[point_1]) not in \
													reduced_edges:
					reduced_edges.append((index_map[point_1], index_map[key]))
					wrappings.append((wrap_x, wrap_y))
					break
			continue
		if point_2 in index_map and point_1 not in index_map:
			# find the vertex that it wraps around to
			# vin = vertex in plane
			vin = np.array(vertices[point_2])
			# vout = vertex out of plane
			vout = np.array(vertices[point_1])
			wrap_x = 0
			wrap_y = 0
			if vout[0] < 0:
				vout[0] = box_size + vout[0]
				wrap_x = -1
			if vout[0] > box_size:
				vout[0] = vout[0] - box_size
				wrap_x = 1
			if vout[1] < 0:
				vout[1] = box_size + vout[1]
				wrap_y = -1
			if vout[1] > box_size:
				vout[1] = vout[1] - box_size
				wrap_y = 1
			# # # find index of this vertex 
			for key in index_map:
				if abs(vertices[key][0] - vout[0]) < 10**-9 and \
				   abs(vertices[key][1] - vout[1]) < 10**-9 and \
				   (index_map[key], index_map[point_2]) not in \
													reduced_edges:
					reduced_edges.append((index_map[point_2], index_map[key]))
					wrappings.append((wrap_x, wrap_y))
					break
			continue
	return np.array(reduced_edges, dtype = int), \
		   np.array(wrappings, dtype = int)

################################################################################

# Add indices mapping to new vertices outside of bounds
def get_new_index_map(vertices, reduced_vertices, index_map, box_size = 1.):
	for i,(x,y) in enumerate(vertices):
		x1 = x
		y1 = y
		if x < 0 and x > -box_size:
			x1 = x + box_size
		if x > box_size and x < 2*box_size:
			x1 = x - box_size
		if y < 0 and y > -box_size:
			y1 = y + box_size
		if y > box_size and y < 2*box_size:
			y1 = y - box_size
		# look up new x,y in list
		for j,(x2,y2) in enumerate(reduced_vertices):
			if abs(x1 - x2) < 10**-9 and \
			   abs(y1 - y2) < 10**-9:
					index_map[i] = j
	return index_map

################################################################################

# Want to check orientation. Is it counterclock-wise (right-handed)?
def check_counter(vertices, edges, wrappings, face, box_size = 1.):
	sumEdges = 0
	wrapped_x, wrapped_y = 0,0
	for edge_index in face:
		if edge_index > 0:
			point_1, point_2 = vertices[edges[edge_index-1]]
			wrap_x, wrap_y = wrappings[edge_index-1]
		else:
			point_2, point_1 = vertices[edges[-edge_index-1]]
			wrap_x, wrap_y = -wrappings[-edge_index-1]
		point_1[0] += wrapped_x * box_size
		point_1[1] += wrapped_y * box_size
		point_2[0] += wrapped_x * box_size
		point_2[1] += wrapped_y * box_size
		if wrap_x == 1:
			point_2[0] += box_size
			wrapped_x += 1
		elif wrap_x == -1:
			point_2[0] -= box_size
			wrapped_x -= 1
		if wrap_y == 1:
			point_2[1] += box_size
			wrapped_y += 1
		elif wrap_y == -1:
			point_2[1] -= box_size
			wrapped_y -= 1
		sumEdges += (point_2[0] - point_1[0]) * (point_2[1] + point_1[1])
	if sumEdges < 0:
		return True
	if sumEdges >= 0:
		return False
	return True

################################################################################

# Check if lists are the same order, but not necessarily indexed the same
def same_list(list_1, list_2):
	list_1_index = list_1[0]
	list_2_index = -1
	for i,check_index in enumerate(list_2):
		if list_1_index == check_index:
			list_2_index = i
	if list_2_index == -1:
		return False
	if list_2_index != -1:
		i2 = list_2_index
		for i1 in range(0,len(list_1)):
			if list_1[i1] != list_2[i2]:
				return False
			i2 += 1
			if i2 >= len(list_2):
				i2 = 0
		return True

################################################################################

# Get faces in surface evolver notation based on edges list
def get_faces(regions, edges, wrappings, index_map, vertices):
	faces = []
	for region in regions:
		if len(region) == 0:
			continue
		count = 0
		face = []
		region.append(region[0])
		for i in range(len(region)-1):
			if region[i] in index_map and \
			   region[i+1] in index_map:
				edge = (index_map[region[i]], index_map[region[i+1]])
				forward = np.where(np.all(edges == edge, axis=1))[0]
				if len(forward) == 1:
					face.append(forward[0]+1)
					count += 1
				else:
					backward = np.where(np.all(edges == edge[::-1], axis=1))[0]
					if len(backward) == 1:
						face.append(-backward[0]-1)
						count += 1
		if count == 0 or count != len(region)-1:
			continue
		face = np.array(face)
		if not check_counter(vertices, edges, wrappings, face):
			face = -face[::-1]
		add_face = True
		for check_face in faces:
			if len(face) == len(check_face):
				if same_list(face, check_face):
					add_face = False
					break
		if add_face:
			faces.append(face)
	return np.array(faces, dtype = object)

################################################################################

if __name__ == '__main__':
	pass

################################################################################
# EOF
