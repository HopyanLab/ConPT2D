#!/usr/bin/env /usr/bin/python3
import numpy as np

################################################################################
#===============================================================================
# periodic.py
#===============================================================================
################################################################################

# Difference with respect to periodic boundaries
def periodic_diff(point1, point2, L = np.array([1., 1.])):
	return ((point1 - point2 + L/2.) % L) - L/2.

################################################################################

# Get vertices inside original bounds
#   assumes box length is 1
def get_vertices(tile_vertices, box_size = 1.):
	vertices = []
	index_map = {}
	for i,(x,y) in enumerate(tile_vertices):
		if x >= 0 and x <= box_size:
			if y >= 0 and y <= box_size:
				curr_indx = len(vertices)
				index_map[i] = curr_indx
				vertices.append((x,y))
	return np.array(vertices), index_map

################################################################################

# Get edges with respect to periodic boundaries
#   vertices and edges are from voronoi with tiled coordinates
def get_edges(vertices, edges, index_map, box_size = 1.):
	reduced_edges = []
	for edge in edges:
		point1 = edge[0]
		point2 = edge[1]
		# Case 1: point1 and point2 in index map
		#   Append indices for vertex list wrt periodic bounds
		if point1 in index_map and point2 in index_map:
			reduced_edges.append((index_map[point1],index_map[point2]))
		# Case 2: neither in index map
		#   Do nothing (because this is an edge totally outside the region)
		# Case 3: point1 or point2 in index map, but not both.
		#   This edge wraps the boundary and so need to be associated
		#   with 
		if point1 in index_map and point2 not in index_map:
			# find the vertex that it wraps around to
			# vin = vertex in plane
			vin = np.array(vertices[point1])
			# vout = vertex out of plane
			vout = np.array(vertices[point2])
			if vout[0] < 0:
				vout[0] = box_size + vout[0]
			if vout[0] > box_size:
				vout[0] = vout[0] - box_size
			if vout[1] < 0:
				vout[1] = box_size + vout[1]
			if vout[1] > box_size:
				vout[1] = vout[1] - box_size
			# # # find index of this vertex 
			for key in index_map:
				if abs(vertices[key][0] - vout[0]) < 10**-6 and \
							not ((index_map[key], index_map[point1]) in \
														reduced_edges):
					reduced_edges.append((index_map[point1], index_map[key]))
		if point2 in index_map and point1 not in index_map:
			# find the vertex that it wraps around to
			# vin = vertex in plane
			vin = np.array(vertices[point2])
			# vout = vertex out of plane
			vout = np.array(vertices[point1])
			if vout[0] < 0:
				vout[0] = box_size + vout[0]
			if vout[0] > box_size:
				vout[0] = vout[0] - box_size
			if vout[1] < 0:
				vout[1] = box_size + vout[1]
			if vout[1] > box_size:
				vout[1] = vout[1] - box_size
			# # # find index of this vertex 
			for key in index_map:
				if abs(vertices[key][0] - vout[0]) < 10**-6 and \
							not ((index_map[key], index_map[point2]) in \
														reduced_edges):
					reduced_edges.append((index_map[point2], index_map[key]))
	return reduced_edges

################################################################################

# Add indices mapping to new vertices outside of bounds
def get_new_index_map(vertices, v, index_map, box_size = 1.):
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
		for j,(x2,y2) in enumerate(v):
			if abs(x1 - x2) < 10**-6:
				if abs(y1 - y2) < 10**-6:
					index_map[i] = j
	return index_map

################################################################################

# Want to check orientation. Is it counterclock-wise (right-handed)?
def check_counter(vertices, face, box_size = 1.):
	L = np.array([box_size, box_size])
	points = vertices[face]
	point0 = points[0]
	sumEdges = 0
	for i,point1 in enumerate(points):
		point1 = point0 + periodic_diff(point1, point0, L)
		if i == len(points) - 1:
			point2 = points[0]
		else:
			point2 = points[i+1]
		point2 = point0 + periodic_diff(point2, point0, L)
		sumEdges += (point2[0] - point1[0]) * (point2[1] + point1[1])
	if sumEdges < 0:
		return True
	if sumEdges >= 0:
		return False
	return True

################################################################################

# Check if lists are the same order, but not necessarily indexed the same
def same_list(list1, list2):
	list1_index = list1[0]
	list2_index = -1
	for i,check_index in enumerate(list2):
		if list1_index == check_index:
			list2_index = i
	if list2_index == -1:
		return False
	if list2_index != -1:
		i2 = list2_index
		for i1 in range(0,len(list1)):
			if list1[i1] != list2[i2]:
				return False
			i2 += 1
			if i2 >= len(list2):
				i2 = 0
		return True

################################################################################

# Get faces
def get_faces(regions, edges, index_map, vertices):
	faces = []
	for region in regions:
		count = 0
		face = []
		for index in region:
			if index in index_map:
				count += 1
				face.append(index_map[index])
		if count == len(region) and count != 0:
			faces.append(face)
	for i,face in enumerate(faces):
		cc = check_counter(vertices, face)
		if cc == False:
			rev_face = []
			for index in reversed(face):
				rev_face.append(index)
			faces[i] = rev_face
	# Remove duplicates
	new_faces = []
	for face in faces:
		add = True
		for face2 in new_faces:
			if len(face) == len(face2):
				# Check if they are the same ordered list
				#same = same_list(face, face2)
				# Was maybe overthinking things.
				# This is probably good enough.
				same = (set(face) == set(face2))
				if same == True:
					add = False
		if add == True:
			new_faces.append(face)
#	return new_faces
	faces = new_faces
	edges = np.array(edges)
	new_faces = []
	# Switch from index numbering vertices to edges
	for face in faces:
		new_face = []
		# Last vertex to first vertex defines the last edge,
		#   so we want that first vertex at the end also.
		face.append(face[0])
		for vi, vertex in enumerate(face[:-1]):
			for ei, edge in enumerate(edges):
				if (edge == [vertex, face[vi+1]]).all():
					new_face.append(ei+1)
				if (edge == [face[vi+1], vertex]).all():
					new_face.append(-ei-1)
		new_faces.append(new_face)
	return new_faces

################################################################################

if __name__ == '__main__':
	pass

################################################################################
# EOF
