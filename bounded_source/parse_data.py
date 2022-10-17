#!/usr/bin/env python3
import numpy as np

################################################################################
#===============================================================================
# parse_data
#===============================================================================
################################################################################

def cast_list(test_list, data_type):
	return list(map(data_type, test_list))

################################################################################

def parse_data (data_filepath):
	vertices = np.empty((0,2), dtype=float)
	edges = np.empty((0,2), dtype=int)
	faces = np.empty(0, dtype=object)
	types = np.empty(0, dtype=int)
	with open(data_filepath) as instream:
		output_type = 'none'
		for line in instream:
			if line == '# vertices\n':
				output_type = 'vertices'
			elif line == '# edges\n':
				output_type = 'edges'
			elif line == '# faces\n':
				output_type = 'faces'
			elif line == '# centroids\n':
				output_type = 'centroids'
			elif output_type == 'vertices' and line != '\n':
				new_vertex = np.array(cast_list(line.split()[1:3], float))
				vertices = np.append(vertices, new_vertex[np.newaxis,:], axis=0)
			elif output_type == 'edges' and line != '\n':
				new_edge = np.array(cast_list(line.split()[1:3], int))
				edges = np.append(edges, new_edge[np.newaxis,:], axis=0)
			elif output_type == 'faces' and line != '\n':
				new_face = cast_list(line.split()[1:-1], int)
				faces = np.append(faces, [None])
				faces[-1] = new_face
				types = np.append(types, int(line.split()[-1]))
	return vertices, edges, faces, types

################################################################################

if __name__ == '__main__':
	pass

################################################################################
# EOF
