#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
import argparse
from pathlib import Path

################################################################################
#===============================================================================
# plot_system.py
#===============================================================================
################################################################################

def plot_system (ax, vertices, edges, faces = None,
				 boundary_edges = None, points = None):
	if faces is not None:
		for face in faces:
			if face is None:
				continue
			face_vertices = np.zeros((len(face),2), dtype = float)
			for index in range(len(face)):
				if face[index] > 0:
					face_vertices[index] = vertices[edges[face[index]-1,0]-1]
				else:
					face_vertices[index] = vertices[edges[-face[index]-1,1]-1]
			ax.fill(face_vertices[:, 0], face_vertices[:, 1],
					alpha = 0.3, #facecolor = 'tab:grey',
					edgecolor = 'tab:grey', linewidth = 3)
	edge_segments = LineCollection(
							np.stack([vertices[edges[:,0]-1],
									  vertices[edges[:,1]-1]],
										axis=1),
									color = 'black',
									linestyles = 'solid',
									linewidth = 1)
	ax.add_collection(edge_segments)
	ax.plot(vertices[:,0], vertices[:,1],
				marker = '.', color = 'black',
				linestyle = '', markersize = 7)
	ax.plot(vertices[:,0], vertices[:,1],
				marker = '.', color = 'tab:red',
				linestyle = '', markersize = 3)
	if points is not None:
		ax.plot(points[:,0], points[:,1],
					marker = '.', color = 'tab:blue',
					linestyle = '', markersize = 3)
	if boundary_edges is not None:
		if len(boundary_edges.shape) > 1:
			boundary_edges = np.any(boundary_edges, axis=1)
		bound_segments = LineCollection(
							np.stack([vertices[edges[boundary_edges,0]-1],
									  vertices[edges[boundary_edges,1]-1]],
										axis=1),
									color = 'black',
									linestyles = 'solid',
									linewidth = 2)
		ax.add_collection(bound_segments)
		ax.plot(vertices[edges[boundary_edges]-1,0],
				vertices[edges[boundary_edges]-1,1],
					marker = '.', color = 'black',
					linestyle = '', markersize = 7)
		ax.plot(vertices[edges[boundary_edges]-1,0],
				vertices[edges[boundary_edges]-1,1],
					marker = '.', color = 'tab:orange',
					linestyle = '', markersize = 3)
	ax.set_box_aspect(1)
	furthest = np.amax(np.linalg.norm(vertices,axis=1))+5e-2
	ax.set_xlim(-furthest,furthest)
	ax.set_ylim(-furthest,furthest)
	ax.tick_params(axis='both', which='both',
				   bottom=False, top=False, labelbottom=False,
				   right=False, left=False, labelleft=False)
	return ax

################################################################################

def parse_datafile (datafile):
	vertices = np.zeros((0,2), dtype = float)
	edges = np.zeros((0,2), dtype = int)
	faces = np.empty(0, dtype = object)
	input_type = 'none'
	with open(datafile) as instream:
		for line in instream:
			data = None
			if line == '\n':
				continue
			elif line == '# vertices\n':
				input_type = 'vertices'
			elif line == '# edges\n':
				input_type = 'edges'
			elif line == '# faces\n':
				input_type = 'faces'
			else:
				data = line.split('\t')[1:]
			if data is None:
				continue
			elif input_type == 'vertices':
				vertices = np.append(vertices,
							[[float(data[0]),float(data[1])]],
								axis=0)
			elif input_type == 'edges':
				edges = np.append(edges,
							[[int(data[0]),int(data[1])]],
								axis=0)
			elif input_type == 'faces':
				if int(data[-1]) == -1:
					continue
				else:
					faces = np.append(faces, None)
					faces[-1] = np.array([int(edge) for edge in data[:-1]],
											dtype = int)
	edge_count = np.zeros(edges.shape[0], dtype = int)
	for face in faces:
		for edge in face:
			edge_count[np.abs(edge)-1] += 1
	boundary_edges = (edge_count == 1)
	return vertices, edges, faces, boundary_edges

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
						description = 'Generate initial condition')
	parser.add_argument('datafile',
						nargs = '?',
						default = '../bounded_output/state.txt',
						type = str,
						help = 'File to plot.')
	args = parser.parse_args()
	datafile = Path(args.datafile)
	if datafile.exists():
		vertices, edges, faces, boundary_edges = parse_datafile(datafile)
		from matplotlib import pyplot as plt
		fig = plt.figure()
		ax = plt.axes()
		ax = plot_system(ax, vertices, edges, faces, boundary_edges)
		plt.savefig(datafile.with_suffix('.svg'))

################################################################################
# EOF
