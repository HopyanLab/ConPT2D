#!/usr/bin/env /usr/bin/python3
import numpy as np 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

################################################################################
#===============================================================================
# plot_system.py
#===============================================================================
################################################################################

# Difference with respect to periodic boundaries
def periodic_diff (point1, point2, L = np.array([1., 1.])):
	return ((point1 - point2 + L/2.) % L) - L/2.

################################################################################

def periodic_plot (point1, point2, color, box_size = 1.):
	L = np.array([box_size, box_size])
	point2_wrapped = point1 + periodic_diff(point2, point1, L)
	point1_wrapped = point2 + periodic_diff(point1, point2, L)
	plt.plot((point1_wrapped[0], point2[0]), (point1_wrapped[1], point2[1]),
				c = color, zorder = 0)
	plt.plot((point1[0], point2_wrapped[0]), (point1[1], point2_wrapped[1]),
				c = color, zorder = 0)
	return

################################################################################

def plot_points(vertices):
	for x, y in vertices:
		plt.scatter(x, y,
					s = 16.,
					edgecolor = 'xkcd:fuchsia',
					facecolor = 'xkcd:light pink',
					zorder = 10)
	return

################################################################################

def plot_edges(vertices, edges):
	for i1, i2 in edges:
		if i1 != -1 and i2 != -1:
			x1, y1 = vertices[i1-1]
			x2, y2 = vertices[i2-1]
			periodic_plot(np.array([x1, y1]), np.array([x2, y2]),
							color = 'xkcd:cerulean')
	return

################################################################################
# TODO: Fix this.
#def plot_faces(vertices, faces):
#	for face in faces:
#		for i,index in enumerate(face):
#			if index == face[-1]:
#				if index != -1 and face[0] != -1:
#					x1,y1 = vertices[index]
#					x2,y2 = vertices[face[0]]
#					periodic_plot(x1, y1, x2, y2,
#									color = 'xkcd:magenta')
#			else:
#				if index != -1 and face[i+1] != -1:
#					x1,y1 = vertices[index]
#					x2,y2 = vertices[face[i+1]]
#					periodic_plot(x1, y1, x2, y2,
#									color = 'xkcd:magenta')
#	return

################################################################################

def save_plot(outfile, box_size = 1.):
	# remove tick marks
	frame = plt.gca()
	frame.spines['bottom'].set_color('xkcd:goldenrod')
	frame.spines['bottom'].set_linewidth(1.6)
	frame.spines['top'].set_color('xkcd:goldenrod')
	frame.spines['top'].set_linewidth(1.6)
	frame.spines['right'].set_color('xkcd:seafoam')
	frame.spines['right'].set_linewidth(1.6)
	frame.spines['left'].set_color('xkcd:seafoam')
	frame.spines['left'].set_linewidth(1.6)
	frame.axes.get_xaxis().set_ticks([0,1])
	frame.set_xticklabels(['',''])
	frame.xaxis.set_ticks_position('both')
	frame.tick_params(axis='x',
						direction='in',
						length = 4.,
						width = 3.,
						colors='xkcd:red')
	frame.axes.get_yaxis().set_ticks([0,1])
	frame.set_yticklabels(['',''])
	frame.yaxis.set_ticks_position('both')
	frame.tick_params(axis='y',
						direction='in',
						length = 4.,
						width = 3.,
						colors='xkcd:red')
	frame.set_facecolor('none')
	plt.axis([0, box_size, 0, box_size])
	plt.scatter([0., 0., box_size, box_size],
				[0., box_size, 0., box_size],
				s = 64.,
				edgecolor = 'xkcd:red',
				facecolor = 'none')
	# save and close plot
	plt.savefig(outfile, format="svg")
	plt.close()
	return

################################################################################
# TODO: Fix this.
#def read_face(file):
#	indices = []
#	with open(file) as infile:
#		for line in infile:
#			cell_indices = []
#			linesplit = line.strip().split('\t')
#			for i in linesplit:
#				cell_indices.append(int(i))
#			indices.append(cell_indices)
#	return indices

################################################################################

def plot_system (
					vertfile = Path('../flat_output/vertices.txt'),
					edgefile = Path('../flat_output/edges.txt'),
					facefile = Path('../flat_output/faces.txt'),
					savefile = Path('../plots/system_plot.svg'),
					box_size = 1
				):
	vertices = np.loadtxt(vertfile, delimiter='\t')
	vertices[:,0] += (vertices[:,0] < 0) * box_size
	vertices[:,0] -= (vertices[:,0] > box_size) * box_size
	vertices[:,1] += (vertices[:,1] < 0) * box_size
	vertices[:,1] -= (vertices[:,1] > box_size) * box_size
	edges = np.loadtxt(edgefile, delimiter='\t', dtype=str)
	edges = edges[:,0:2].astype(int)
	#faces = read_face(facefile)
	#
	# plot_faces(vertices, faces)
	plot_edges(vertices, edges)
	plot_points(vertices)
	#
	save_plot(savefile, box_size = box_size)

################################################################################

if __name__ == '__main__':
	# Filenames can be passed as arguements so
	#  we use argparse module to parse them.
	parser = argparse.ArgumentParser(
							description = 'Plot periodic Voronoi diagram')
	# outfile = sys.argv[1]
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = ['../plots/system_plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('-vf', '--vertfile',
						nargs = 1,
						default = ['../flat_output/vertices.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	parser.add_argument('-ef', '--edgefile',
						nargs = 1,
						default = ['../flat_output/edges.txt'],
						type = str,
						required = False,
						help = 'file to put edge data in')
	parser.add_argument('-pf', '--facefile',
						nargs = 1,
						default = ['../flat_output/faces.txt'],
						type = str,
						required = False,
						help = 'file to put face data in')
	args = parser.parse_args()
	plot_system(vertfile = Path(args.vertfile[0]),
				edgefile = Path(args.edgefile[0]),
				facefile = Path(args.facefile[0]),
				savefile = Path(args.outfile[0]),
				box_size = 1.
				)

################################################################################
# EOF
