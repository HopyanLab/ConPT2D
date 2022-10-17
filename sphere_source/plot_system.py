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

################################################################################

def plot_points(vertices):
	for x, y in vertices:
		plt.scatter(x, y,
					s = 16.,
					edgecolor = 'xkcd:fuchsia',
					facecolor = 'xkcd:light pink',
					zorder = 10)

################################################################################

def plot_edges(vertices, edges):
	for i1, i2 in edges:
		if i1 != -1 and i2 != -1:
			x1, y1 = vertices[i1-1]
			x2, y2 = vertices[i2-1]
			periodic_plot(np.array([x1, y1]), np.array([x2, y2]),
							color = 'xkcd:cerulean')

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

################################################################################

def plot_system (
					savefile = Path('../flat_output/system_state.txt'),
					plotfile = Path('../plots/system_plot.svg'),
					box_size = 1.0
				):
	#TODO: fix this!
	vertices = np.loadtxt(vertfile, delimiter='\t')
	vertices[:,0] += (vertices[:,0] < 0) * box_size
	vertices[:,0] -= (vertices[:,0] > box_size) * box_size
	vertices[:,1] += (vertices[:,1] < 0) * box_size
	vertices[:,1] -= (vertices[:,1] > box_size) * box_size
	edges = np.loadtxt(edgefile, delimiter='\t', dtype=str)
	edges = edges[:,0:2].astype(int)
	plot_edges(vertices, edges)
	plot_points(vertices)
	#
	save_plot(plotfile, box_size = box_size)

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
	parser.add_argument('-i', '--infile',
						nargs = 1,
						default = ['../flat_output/system_state.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	args = parser.parse_args()
	plot_system(savefile = Path(args.infile[0]),
				plotfile = Path(args.outfile[0]),
				box_size = 1.
				)

################################################################################
# EOF
