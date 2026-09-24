#!/usr/bin/env python3
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import geometric_slerp
import argparse
from pathlib import Path
from parse_data import parse_data

################################################################################
#===============================================================================
# plot_system.py
#===============================================================================
################################################################################

def subdivide(verts, faces):
	for index, face in enumerate(faces):
		# Create three new verts at the midpoints of each edge:
		new_verts = (verts[face] + verts[np.roll(face,-1)])/2
		new_verts = new_verts / np.linalg.norm(new_verts, axis=1)[:,np.newaxis]
		verts = np.append(verts, new_verts, axis=0)
		# Split the current triangle into four smaller triangles:
		i = len(verts) - 3
		j, k = i+1, i+2
		faces = np.append(faces, [[i, j, k]], axis=0)
		faces = np.append(faces, [[face[0], i, k]], axis=0)
		faces = np.append(faces, [[i, face[1], j]], axis=0)
		faces[index] = [k, j, face[2]]
	return verts, faces

################################################################################

def plot_points(vertices, ax):
#	alphas = (vertices[:,2] >= 0) + \
#			 (vertices[:,2] < 0) * (0.7 + 0.5*vertices[:,2])
#	ax.plot(vertices[:,0],
#			vertices[:,1],
#			vertices[:,2],
#			marker = 'o',
#			markersize = 7.,
#			linestyle = '',
#			color = 'black',
#			zorder = 10)
#	ax.plot(vertices[:,0],
#			vertices[:,1],
#			vertices[:,2],
#			marker = 'o',
#			markersize = 3.,
#			linestyle = '',
#			color = 'tab:red',
#			zorder = 11)
	ax.scatter( vertices[:,0],
				vertices[:,1],
				vertices[:,2],
				marker = 'o',
				s = 36.,
				linewidths = 2.,
				edgecolor = 'black',
				facecolor = 'tab:red',
				zorder = 10)

################################################################################

def plot_edges(vertices, edges, ax,  subdiv = 2):
	t_vals = np.linspace(0, 1, 2**subdiv+1)
	segments = np.empty((len(edges), 2**subdiv+1, 3), dtype = float)
	alphas = np.empty((len(edges)), dtype = float)
	for index, edge in enumerate(edges):
		segments[index] = geometric_slerp(vertices[edge[0]],
										  vertices[edge[1]],
										  t = t_vals)
		alphas[index] = (np.mean(segments[index,:,2]) >= 0) + \
						(np.mean(segments[index,:,2]) < 0) * \
									(0.6+0.5*np.mean(segments[index,:,2]))
	lc = Line3DCollection(segments,
						  colors='black',
						  alpha = alphas,
						  linewidth = 2.,
						  zorder = 8)
	ax.add_collection(lc)

################################################################################

def plot_faces(vertices, faces, ax, subdiv = 2):
	verts = vertices
	tri_faces = np.empty((0,3), dtype = int)
	for face in faces:
		centroid = np.mean(verts[face], axis=0)
		tri_faces = np.append(tri_faces,
								[[face[-1], face[0], len(verts)]],
								axis=0)
		for index in range(len(face)-1):
			tri_faces = np.append(tri_faces,
								[[face[index], face[index+1], len(verts)]],
								axis=0)
		verts = np.append(verts, centroid[np.newaxis,:], axis=0)
	for level in range(subdiv):
		verts, tri_faces = subdivide(verts, tri_faces)
	polygons = verts[tri_faces]
	pc = Poly3DCollection(polygons,
						  linewidth = 0,
						  color = np.array([0x00, 0xbe, 0xe0],
											dtype = float)/0xff,
#						  color = np.array([0x00, 0x66, 0xb5],
#											dtype = float)/0xff,
						  zorder = 6)
	ax.add_collection(pc)

################################################################################

def save_plot(outfile):
	plt.savefig(outfile, format="svg")
	plt.close()

################################################################################

def plot_system (
					savefile = Path('../sphere_output/system_state.txt'),
					plotfile = Path('../plots/sphere_system_plot.svg'),
					subdiv = 2
				):
	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(111, projection='3d', computed_zorder=False)
	vertices, edges, faces = parse_data(savefile)
	vertices = vertices / np.linalg.norm(vertices, axis=1)[:,np.newaxis]
	plot_points(vertices, ax)
	plot_edges(vertices, edges, ax, subdiv = subdiv)
#	plot_faces(vertices, faces, ax, subdiv = subdiv)
	ax.set_box_aspect((1,1,1))
	ax.set_axis_off()
	ax.elev = 90
	ax.azim = 0
	ax.dist = 6
#	ax.view_init(elev = 90, azim = 0)
	plt.tight_layout()
#	plt.show()
	if plotfile.suffix == '.svg':
		plt.savefig(plotfile, format='svg')
	elif plotfile.suffix == '.png':
		plt.savefig(plotfile, format='png')
	elif plotfile.suffix == '.pgf':
		plt.rc('pgf', texsystem='pdflatex')
		plt.savefig(plotfile.with_suffix('.pgf'))
	plt.close()

################################################################################

if __name__ == '__main__':
	# Filenames can be passed as arguements so
	#  we use argparse module to parse them.
	parser = argparse.ArgumentParser(
							description = 'Plot spherical Voronoi diagram')
	# outfile = sys.argv[1]
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = ['../plots/sphere_system_plot.svg'],
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('-i', '--infile',
						nargs = 1,
						default = ['../sphere_output/system_40.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	args = parser.parse_args()
	subdiv = 2
	if Path(args.infile[0]).exists():
		plot_system(
					savefile = Path(args.infile[0]),
					plotfile = Path(args.outfile[0]),
					subdiv = subdiv
					)
	for statefile in Path(args.infile[0]).parent.glob('*.txt'):
		if statefile.name.split('_')[0] == 'system':
			plot_system(statefile, statefile.with_suffix('.png'),
						subdiv = subdiv)

################################################################################
# EOF
