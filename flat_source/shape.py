#!/usr/bin/env /usr/bin/python3
import numpy as np 
import argparse
from geometry import *

################################################################################
#===============================================================================
# shape.py
#===============================================================================
################################################################################

def analyse_shape (vertexfile = './output/vertices.txt',
					facefile = './output/faces.txt'):
	# Load verticies
	vertices = np.loadtxt(vertexfile)
	# We need faces in terms of actual verticies rather than references
	faces = []
	# So we associate those references to the actual verticies
	with open(facefile) as infile:
		for line in infile:
			linesplit = line.strip().split()
			n = len(linesplit)
			face = np.zeros((n,2))
			for i,x in enumerate(linesplit):
				face[i,:] = vertices[int(x),:]
			faces.append(face)
	# Then work out the shape indices for the faces
	shapes = np.zeros((0, 3))
	for face in faces:
		n = len(face)
		a = abs(get_area(n,face))
		p = get_perimeter(n,face)
		s = p**2 / np.sqrt(a)
		print('perimeter = ', p,
			'\narea = ', a,
			# '\nshape [P^2/(4*pi*A)] = ', p**2 / (4.*np.pi*a),
			'\nshape [P/sqrt(A)] = ', s,
			'\n')
		np.append(shapes,np.array([[p,a,s]]),axis=0)
	return shapes


if __name__ == '__main__':
	# Check arguements
	parser = argparse.ArgumentParser(description='Shape indices for faces')
	# vertexfile = sys.argv[1]
	parser.add_argument('-v', '--vertexfile',
						metavar = 'vf',
						nargs = 1,
						default = './output/vertices.txt',
						type = str,
						required = False,
						help = 'File with vertex data')
	# facefile = sys.argv[2]
	parser.add_argument('-p', '--facefile',
						metavar = 'pf',
						nargs = 1,
						default = './output/faces.txt',
						type = str,
						required = False,
						help = 'File with face data')
	args = parser.parse_args()
	analyse_shape(args.vertexfile, args.facefile)

################################################################################
# EOF
