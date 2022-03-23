#!/usr/bin/env /usr/bin/python3
import numpy as np

################################################################################
#===============================================================================
# geometry.py
#===============================================================================
################################################################################

# Area of a polygon
def get_area(Nv,vertices):
	#
	cross = 0.
	for i in range(0,Nv):
		x1,y1 = vertices[i,:]
		if i != Nv - 1:
			x2,y2 = vertices[i+1,:]
		if i == Nv - 1:
			x2,y2 = vertices[0,:]
		cross += ((x1 * y2) - (x2 * y1))
	return 0.5 * cross

################################################################################

# Perimeter of a polygon
def get_perimeter(Nv,vertices):
	#
	lengths = np.zeros(Nv)
	for i in range(0,Nv):
		x1,y1 = vertices[i,:]
		if i !=  Nv - 1:
			x2,y2 = vertices[i+1,:]
		if i == Nv - 1:
			x2,y2 = vertices[0,:]
		d = np.sqrt((x2-x1)**2 + (y2-y1)**2)
		lengths[i] = d
	return np.sum(lengths)

################################################################################

if __name__ == '__main__':
	pass

################################################################################
# EOF
