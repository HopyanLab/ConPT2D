#!/usr/bin/env /usr/bin/python3
import numpy as np

################################################################################
#===============================================================================
# tile.py
#===============================================================================
################################################################################

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

if __name__ == '__main__':
	pass

################################################################################
# EOF
