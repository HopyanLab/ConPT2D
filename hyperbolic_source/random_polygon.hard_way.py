#!/usr/bin/env /usr/bin/python3
import numpy as np

################################################################################
#===============================================================================
# random_polygon.py
#===============================================================================
################################################################################

def random_polygon (N = 64, polygon_sides = 8):
	# Generate random points in our polygon. (Poisson Process)
	np.random.seed()
	points = np.zeros((N,2))
	mirror_points = np.zeros((N,2))
	x1, y1 = 1, 0
	x2, y2 = np.cos(2*np.pi/polygon_sides), np.sin(2*np.pi/polygon_sides)
	slope = (y2-y1)/(x2-x1)
	inter = -slope*x1
	slope_perp = -1/slope
	for index in range(N):
		# pick in a triangle defined by origin, (x1,y1), and (x2,y2)
		# polygon is just N of those rotated around origin
		# so we can rotate the point into a random triangle
		angle_add = np.random.randint(0, polygon_sides, size = None) * \
						2*np.pi/polygon_sides
		repick = True
		while repick:
			angle = np.random.uniform(0, 2*np.pi/polygon_sides, size = None)
			length = np.sqrt(np.random.uniform(0, 1, size = None))
			x0, y0 = length*np.cos(angle), length*np.sin(angle)
			inter_perp = y0 - slope_perp*x0
			x_cross = (inter_perp - inter)/(slope - slope_perp)
			if x_cross > x0: # point is inside triangle
			# second condition only needed for N = 3
				repick = False
				x3 = 2*x_cross - x0
				y3 = slope_perp*x3 + inter_perp
				# rotate to proper triangle
				x_point = x0*np.cos(angle_add) - y0*np.sin(angle_add)
				y_point = x0*np.sin(angle_add) + y0*np.cos(angle_add)
				x_mirror = x3*np.cos(angle_add) - y3*np.sin(angle_add)
				y_mirror = x3*np.sin(angle_add) + y3*np.cos(angle_add)
				points[index] = np.array([x_point,y_point])
				mirror_points[index] = np.array([x_mirror,y_mirror])
	boundary_points = np.zeros((polygon_sides,2))
	for index in range(polygon_sides):
		boundary_points[index,0] = np.cos(2*np.pi*index/polygon_sides)
		boundary_points[index,1] = np.sin(2*np.pi*index/polygon_sides)
	return points, mirror_points, boundary_points

################################################################################

if __name__ == '__main__':
	points, mirror_points, boundary_points = random_polygon(600, 3)
	boundary_points = np.append(boundary_points,
								boundary_points[0][np.newaxis,:],
									axis=0)
	from matplotlib import pyplot as plt
	plt.plot(points[:,0], points[:,1],
				marker = '.', color = 'tab:blue',
				linestyle = '')
	plt.plot(mirror_points[:,0], mirror_points[:,1],
				marker = '.', color = 'tab:red',
				linestyle = '')
	plt.plot(boundary_points[:,0], boundary_points[:,1],
				marker = '.', color = 'black',
				linestyle = '-')
	plt.show()
