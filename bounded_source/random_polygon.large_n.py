#!/usr/bin/env /usr/bin/python3
import numpy as np

################################################################################
#===============================================================================
# random_polygon.py
#===============================================================================
################################################################################

def random_polygon (N = 64, polygon_sides = 8):
	# Generate random points in our polygon.
	np.random.seed()
	points = np.zeros((N,2), dtype = float)
	mirror_points = np.zeros((3*N,2), dtype = float)
	x_cross = np.cos(-np.pi/polygon_sides)
	for index in range(N):
		# pick in a triangle defined by origin, (x1,y1), and (x2,y2)
		# polygon is just N of those rotated around origin
		# so we can rotate the point into a random triangle
		angle_add = (np.random.randint(0, polygon_sides, size = None) + 1/2) * \
						2*np.pi/polygon_sides
		repick = True
		while repick:
			angle = np.random.uniform(-np.pi/polygon_sides,
									   np.pi/polygon_sides, size = None)
			length = np.sqrt(np.random.uniform(0, 1, size = None))
			x_point, y_point = length*np.cos(angle), length*np.sin(angle)
			if x_cross > x_point: # point is inside triangle
				repick = False
				# rotate to proper triangle
				points[index] = rotate_point(x_point, y_point, angle_add)
				# Make points to build boundary
				mirror_points[3*index] = rotate_point(
								2*x_cross-x_point, y_point, angle_add)
				x_temp, y_temp = rotate_point(x_point, y_point,
												-2*np.pi/polygon_sides)
				mirror_points[3*index+1] = rotate_point(
										2*x_cross-x_temp, y_temp,
										angle_add + 2*np.pi/polygon_sides)
				x_temp, y_temp = rotate_point(x_point, y_point,
												2*np.pi/polygon_sides)
				mirror_points[3*index+2] = rotate_point(
										2*x_cross-x_temp, y_temp,
										angle_add - 2*np.pi/polygon_sides)
	return points, mirror_points

################################################################################

def rotate_point (x, y, theta = 0):
	return x*np.cos(theta) - y*np.sin(theta), \
		   x*np.sin(theta) + y*np.cos(theta)

################################################################################

# area of polygon is n*cos(2*pi/n)*sin(2*pi/n)
def length_correction (polygon_sides = 8):
	return np.sqrt(polygon_sides * np.sin(2*np.pi/polygon_sides) * \
								   np.cos(2*np.pi/polygon_sides))

################################################################################

if __name__ == '__main__':
	points, mirror_points = random_polygon()
	from matplotlib import pyplot as plt
	fig = plt.figure()
	ax = plt.axes()
	ax.plot(points[:,0], points[:,1],
				marker = '.', color = 'tab:blue',
				linestyle = '')
	ax.plot(mirror_points[:,0], mirror_points[:,1],
				marker = '.', color = 'tab:red',
				linestyle = '')
	ax.set_box_aspect(1)
	plt.show()

################################################################################
# EOF
