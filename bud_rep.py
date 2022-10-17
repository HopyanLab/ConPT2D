#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from sympy.abc import x, y
from sympy import ordered, Matrix, hessian, exp, sqrt, lambdify

################################################################################

#def gaussian (a, mu=0, sigma=1):
#	return np.exp(-(a-mu)**2/sigma**2/2)/(sigma*np.sqrt(2*np.pi))

#def f(a, b, mu_a=0, sigma_a=1,
#			mu_b=0, sigma_b=1):
#	return gaussian(a,mu_a,sigma_a) * gaussian(b,mu_b,sigma_b)

#
#def K(x, y, mu_x=0, sigma_x=1,
#			mu_y=0, sigma_y=1):
#	return

#if __name__ == "__main__":
#	#
#	x,y = np.meshgrid(np.linspace(-4, 4, 64), np.linspace(-4, 4, 64))
#	z = f(x,y,0,1.5,0,1)*10
#	#
#	fig = plt.figure(figsize =(14, 9))
#	ax = plt.axes(projection ='3d')
#	# Creating plot
#	surf = ax.plot_surface(x, y, z,
#					rstride=1, cstride=1,
#					shade=False,
#					cmap=mpl.colors.ListedColormap(['white']), #'jet',
#					linewidth=1)
#	ax.set_aspect('auto')
#	plt.draw()
#	m = plt.cm.ScalarMappable(surf.norm)
#	surf.set_edgecolors(m.to_rgba(surf.get_array()))
#	# surf.set_edgecolors(surf.to_rgba(surf._A))
#	# surf.set_facecolors("white")
#	# show plot
#	plt.show()

def gaussian (a, mu=0, sigma=1):
	return np.exp(-(a-mu)**2/sigma**2/2)/(sigma*np.sqrt(2*np.pi))

def sym_gaussian (a, mu=0, sigma=1):
	return exp(-(a-mu)**2/sigma**2/2)/(sigma*sqrt(2*np.pi))

if __name__ == "__main__":
	gradient = lambda f, v: Matrix([f]).jacobian(v)
	z = sym_gaussian(x,0,1.5)*sym_gaussian(y,0,1.)
	v = list(ordered(z.free_symbols));
	H_det = hessian(z,v).det()
	K = H_det/(1+gradient(z,v).dot(gradient(z,v)))**2
	# det(hessian(f))/(1+grad(f)^2)^2
	curvature = lambdify(v,K)
	#
	xx,yy = np.meshgrid(np.linspace(-4, 4, 64), np.linspace(-4, 4, 64))
	zz = gaussian(xx,0,1.5)*gaussian(yy,0,1.)*30
	kk = curvature(xx,yy)
	fig = plt.figure(figsize =(14, 9))
	ax = plt.axes(projection ='3d')
	surf = ax.plot_surface(xx, yy, zz,
					rstride=1, cstride=1,
					shade=False,
					cmap=mpl.colors.ListedColormap(['lightgrey']), #'jet',
					linewidth=1)
	ax.set_aspect('auto')
	ax.set_box_aspect((np.amax(xx) - np.amin(xx),
					   np.amax(yy) - np.amin(yy),
					   np.amax(zz) - np.amin(zz)))
	ax.xaxis.set_ticklabels([])
	ax.yaxis.set_ticklabels([])
	ax.zaxis.set_ticklabels([])
	plt.draw()
	#m = plt.cm.ScalarMappable(surf.norm)
	#surf.set_edgecolors(m.to_rgba(surf.get_array()))
	#colormap = 'jet'
	#colormap = 'hot'
	#colormap = 'viridis'
	#colormap = 'RdBu'
	k_min = np.amin(kk)
	k_max = np.amax(kk)
	full_range = k_max - k_min
	num_pos = int(np.floor(256*k_max/full_range))
	num_neg = int(np.ceil(-256*k_min/full_range))
	colors_pos = plt.cm.Blues(np.linspace(0.,1.,num_pos))
	colors_neg = plt.cm.Reds(np.linspace(0.,0.5,num_neg))
	colors_full = np.vstack((colors_neg[::-1], colors_pos))
	colormap = mpl.colors.LinearSegmentedColormap.from_list(
						'custom_colormap', colors_full)
	cm = plt.get_cmap(colormap)#.reversed()
	cnorm = mpl.colors.Normalize(vmin=np.amin(kk), vmax=np.amax(kk))
	scalar_map = plt.cm.ScalarMappable(norm=cnorm, cmap=cm)
	scalar_map.set_array(kk)
	colors = scalar_map.to_rgba(kk[:-1,:-1].reshape(surf.get_array().shape))
	surf.set_edgecolors(colors)
#	cbar = plt.colorbar(scalar_map, ticks=[0.], shrink = 0.5)
#	cbar.ax.set_yticklabels(['0'])
#	cbar.set_label('Gaussian curvature', rotation=270)
	plt.axis('off')
	plt.tight_layout()
	plt.savefig('./bud_rep.svg', bbox_inches='tight')
	plt.show()

################################################################################
# EOF


