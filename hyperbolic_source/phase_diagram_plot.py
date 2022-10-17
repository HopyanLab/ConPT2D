#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

n_params = np.array([24, 32, 48, 64, 128])
curvatures = 4*np.pi/n_params
curvatures = np.append(curvatures,0)
high_errors = np.array([0.016, 0.008, 0.004, 0.002, 0.002])
low_errors = np.ones_like(high_errors)*0.002
sigmas = np.array([4.5,2.5,1.5,1,1,0.01]) #([8,4,2,1,1,0.01])
markers = np.array(['o', 'v', 's', 'P', 'X','.'])
colors = np.array(['purple', 'blue', 'green',
					'darkorange', 'red', 'black'])
area_params = np.pi*np.array([2,1,1/2,1/4,1/8])
area_labels = np.array(['2$\pi$','$\pi$','$\pi/2$','$\pi/4$','$\pi/8$'])
neg_curvatures = -area_params/96#(4*np.pi)
neg_markers = np.array(['h', 'H', 'p', 'd', '^'])
neg_colors = np.array(['tab:purple', 'tab:blue', 'tab:green',
					'tab:orange', 'tab:red'])
neg_sigmas = np.array([1,1,1,1,1])
p0_flat = 3.8119822257016356
p0_fits = np.array([
					3.7474771420583752,
					3.761737806413651,
					3.7791878187623595,
					3.787693819705085,
					3.7993976597594026
					])

p0_neg_fits = p0_flat + np.array([
					0.051411663141728604,
					0.029691876418361636,
					0.015751701525123982,
					0.008985987951775805,
					0.004768027473821682
					])
#p0_neg_fits = np.array([3.95895526,
#						3.92341919,
#						3.89729427,
#						3.89160529,
#						3.88855751]) - 3.88280595 + p0_flat
p0_fits_corrected = p0_fits + (high_errors - low_errors)/2
p0_fits = np.append(p0_fits, p0_flat)
p0_fits_corrected = np.append(p0_fits_corrected, p0_flat)

fig, ax = plt.subplots(1)
x_values = np.linspace(3.,4.,11)
ax.plot(x_values, np.zeros_like(x_values),
		color = 'darkgrey', linestyle = '-', marker = '',
		linewidth = 2.0)
ax.grid(True)

def SRP (R, n, A = 1): #spherical_regular_perimeter
	# (radius, number_sides, area)
	SL = 2*R*np.arccos( np.cos(np.pi/n) / np.cos(A/2/R**2/n - np.pi/n) )
	# side_length
	return n*SL

R_params = np.sqrt(n_params/4/np.pi)
SRP_5 = SRP(R_params, 5)
SRP_4 = SRP(R_params, 4)
#estimates = ((n_params-12)*SRP_5 + 12*SRP_4)/n_params
estimates = SRP_5
estimates = np.append(estimates, 3.812)
neg_SRP_5 = 2*p0_flat + SRP(-np.sqrt(1/-neg_curvatures), 5)
estimates = np.append(estimates, neg_SRP_5)
#print(estimates)
#ax.plot(estimates, np.append(curvatures, neg_curvatures),
#		color = 'purple', linestyle = 'dotted', marker = '')

################################################################################
#def F (x,a,b,c):
#	return a*x**2+b*x+c
#test_point = [1,1,-3.813]

#def F (x,a,b):
#	return a*x+b
#test_point = [-1,3.813]

#params, covar = curve_fit(F, curvatures, p0_fits_corrected,
#							p0 = test_point, sigma = sigmas)

def F(x,a_p,a_n,b,c):
	return (a_p*x+b)*((np.tanh((-x)/c) + 1)/2) + \
		   (a_n*x+b)*((np.tanh((x)/c) + 1)/2)
test_point = [-1,-3,3.812,1]

params, covar = curve_fit(F,
							np.append(curvatures, neg_curvatures),
							np.append(p0_fits_corrected, p0_neg_fits),
							p0 = test_point,
							sigma = np.append(sigmas, neg_sigmas)
							)

print(params)
x = np.linspace(-0.8,0.8,1000)
ax.fill(np.append(F(x,*params),F(x[-1],*params)), np.append(x,x[0]),
	#	facecolor = 'lightblue',
		facecolor = 'lightgray',
		alpha = 0.6)
#ax.fill(np.append(F(x,*params),(F(x[0],*params))), np.append(x,x[-1]),
#		facecolor = 'pink',
#		alpha = 0.3)
ax.plot(F(x,*params), x, color = 'gray', linestyle = 'dashed', marker = '')
################################################################################


#x = np.linspace(3.7,3.813,100)
black_handles = np.empty(len(n_params)+len(area_params), dtype=object)
color_handles = np.empty(len(n_params)+len(area_params), dtype=object)
curve_handles = np.empty(len(n_params)+len(area_params), dtype=object)
for index, n_param in enumerate(n_params):
	black_handles[index],_,curve_handles[index] = \
				ax.errorbar(p0_fits[index],
							curvatures[index],
							xerr = np.array([[low_errors[index]],
											 [high_errors[index]]]),
							color = 'black', #colors[index],
							linestyle = '',
							marker = markers[index],
							markersize = 6.,
						#	alpha = 0.5,
							zorder = 4)
	color_handles[index], = ax.plot(p0_fits[index],
							curvatures[index],
							color = colors[index],
							linestyle = '',
							marker = markers[index],
							markersize = 3.,
							label = 'N = {0:d}'.format(n_param),
							zorder = 5)
for index, area_param in enumerate(area_params[::-1]):
	black_handles[len(n_params)+index],_,curve_handles[len(n_params)+index] = \
				ax.errorbar(p0_neg_fits[index],
							neg_curvatures[index],
							xerr = np.array([[low_errors[index]],
											 [low_errors[index]]])*2,
							color = 'black', #neg_colors[index],
							linestyle = '',
							marker = neg_markers[index],
							markersize = 6.,
						#	alpha = 0.5,
							zorder = 4)
	color_handles[len(n_params)+index], = ax.plot(p0_neg_fits[index],
							neg_curvatures[index],
							color = neg_colors[index],
							linestyle = '',
							marker = neg_markers[index],
							markersize = 3.,
							label = r'A = {0:s}'.format(area_labels[index]),
							zorder = 5)
flat_black_handle, = ax.plot(p0_flat, 0,
		color = 'black', linestyle = '',
		marker = '.', markersize = 8., zorder = 4)
flat_grey_handle, = ax.plot(p0_flat, 0,
		color = 'grey', linestyle = '',
		marker = '.', markersize = 3., zorder = 5,
		label = r'N $\rightarrow \infty$ (flat)')
ax.text(3.77,0.00,'solid-like', color='black',
		bbox=dict(facecolor='white', edgecolor='black',
					boxstyle='round,pad=0.4'))
ax.text(3.81,0.36,'fluid-like', color='black',
		bbox=dict(facecolor='white', edgecolor='black',
					boxstyle='round,pad=0.4'))
ax.set_ylim((-0.22,0.62))
ax.set_xlim((3.738,3.882))
ax.set_xlabel('$p_0$')
ax.set_ylabel('curvature (4$\pi$/N or A/96)')
handles,labels = ax.get_legend_handles_labels()
combined_handles = handles
combined_handles[:-1] = zip(curve_handles,
							black_handles,
							color_handles)
combined_handles[-1] = (flat_black_handle, flat_grey_handle)
order = np.array([0,1,2,3,4,10,9,8,7,6,5])
combined_handles = np.array(combined_handles,dtype=object)[order]
labels = np.array(labels,dtype=object)[order]
ax.legend(combined_handles, labels,
		  loc = 'best', fancybox = True, framealpha = 1.)
plt.savefig('../plots/phase_diagram.svg')
plt.rc('pgf', texsystem='pdflatex')
plt.savefig('../plots/phase_diagram.pgf')
#plt.show()
