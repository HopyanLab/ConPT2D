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
p0_flat = 3.812
#p0_fits = np.array([3.7451,3.7711,3.7898,3.8016,3.8098])
p0_fits = np.array([
					3.7474771420583752,
					3.761737806413651,
					3.7791878187623595,
					3.787693819705085,
					3.7993976597594026
					])
p0_fits_corrected = p0_fits + (high_errors - low_errors)/2
p0_fits = np.append(p0_fits, p0_flat)
p0_fits_corrected = np.append(p0_fits_corrected, p0_flat)

def SRP (R, n, A = 1): #spherical_regular_perimeter
	# (radius, number_sides, area)
	SL = 2*R*np.arccos( np.cos(np.pi/n) / np.cos(A/2/R**2/n - np.pi/n) )
	# side_length
	return n*SL

R_params = np.sqrt(n_params/4/np.pi)
SRP_5 = SRP(R_params, 5)
SRP_4 = SRP(R_params, 4)
estimates = ((n_params-12)*SRP_5 + 12*SRP_4)/n_params
estimates = np.append(estimates, 3.812)

#def F (x,a,b,c):
##	return a*np.sqrt(b-x)+c
#	return a*x**2+b*x+c
#params, covar = curve_fit(F, curvatures, p0_fits_corrected,
##							p0=[1,3.813,1], sigma=sigmas)
#							p0=[1,1,-3.813], sigma=sigmas)
def F (x,a,b):
	return a*x+b
params, covar = curve_fit(F, curvatures, p0_fits_corrected,
							p0=[-1,-3.812], sigma=sigmas)
fig, ax = plt.subplots(1)
x_values = np.linspace(3.,4.,11)
ax.plot(x_values, np.zeros_like(x_values),
		color = 'darkgrey', linestyle = '-', marker = '')
#x = np.linspace(3.7,3.813,100)
x = np.linspace(0,0.8,100)
ax.fill(np.append(F(x,*params),F(x[-1],*params)), np.append(x,x[0]),
		facecolor = 'lightgray', alpha = 0.5)
#ax.fill(np.append(F(x,*params),(F(0.7,*params),3.85,3.85)),
#		np.append(x,(0.7,0.7,0)),
#		facecolor = 'pink', alpha = 0.5)
ax.plot(F(x,*params), x, color = 'gray', linestyle = '--', marker = '')
#ax.plot(estimates, curvatures)
ax.grid(True)
for index, n_param in enumerate(n_params):
	ax.errorbar(p0_fits[index],
				curvatures[index],
				xerr = np.array([[low_errors[index]],[high_errors[index]]]),
				color = 'black', #colors[index],
				linestyle = '',
				marker = markers[index],
			#	alpha = 0.5,
				zorder = 4)
	ax.plot(p0_fits[index],
			curvatures[index],
			color = colors[index],
			linestyle = '',
			marker = markers[index],
			label = 'N = {0:d}'.format(n_param),
			zorder = 5)
ax.plot(p0_flat, 0,
		color = 'black', linestyle = '',
		marker = '.', label = r'N $\rightarrow \infty$ (flat)',
		zorder = 5)
ax.text(3.755,0.144,'solid-like', color='black',
		bbox=dict(facecolor='white', edgecolor='black',
					boxstyle='round,pad=0.4'))
ax.text(3.775,0.444,'fluid-like', color='black',
		bbox=dict(facecolor='white', edgecolor='black',
					boxstyle='round,pad=0.4'))
ax.set_ylim((-0.02,0.62))
ax.set_xlim((3.738,3.822))
ax.set_xlabel('$p_0$')
ax.set_ylabel('curvature (4$\pi$/N)')
handles,labels = ax.get_legend_handles_labels()
handles = np.array(handles,dtype=object)
labels = np.array(labels,dtype=object)
ax.legend(handles,labels,
		  loc = 'best', fancybox = True, framealpha = 1.)
plt.savefig('../plots/phase_diagram.svg')
plt.rc('pgf', texsystem='pdflatex')
plt.savefig('../plots/phase_diagram.pgf')
#plt.show()
