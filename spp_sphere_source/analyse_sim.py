#!/usr/bin/env python3
import numpy as np
import os
import subprocess
import glob
from pathlib import Path
import argparse
from parse_data import parse_data
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import matplotlib as mpl
from matplotlib import pyplot as plt

from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D

mpl.rcParams["font.size"] = 14

################################################################################
#===============================================================================
# analyse_spherical_sim.py
#===============================================================================
################################################################################

def rotve(v,erot,angle):
	rotmeasure=np.linalg.norm(erot)
	erot=erot/rotmeasure;
	norme=np.dot(v,erot)
	vplane=v-norme*erot
	plnorm=np.linalg.norm(vplane)
	ep=vplane/plnorm
	eo=np.cross(erot,ep)
	vrot=(np.cos(angle)*ep+np.sin(angle)*eo)*plnorm+norme*erot
	return(vrot)

def rot2d (v,angle):
	rotMatrix = np.array([[np.cos(angle), -np.sin(angle)],
						  [np.sin(angle),  np.cos(angle)]])
	return(np.dot(rotMatrix,v))

################################################################################

def assemble_positions (out_dir, run_values, n_values, p0_values, time_values):
	run_values = np.delete(run_values, 5)
	data = np.zeros(0, dtype = object)
	for n_index, n_value in enumerate(n_values):
		temp_data = np.zeros((len(run_values),
							  len(p0_values),
							  n_value,
							  len(time_values),
							  3),
							dtype = float)
		for run_index, run_value in enumerate(run_values):
			run_dir = out_dir / f'run_{run_value:02d}' / f'n_{n_value:d}'
			for p0_index, p0_value in enumerate(p0_values):
				p0_dir = run_dir / f'p0_{p0_value:1.3f}'
				for time_index, time_value in enumerate(time_values):
					time_file = p0_dir / f'system_{time_value:d}.txt'
					vertices, edges, faces = parse_data(time_file)
					for face_index, face in enumerate(faces):
						position = np.mean(vertices[face], axis=0)
						position /= np.linalg.norm(position)
						temp_data[run_index, p0_index,
								  face_index, time_index] = position
		data = np.append(data,None)
		data[-1] = temp_data
	return data

################################################################################

def assemble_msd_flattened (data, run_values, n_values, p0_values, time_values):
	full_msd = np.zeros((len(n_values), len(p0_values),len(time_values)))
	for n_index, n_value in enumerate(n_values):
		temp_data = data[n_index]
		disp = np.zeros((len(run_values),
						 len(p0_values),
						 n_value,
						 len(time_values),
						 2), dtype = float)
		for run_index, run_value in enumerate(run_values):
			for p0_index, p0_value in enumerate(p0_values):
				for cell_index in range(n_value):
					for time_index, time_value in enumerate(time_values):
						if time_index == 1:
							position = temp_data[run_index,p0_index,
												 cell_index,time_index]
							previous = temp_data[run_index,p0_index,
												 cell_index,time_index-1]
							disp[run_index, p0_index,
								 cell_index, time_index,0] = \
								np.linalg.norm(position - previous)
						if time_index > 1:
							position = temp_data[run_index,p0_index,
												 cell_index,time_index]
							previous = temp_data[run_index,p0_index,
												 cell_index,time_index-1]
							two_ago = temp_data[run_index,p0_index,
												 cell_index,time_index-2]
							v1 = position - previous
							v0 = previous - two_ago
							length = np.linalg.norm(v1)
							v1 /= length
							v0 /= np.linalg.norm(v0)
							angle = np.arccos(np.clip(np.dot(v0, v1),
														-1.0, 1.0))
							prev_flat = disp[run_index, p0_index,
											 cell_index, time_index-1]
							two_ago_flat = disp[run_index, p0_index,
												cell_index, time_index-2]
							u0 = prev_flat - two_ago_flat
							u0 /= np.linalg.norm(u0)
							u1 = rot2d(u0, angle)*length
							disp[run_index, p0_index,
								 cell_index, time_index] = prev_flat + u1
		square_disp = disp[:,:,:,:,0]**2 + disp[:,:,:,:,1]**2
		msd = np.mean(square_disp, axis = (0,2)) # average over runs and faces
		full_msd[n_index] = msd
	return full_msd

################################################################################

def assemble_msd (data, run_values, n_values, p0_values, time_values):
	full_msd = np.zeros((len(n_values), len(p0_values),len(time_values)))
	for n_index, n_value in enumerate(n_values):
		r_value = np.sqrt(n_value/2/np.pi)
		temp_data = data[n_index]
		disp = np.zeros((len(run_values),
						 len(p0_values),
						 n_value,
						 len(time_values),
						 3), dtype = float)
		disp = temp_data - temp_data[: ,: ,: , 0, np.newaxis, :]
		scalar_disp = np.linalg.norm(disp , axis=-1)
		length_square = (scalar_disp*2)**2
		msd = np.mean(length_square, axis = (0,2)) # average over runs and faces
		full_msd[n_index] = msd
	return full_msd

################################################################################

def plot_msd (msd, n_values, p0_values, time_values):
	fig, ax = plt.subplots(1)#, figsize=(8, 6))
	cmap = plt.get_cmap('jet')
	p0_value = 3.78
	p0_index = int(np.argwhere(p0_values == p0_value)[0,0])
	lower = 10
	upper = int(np.argwhere(time_values > 1e4)[0,0])
	colors = np.array(['purple', 'blue', 'green',
						'darkorange', 'red'])
	for n_index, n_value in enumerate(n_values):
		y_values = msd[:, p0_index, :]
		ax.plot(time_values[lower:upper],
				y_values[n_index, lower:upper],
				color = colors[n_index],
		#		color = cmap((n_index)/(len(n_values)-1)),
				label = 'N = {0:d}'.format(n_value),
				zorder = 8-n_index)

	ax.set_title(f'$p_0 = {p0_value:1.3f}$')

	ax.text(1.6e3,4.8e-2,'Increasing',
				bbox = dict(facecolor='white', alpha=1,
						boxstyle='round', linestyle=''))
	ax.text(1.6e3,2.3e-2,'Curvature',
				bbox = dict(facecolor='white', alpha=1,
						boxstyle='round', linestyle=''))
	ax.annotate("", xy=(3e3, 1.3), xytext=(3e3, 8e-2),
					arrowprops=dict(arrowstyle='->',
									color = 'white',
									linewidth=4),
					zorder = 9)
	ax.annotate("", xy=(3e3, 1.2), xytext=(3e3, 8e-2),
					arrowprops=dict(arrowstyle='->',
									color = 'black',
									linewidth=2),
					zorder = 10)

	ax.set_xlabel('Time Steps')
	ax.set_ylabel('Mean Square Displacement')
	ax.grid(True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.legend(loc = 'best', fancybox = True, framealpha = 1.)
	fig.tight_layout()
	plt.savefig(f'../plots/sphere_msd_{p0_value:1.3f}.svg')
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig(f'../plots/sphere_msd_{p0_value:1.3f}.pgf')
#	plt.show()
	plt.close()
	return msd

################################################################################

def plot_Deff (msd, run_values, p0_values, time_values):
	fig, ax = plt.subplots(1)
	start = 24
	end = 36
	def F(x,a,b):
		return a*x+b
	D_eff = np.zeros_like(p0_values)
	for p0_index, p0_value in enumerate(p0_values):
		results = np.zeros(end+1-start)
		for index, point in enumerate(np.arange(start, end+1)):
			fit, cov = curve_fit(F, np.log(time_values[point:]),
									np.log(msd[p0_index,point:]))
			results[index] = fit[0] / 4 / ((0.1/0.02)**2/2/4.)
		D_eff[p0_index] = np.mean(results)
	ax.plot(p0_values, D_eff,
			marker = '',
			color = 'gray',
			linestyle = '-',
			zorder = 4)
	ax.plot(p0_values, D_eff,
			marker = '.',
			color = 'black',
			linestyle = '',
			markersize = 8,
			zorder = 5)
	ax.plot(p0_values, D_eff,
			marker = '.',
			color = 'gray',
			linestyle = '',
			markersize = 5,
			zorder = 6)
	fig.tight_layout()
	plt.show()

################################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
							description = '')
	parser.add_argument('datadir',
						nargs = '?',
						default = '../spp_sphere_output',
						type = str,
						help = 'directory of simulation output data')
	args = parser.parse_args()
	out_dir = Path(args.datadir)
	run_values = np.zeros(0, dtype = int)
	for run_dir in out_dir.glob('run_*'):
		run_values = np.append(run_values, int(run_dir.name.split('_')[-1]))
	run_values = np.sort(run_values)
	run_dir = out_dir / (f'run_{run_values[0]:02d}')
	n_values = np.zeros(0, dtype = int)
	for n_dir in run_dir.glob('n_*'):
		n_values = np.append(n_values, int(n_dir.name.split('_')[-1]))
	n_values = np.sort(n_values)
	n_dir = run_dir / (f'n_{n_values[0]:d}')
	p0_values = np.zeros(0, dtype = float)
	for p0_dir in n_dir.glob('p0_*'):
		p0_values = np.append(p0_values, float(p0_dir.name.split('_')[-1]))
	p0_values = np.sort(p0_values)
	p0_dir = n_dir / (f'p0_{p0_values[0]:1.3f}')
	time_values = np.zeros(0, dtype = int)
	for time_file in p0_dir.glob('system_*.txt'):
		time_value = int((time_file.name.split('.')[0]).split('_')[-1])
		time_values = np.append(time_values, time_value)
	time_values = np.sort(time_values)
	data_file = out_dir/'data.npy'
	if data_file.exists():
		data = np.load(data_file, allow_pickle = True)
	else:
		data = assemble_positions(out_dir, run_values, n_values,
										p0_values, time_values)
		np.save(data_file, data)
	msd_file = out_dir/'msd.npy'
	if msd_file.exists():
		msd = np.load(msd_file, allow_pickle = True)
	else:
		msd = assemble_msd(data, run_values, n_values, p0_values, time_values)
		np.save(msd_file, msd)
	plot_msd(msd, n_values, p0_values, time_values)
#	plot_Deff(msd, run_values, p0_values, time_values)

################################################################################
# EOF
