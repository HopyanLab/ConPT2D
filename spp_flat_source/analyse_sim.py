#!/usr/bin/env python3
import numpy as np
import os
import subprocess
import glob
from pathlib import Path
import argparse
from parse_data import parse_data
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib import pyplot as plt

mpl.rcParams["font.size"] = 14

################################################################################
#===============================================================================
# analyse_sim.py
#===============================================================================
################################################################################

def assemble_positions (out_dir, run_values, p0_values, time_values):
	data = []
	num_faces = 0
	size = 0
	for run_index, run_value in enumerate(run_values):
		run_dir = out_dir / f'run_{run_value:02d}'
		for p0_index, p0_value in enumerate(p0_values):
			p0_dir = run_dir / f'p0_{p0_value:1.3f}'
			for time_index, time_value in enumerate(time_values):
				time_file = p0_dir / f'time_{time_value:d}.txt'
				vertices, edges, faces = parse_data(time_file)
				if len(data) == 0:
					num_faces = len(faces)
					size = int(np.floor(np.sqrt(num_faces)))
					data = np.zeros((len(run_values),
									 len(p0_values),
									 num_faces,
									 len(time_values),
									  2),
								dtype = float)
				for face_index, face in enumerate(faces):
					unwrap_verts = vertices[face]
					for vert_index, vertex in enumerate(unwrap_verts):
						if vert_index > 0:
							last_vert = unwrap_verts[vert_index-1]
							if vertex[0] - last_vert[0] > size/2:
								unwrap_verts[vert_index, 0] -= size
							if last_vert[0] - vertex[0] > size/2:
								unwrap_verts[vert_index, 0] += size
							if vertex[1] - last_vert[1] > size/2:
								unwrap_verts[vert_index, 1] -= size
							if last_vert[1] - vertex[1] > size/2:
								unwrap_verts[vert_index, 1] += size
					data[run_index, p0_index, face_index, time_index] = \
						np.mean(unwrap_verts, axis=0)
				if time_index > 0:
					for face_index in range(num_faces):
						position = data[run_index, p0_index,
										face_index, time_index]
						last_pos = data[run_index, p0_index,
										face_index, time_index-1]
						if position[0] - last_pos[0] > size/2:
							data[run_index, p0_index, face_index,
									time_index, 0] -= size
						if last_pos[0] - position[0] > size/2:
							data[run_index, p0_index, face_index,
									time_index, 0] += size
						if position[1] - last_pos[1] > size/2:
							data[run_index, p0_index, face_index,
									time_index, 1] -= size
						if last_pos[1] - position[1] > size/2:
							data[run_index, p0_index, face_index,
									time_index, 1] += size
	return data

################################################################################

def plot_msd (data, run_values, p0_values, time_values):
	disp = data - data[:,:,:,0,np.newaxis,:]
	mean_disp = np.mean(disp, axis = 2)
	disp -= mean_disp[:,:,np.newaxis,:,:]
	square_disp = disp[:,:,:,:,0]**2 + disp[:,:,:,:,1]**2
	msd = np.mean(square_disp, axis = (0,2)) # average over runs and faces
	fig, ax = plt.subplots(1)#, figsize=(8, 6))
	cmap = plt.get_cmap('jet_r')
	lower = 10
	upper = 51
	for p0_index, p0_value in enumerate(p0_values):
		ax.plot(time_values[lower:upper], msd[p0_index,lower:upper],
				color = cmap(p0_index/(len(p0_values)-1)))
	ax.text(1.2e3,1.3e-2,'Increasing $p_0$',
			bbox = dict(facecolor='white', alpha=1,
						boxstyle='round', linestyle=''))
	ax.annotate("", xy=(3e3, 1.4), xytext=(3e3, 2e-2),
					arrowprops=dict(arrowstyle='->',
									color = 'white',
									linewidth=4),
					zorder = 9)
	ax.annotate("", xy=(3e3, 1.3), xytext=(3e3, 2e-2),
					arrowprops=dict(arrowstyle='->',
									color = 'black',
									linewidth=2),
					zorder = 10)
#	ax.arrow(3e3, 2e-2, 0, 0.9,
#				color = 'black',
#				zorder = 20,
#				width = 40,
#				head_width = 500,
#				head_length = 0.3)
	ax.set_xlabel('Time Steps')
	ax.set_ylabel('Mean Square Displacement')
	ax.grid(True)
	ax.set_xscale('log')
	ax.set_yscale('log')
	fig.tight_layout()
	plt.savefig('../plots/flat_msd.svg')
	plt.rc('pgf', texsystem='pdflatex')
	plt.savefig('../plots/flat_msd.pgf')
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
#	D_eff = (msd[:,-1] - msd[:,start])/(time_values[-1] - time_values[start])/4
#	D_eff /= 0.1**2/2/4.
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
						default = '../spp_flat_output',
						type = str,
						help = 'directory of simulation output data')
	args = parser.parse_args()
	out_dir = Path(args.datadir)
	run_values = np.zeros(0, dtype = int)
	for run_dir in out_dir.glob('run_*'):
		run_values = np.append(run_values, int(run_dir.name.split('_')[-1]))
	run_values = np.sort(run_values)
	run_dir = out_dir / (f'run_{run_values[0]:02d}')
	p0_values = np.zeros(0, dtype = float)
	for p0_dir in run_dir.glob('p0_*'):
		p0_values = np.append(p0_values, float(p0_dir.name.split('_')[-1]))
	p0_values = np.sort(p0_values)
	p0_dir = run_dir / (f'p0_{p0_values[0]:1.3f}')
	time_values = np.zeros(0, dtype = int)
	for time_file in p0_dir.glob('time_*.txt'):
		time_value = int((time_file.name.split('.')[0]).split('_')[-1])
		time_values = np.append(time_values, time_value)
	time_values = np.sort(time_values)
	data_file = out_dir/'data.npy'
	if data_file.exists():
		data = np.load(data_file)
	else:
		data = assemble_positions(out_dir, run_values, p0_values, time_values)
		np.save(data_file, data)
	msd = plot_msd(data, run_values, p0_values, time_values)
#	plot_Deff(msd, run_values, p0_values, time_values)

################################################################################
# EOF
