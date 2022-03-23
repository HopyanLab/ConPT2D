#!/usr/bin/env /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

################################################################################
#===============================================================================
# plot_energy.py
#===============================================================================
################################################################################

def plot_energy (energy_file = Path('../flat_output/edge_energy.csv'),
				 output_file = Path('../plots/energy_plot.svg')):
	energy_data = np.loadtxt(energy_file, comments = '#', delimiter = ',')
	lengths = -energy_data[:,0] / energy_data[0,0]
	energies = energy_data[:,1]
	energies -= min(energies)
	plt.plot(lengths, energies,'-')
	plt.grid(True)
	plt.title('Energy Barrier versus Edge Length')
	plt.xlabel('Normalised Edge Length')
	plt.ylabel('System Energy ($\Delta \epsilon$)')
	plt.savefig(output_file)

################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'Generate periodic Voronoi diagram')
	parser.add_argument('-o', '--outfile',
						nargs = 1,
						default = '../plots/energy_plot.svg',
						type = str,
						required = False,
						help = 'file to put plot in')
	parser.add_argument('energyfile',
						nargs = '?',
						default = '../flat_output/edge_energy.csv',
						type = str,
						help = 'file with length vs energy data')
	args = parser.parse_args()
	plot_energy(energy_file = Path(args.energyfile),
				output_file = Path(args.outfile))

################################################################################
# EOF
