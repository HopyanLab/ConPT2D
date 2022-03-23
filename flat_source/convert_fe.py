#!/usr/bin/env /usr/bin/python3
import argparse
from pathlib import Path

################################################################################
#===============================================================================
# convert_fe.py
#===============================================================================
################################################################################

if __name__ == '__main__':
	# Use argparse to get arguements from commandline call
	parser = argparse.ArgumentParser(
							description = 'convert Surface Evolver fe file')
	# fe_filename = int(sys.argv[1])
	parser.add_argument('fe_filename',
						nargs = 1,
						default = ['../flat_output/initial_state.fe'],
						type = str,
						help = 'fe file to convert')
	parser.add_argument('-v', '--vertfile',
						nargs = 1,
						default = ['../flat_output/vertices.txt'],
						type = str,
						required = False,
						help = 'file to put vertex data in')
	parser.add_argument('-e', '--edgefile',
						nargs = 1,
						default = ['../flat_output/edges.txt'],
						type = str,
						required = False,
						help = 'file to put edge data in')
	args = parser.parse_args()
	vertfile = Path(args.vertfile[0])
	vertfile.parent.mkdir(exist_ok = True)
	edgefile = Path(args.edgefile[0])
	edgefile.parent.mkdir(exist_ok = True)
	fe_file = Path(args.fe_filename[0])
	output_type = 'none'
	size = 0
	with open(fe_file) as instream, \
			open(vertfile, 'w+') as vertstream, \
			open(edgefile, 'w+') as edgestream:
		for line in instream:
			if line == 'PERIODS\n':
				output_type = 'period'
			elif output_type == 'period' and size == 0:
				size = float(line.split()[0])
			elif line == 'vertices        /*  coordinates  */    \n':
				output_type = 'vertices'
			elif line == 'edges  \n':
				output_type = 'edges'
			elif line.startswith('faces'):
				output_type = 'none'
			elif output_type == 'vertices' and len(line.split()) >= 3:
				vertstream.write('\t'.join(
						[str(float(line.split()[1])/size),
						 str(float(line.split()[2])/size)]
								) + '\n')
			elif output_type == 'edges':
				edgestream.write('\t'.join(line.split()[1:5]) + '\n')

################################################################################
# EOF
