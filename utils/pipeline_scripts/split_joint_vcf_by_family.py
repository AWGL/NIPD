"""
A script for splitting a VCF file by family.

Would usually use a PED file for this, but NIPD samples do not fit into the PED file structure.

Usage: 

"""

import argparse
import yaml
import os
from pathlib import Path 

def parse_config(config_location):
	"""
	Parse the YAML config file.
	"""

	with open(config_location, 'r') as stream:

		try:
			return yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
			raise


parser = argparse.ArgumentParser(description='Split a VCF file by family')
parser.add_argument('--input', type=str, nargs=1, required=True,
				help='The input VCF file with many samples')
parser.add_argument('--config', type=str, nargs=1,required=True,
				help='The config file which describes the family structure.')
parser.add_argument('--ref', type=str, nargs=1,required=True,
				help='The reference genome location.')
parser.add_argument('--output_dir', type=str, nargs=1,required=True,
				help='The output directory to put the output files in.')
args = parser.parse_args()


config_dict = parse_config(args.config[0])

# Get worksheet ID from input
worksheet_id = Path(args.input[0]).name
worksheet_id = worksheet_id.split('_')[0]


# Get relevant sample names for each family
for family in config_dict['families'].keys():

	family_samples =[]

	# Loop though each member of the family e.g. mum, dad
	for family_member in config_dict['families'][family]['members'].keys():

		# If the sample name is not None
		if config_dict['families'][family]['members'][family_member] is not None:

			# Append the sample name to the list
			family_samples.append(config_dict['families'][family]['members'][family_member])

	#TO DO - Add checks to see if we have the minimum correct samples for each type of analysis
	#For example do we have mum, plasma and proband for X linked analysis.

	sample_arguments = ('-sn ' +  ' -sn '.join(family_samples))

	output_file = '{output_dir}/{worksheet}_qfiltered_selected_{FAMID}.vcf'.format(
		output_dir = args.output_dir[0],
		worksheet=worksheet_id,
		FAMID = family
		)
	
	command = 'gatk SelectVariants -R {ref} -V {input} -O {output} {samples}'.format(ref=args.ref[0], input=args.input[0], output=output_file, samples=sample_arguments)


	os.system(command)