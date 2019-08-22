"""
A script for fixing the column headers of a CSV file made from a Platypus VCF.

We need a depth (DP) column and allele depths (AD) field for each sample in the CSV.

Platypus does not create a DP or AD field so we create these from other FORMAT fields availible.

Usage: 

python utils/pipeline_scripts/split_joint_vcf_by_family.py 
		--input input.csv
		--output output.csv

This script is NOT designed to be used outside of the NIPD pipeline. 

"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Split a VCF file by family')
parser.add_argument('--input', type=str, nargs=1, required=True,
				help='The input CSV file.')
parser.add_argument('--output', type=str, nargs=1,required=True,
				help='The output CSV file.')

args = parser.parse_args()

# Useful functions

def get_depth(df, sample_id):
	"""
	Create a columns sample_id_DP from the NR format field.
	"""
	total_reads = df[f'{sample_id}.NR']
	
	return total_reads

def get_allele_depths(df, sample_id):
	"""
	Create an AD field from the NR and NV format fields.
	"""
	
	total_reads = df[f'{sample_id}.NR']
	
	alt_reads = df[f'{sample_id}.NV']
	
	ref_reads = total_reads - alt_reads
	
	return f'{ref_reads},{alt_reads}'

# open csv
df = pd.read_csv(args.input[0], sep='\t')

# get all columns as list
columns = list(df.columns)

# remove fields not associated with a sample e.g non FOMMAT fields
columns.remove('CHROM')
columns.remove('POS')
columns.remove('REF')
columns.remove('ALT')
columns.remove('ID')
columns.remove('GENE')

# get the sample id from remaining columns
columns = [x.split('.')[0] for x in columns]

# get unique sample ids
columns = list(set(columns))

# apply functions to get the extra columns
for sample in columns:
	
	df[f'{sample}.DP'] = df.apply(get_depth, axis=1, args=(sample,))
	df[f'{sample}.AD'] = df.apply(get_allele_depths, axis=1, args=(sample,))

# save to csv

df.to_csv(args.output[0], sep='\t', index=False)