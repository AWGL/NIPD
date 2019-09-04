"""
A script for creating a calculating the sex of a sample from a VCF.

Uses the mean heterozygosity on X chromosome.

Only works on platypus VCFs.

Prints to stdout

Usage: 

python utils/pipeline_scripts/calculate_sex.py 
		--input input.vcf
		--sample_id test_sample
		--female_cutoff 0.2
		--min_variants 100
		--min_depth 10
		--min_gq 20

Options female_cutoff, min_variants, min_depth and min_gq are optional and sensible defaults are supplied.

This script is NOT designed to be used outside of the NIPD pipeline. 

"""

from pysam import VariantFile
import numpy as np
import argparse


parser = argparse.ArgumentParser(description='Calculate the sex of a sample in a VCF file')
parser.add_argument('--input', type=str, nargs=1, required=True,
				help='The input Platypus VCF file')
parser.add_argument('--sample_ids', type=str, nargs='+',required=True,
				help='The sample ids to analyse')

parser.add_argument('--min_variants', type=int, nargs=1,required=False,
				help='The min number of calls needed to make a sex determination.')
parser.add_argument('--female_cutoff', type=float, nargs=1,required=False,
				help='samples with more than this level of heterozygosity on X chrom will be called as Female.')
parser.add_argument('--min_depth', type=str, nargs=1,required=False,
				help='The minimum depth required to use a variant.')
parser.add_argument('--min_gq', type=str, nargs=1,required=False,
				help='The minimum GQ required to use a variant.')

args = parser.parse_args()

input_vcf = args.input[0]
sample_ids = args.sample_ids
min_variants = args.min_variants
female_cutoff = args.female_cutoff
min_depth = args.min_depth
min_gq = args.min_gq

if min_variants == None:

	min_variants = 75

if female_cutoff == None:

	female_cutoff = 0.2

if min_depth == None:

	min_depth = 10

if min_gq == None:

	min_gq = 20

# print header
print ('sample_id\tcalculated_sex\tmean_x_heterozygosity')

# loop through each sample
for sample_id in sample_ids:

	het_array = []

	bcf_in = VariantFile(input_vcf)

	# Loop through X chromosome
	for rec in bcf_in.fetch('X'):
		
		ref = rec.ref
		alt = rec.alts
		alt = alt[0]
		ref_and_alt = [ref, alt]


		
		# Get the genotype of a sample
		sample_genotype_data = rec.samples[sample_id]

		gts =[]

		for allele in sample_genotype_data['GT']:

			if allele == None:

				gts.append('.')

			else:

				gts.append(ref_and_alt[allele])
		
		# get the filter status of the variant
		filter_status = rec.filter.keys()
		
		# get the depth of the variant
		depth = sample_genotype_data['NR']
		depth = depth[0]
		
		# get the GQ of the variant
		gq = sample_genotype_data['GQ']
		gq = gq[0]

		# skip rows where we don't have value for either of these
		if depth == None or gq == None:

			break
		
		# if the variant passese QC
		if filter_status[0] == 'PASS' and depth > min_depth and gq > min_gq: 
			
			# don't add missing variants
			if gts[0] == '.' and gts[1] == '.':
				
				pass
			
			# if the variant is het add a 1 to the array else 0 for homs
			elif gts.count(alt) == 1:
				
				het_array.append(1)
				
			else:
				
				het_array.append(0)

	# unknown if we have less than the minimum variants
	if len(het_array) < min_variants:

		calculated_sex = 'Unknown'
		mean_heterozygosity = 'NA'

	else:

		mean_heterozygosity = round(np.array(het_array).mean(), 4)

		if mean_heterozygosity > female_cutoff:

			calculated_sex = 'Female'

		else:

			calculated_sex = 'Male'

	# print results to screen
	print (f'{sample_id}\t{calculated_sex}\t{mean_heterozygosity}')





	

