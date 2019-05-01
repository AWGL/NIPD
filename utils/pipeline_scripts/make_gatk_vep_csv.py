"""
A script for creating a CSV file from the VEP annotated GATK VCF

Usage: 

python utils/pipeline_scripts/split_joint_vcf_by_family.py 
		--input qfiltered_jointvcf_anno_norm_vep/190215_NB551415_0010_AH5772AFXY_all_chr_qfiltered_anno_norm_vep.vcf
		--config config/development_local.yaml
		--output_dir output/family_vcfs 
		--output_prefix 190215_NB551415_0010_AH5772AFXY_all_chr_qfiltered_anno_norm_vep
		--family_id FAM001

This script is NOT designed to be used outside of the NIPD pipeline. 

"""
import argparse
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet
import yaml

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

def filter_func(variant):

	family_members = variant.family.get_members()

	


parser = argparse.ArgumentParser(description='Split a VCF file by family')
parser.add_argument('--input', type=str, nargs=1, required=True,
				help='The input VCF file with many samples')
parser.add_argument('--config', type=str, nargs=1,required=True,
				help='The config file which describes the family structure.')
parser.add_argument('--output_dir', type=str, nargs=1,required=True,
				help='The output directory to put the output files in.')
parser.add_argument('--output_prefix', type=str, nargs=1,required=True,
				help='The output file prefix.')
args = parser.parse_args()

config_dict = parse_config(args.config[0])

# Get args
config = args.config[0]
vcf = args.input[0]



for family_id in config_dict['families'].keys():

	# Get family members for provided family id
	family_info = config_dict['families'][family_id]
	members = family_info['members']

	# Create a family object - affected and sex are not important we are just converting to dataframe.
	my_family = Family(family_id)

	for member in members:

		sample_id = members[member]

		if sample_id != None:
			new_family_member = FamilyMember(sample_id, family_id, 1,False  )
			my_family.add_family_member(new_family_member)

	# Create variant set
	my_variant_set = VariantSet()
	my_variant_set.add_family(my_family)

	# Read vcf
	my_variant_set.read_variants_from_vcf(vcf, proband_variants_only=False )

	# Convert to dataframe
	df = my_variant_set.to_df(add_inheritance=False)

	# Get columns to keep
	columns = ['variant_id',
	  'csq_Consequence',
	  'filter_status',
	  'csq_SYMBOL',
	  'csq_Feature',
	  'csq_HGVSc',
	  'csq_HGVSp',
	  'worst_consequence',
	  'csq_CLIN_SIG',
	  'csq_gnomADg_AF_POPMAX',
	  'csq_gnomADe_AF_POPMAX',
	  'csq_Existing_variation',
	  'csq_IMPACT',
	  'csq_PICK']

	for member in members:
	    
	    sample_id = members[member]

	    if sample_id != None:
	    
	    	columns.append(f'{sample_id}_GT')
	    	columns.append(f'{sample_id}_AD')
	    	columns.append(f'{sample_id}_DP')
	    	columns.append(f'{sample_id}_GQ')


	df[columns].to_csv(f'{args.output_dir[0]}/{args.output_prefix[0]}_{family_id}.csv', index=False)