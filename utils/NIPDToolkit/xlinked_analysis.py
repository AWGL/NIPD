import pandas as pd
import numpy as np
import collections
from NIPDToolkit.nipd_tools import NIPDTools


class XLinkedAnalyser():
	
	def __init__(self,
		family_csv,
		maternal_sample,
		proband_sample,
		plasma_sample,
		gene,
		min_genotype_depth=30,
		min_distance_between_snps=200,
		min_snps_per_block=25,
		min_fetal_fraction=0.01):

		self.family_csv = family_csv
		self.maternal_sample = maternal_sample
		self.proband_sample = proband_sample
		self.plasma_sample = plasma_sample
		self.gene = gene
		self.min_genotype_depth = min_genotype_depth
		self.min_distance_between_snps = min_distance_between_snps
		self.min_snps_per_block = min_snps_per_block
		self.min_fetal_fraction = min_fetal_fraction
		
		self.df = None
		
		self.informative_snp_count = None
		self.xlinked_d_value = None
		self.xlinked_g_value = None
		self.xlinked_df_fwd = None
		self.xlinked_df_fwd_haplotype_blocks = None
		self.xlinked_df_rev = None
		self.xlinked_df_rev_haplotype_blocks = None
		self.xlinked_mean_snps_per_block = None
		
		
		self.initial_variant_count = None
		self.after_depth_filter_variant_count = None
		self.after_distance_filter_variant_count = None
		self.after_gene_filter_variant_count = None

		self.mean_snp_depth = None
		self.fetal_fraction = None
		
	
	
	def read_family_csv(self):
		"""
		Set the initial dataframe
		"""
		
		self.df = pd.read_csv(self.family_csv, sep='\t')
		self.df = self.df.sort_values(['CHROM', 'POS'])
		
		self.df = self.df[self.df['CHROM'] == 'X']

		self.df.sort_values(by='POS', ascending=True)
		
		#Get the initial variant count
		self.initial_variant_count = self.df.shape[0]


	def fix_column_names(self):
		"""
		Change column names to valid python variable names by replacing '.' with '_'

		"""
		
		# Loop through columns and change names
		fixed_columns = []

		for column_name in self.df.columns:

			fixed_column_name = column_name.replace('.', '_').replace('-', '_')

			if fixed_column_name not in ['CHROM', 'POS', 'REF', 'ALT','ID','GENE']:

				fixed_column_name = 'sample_' + fixed_column_name

			fixed_columns.append(fixed_column_name)
		
		# reset columns
		self.df.columns = fixed_columns
		
		#change sample names
		self.maternal_sample = self.maternal_sample.replace('-', '_')
		self.proband_sample = self.proband_sample.replace('-', '_')
		self.plasma_sample = self.plasma_sample.replace('-', '_')

		self.maternal_sample = 'sample_' + self.maternal_sample
		self.proband_sample = 'sample_' + self.proband_sample
		self.plasma_sample = 'sample_' + self.plasma_sample


	def get_average_snp_depth(self):
		"""
		Get the average SNP depth across all samples.
		"""
		
		df = self.df
	
		snp_depths =  df[self.maternal_sample + '_DP'] + \
			df[self.proband_sample + '_DP'] + \
			df[self.plasma_sample + '_DP']
		
		self.mean_snp_depth = snp_depths.mean()
	
		return self.mean_snp_depth

	def filter_on_depth(self):
		"""
		Drop any SNPs with a depth below the minimum in any of the samples.
		"""
		
		df = self.df
		
		df = df[(df[self.maternal_sample + '_DP'] > self.min_genotype_depth) &
		   (df[self.proband_sample + '_DP'] > self.min_genotype_depth) &
		   (df[self.plasma_sample + '_DP'] > self.min_genotype_depth)]
		
		df = df[(df[self.maternal_sample + '_GT'] != './.') &
		   (df[self.proband_sample + '_GT'] != './.') &
		   (df[self.plasma_sample + '_GT'] != './.')]
		
		self.df = df
		self.after_depth_filter_variant_count = self.df.shape[0]


	def filter_on_distance(self):
		"""
		Filter snps so there is a minimum distance between of a certain number of basepairs.
	
		Assumes df is sorted in order of position.
	
		The standard as per most SPRT/RHDO papers is a distance of 200.
	
		"""
		snps_to_keep = []

		last_position = 0

		for snp in self.df.itertuples():

			position =  snp.POS

			if (position - last_position) > self.min_distance_between_snps:

				snps_to_keep.append(snp)

				last_position = position


		self.df = pd.DataFrame(snps_to_keep)
		self.df.drop('Index', inplace=True, axis=1)
		self.after_distance_filter_variant_count = self.df.shape[0]


	def assign_analysis_gene(self):
		"""
		Select gene with which to analyse 
		"""
		self.df = self.df[self.df['GENE'] == self.gene]
		self.after_gene_filter_variant_count = self.df.shape[0]


	def genotype_samples(self):
		"""
		Are each of the samples HET or HOM?
		"""
		
		for sample in [self.maternal_sample, self.proband_sample, self.plasma_sample]:
		
			self.df[sample + '_call'] = self.df.apply(NIPDTools.genotype_sample, axis=1, sample_name = sample)

	def get_allele_depths(self):
		"""
		Get the allele depths for the ref and alt alleles for the ref and alts
		"""
		
		for sample in [self.maternal_sample, self.proband_sample, self.plasma_sample]:
			
			self.df[sample + '_AD_ref'] = self.df.apply(NIPDTools.get_allele_depth_ref, axis=1, sample_name = sample)
			self.df[sample + '_AD_alt'] = self.df.apply(NIPDTools.get_allele_depth_alt, axis=1, sample_name = sample)


	def get_relevant_snps(self):

		self.df = self.df[(self.df[self.maternal_sample + '_call'] == 'HET') & (self.df[self.proband_sample + '_call'] == 'HOM') ]
		self.informative_snp_count = self.df.shape[0]

	def get_haplotype_counts(self):

		self.df['hapA_count'] = self.df.apply(NIPDTools.get_hapA_counts_xlinked, axis=1, proband_sample=self.proband_sample, plasma_sample=self.plasma_sample)
		self.df['hapB_count'] = self.df.apply(NIPDTools.get_hapB_counts_xlinked, axis=1, proband_sample=self.proband_sample, plasma_sample=self.plasma_sample)
		self.df['hapA_count'] = self.df['hapA_count'].astype('int64')
		self.df['hapB_count'] = self.df['hapB_count'].astype('int64')

	def get_fetal_fraction_x_linked(self):
	    """
	    Use the formula from the Birmingham paper to calculate the fetal \
	    fraction for X-linked data.
	    
	    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4864947/
	    
	    cffDNA % = [(|cHapA – cHapB| * 2) / (cHapA + cHapB + |cHapA – cHapB|)] * 100
	    
	    """
	    cHapA = self.df['hapA_count'].sum()
	    cHapB = self.df['hapB_count'].sum()

	    fetal_fraction = (abs(cHapA-cHapB)*2) / (cHapA + cHapB + abs(cHapA-cHapB))
	    
	    self.fetal_fraction = fetal_fraction
	    return fetal_fraction

	def sprt_calculate_d_xlinked(self):
	    """
	    Calculate the d value needed for SPRT analysis:
	    
	    See the following papers for more information:
	    
	    1) https://www.ncbi.nlm.nih.gov/pubmed/21148127
	    
	    2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1934923/bin/pnas_0705765104_index.html
	    
	    The d value is used by the function sprt_perform_test()
	    
	    """
	    
	    # q0 is the allele imbalance if the null hypothesis is correct i.e. the fetus has the normal allele.
	    # the 0.5 value is because we would expect to see this much deviation from the 0.5 allele balance seen \
	    # in HET SNPs and then some more deviancy depending on the fetal fraction.
	    
	    
	    q0 = 0.5 - (self.fetal_fraction / 2)
	        
	    # q1 is the allele imbalance if the null hypothesis is incorrect i.e. the fetus has the mutant allele.

	    q1 = 0.5 + (self.fetal_fraction / 2)
	    
	    
	    d_value = (1 - q1) / (1 - q0)
	    self.xlinked_d_value = d_value
	    
	    return d_value


	def sprt_calculate_g_xlinked(self):
	    """
	    Calculate the g value needed for SPRT analysis:
	    
	    See the following papers for more information:
	    
	    1) https://www.ncbi.nlm.nih.gov/pubmed/21148127
	    
	    2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1934923/bin/pnas_0705765104_index.html
	    
	    The g value is used by the function sprt_perform_test()
	    
	    """
	    
	    # q0 is the allele imbalance if the null hypothesis is correct i.e. the fetus has the normal allele.
	    # the 0.5 value is because we would expect to see this much deviation from the 0.5 allele balance seen \
	    # in HET SNPs and then some more deviancy depending on the fetal fraction.
	    
	    q0 = 0.5 - (self.fetal_fraction/2)

	    # q1 is the allele imbalance if the null hypothesis is incorrect i.e. the fetus has the mutant allele.

	    q1 = 0.5 + (self.fetal_fraction/2)

	    g_value = ( q1*(1 - q0) ) / ( q0*(1 - q1) )

	    self.xlinked_g_value = g_value
	    
	    return g_value

	def analyse_xlinked_snps(self):

		self.sprt_calculate_d_xlinked()
		self.sprt_calculate_g_xlinked()

		self.xlinked_df_fwd = NIPDTools.sprt_perform_test(self.df, self.xlinked_d_value, self.xlinked_g_value, self.min_snps_per_block)
		self.xlinked_df_fwd_haplotype_blocks = NIPDTools.call_haplotype_blocks(self.xlinked_df_fwd)

		self.xlinked_df_rev = NIPDTools.sprt_perform_test(self.df.sort_values(by='POS', ascending=False), self.xlinked_d_value, self.xlinked_g_value, self.min_snps_per_block)
		self.xlinked_df_rev_haplotype_blocks = NIPDTools.call_haplotype_blocks(self.xlinked_df_rev)

	def run_analysis(self):

		self.read_family_csv()
		self.fix_column_names()
		self.get_average_snp_depth()
		self.filter_on_depth()
		self.genotype_samples()
		self.assign_analysis_gene()
		self.filter_on_distance()
		self.get_allele_depths()
		self.get_relevant_snps()
		self.get_haplotype_counts()
		self.get_fetal_fraction_x_linked()
		self.analyse_xlinked_snps()






