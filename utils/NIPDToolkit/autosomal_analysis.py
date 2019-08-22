import pandas as pd
import numpy as np
import collections
from NIPDToolkit.nipd_tools import NIPDTools

np.seterr(divide='ignore', invalid='ignore')

class AutosomalAnalyser():
	
	def __init__(self,
		family_csv,
		maternal_sample,
		paternal_sample,
		proband_sample,
		plasma_sample,
		gene,
		min_genotype_depth=30,
		min_distance_between_snps=200,
		min_snps_per_block=25,
		min_ks_p_value=0.001,
		min_fetal_fraction=0.01):

		self.family_csv = family_csv
		self.maternal_sample = maternal_sample
		self.paternal_sample = paternal_sample
		self.proband_sample = proband_sample
		self.plasma_sample = plasma_sample
		self.gene = gene
		self.min_genotype_depth = min_genotype_depth
		self.min_distance_between_snps = min_distance_between_snps
		self.min_snps_per_block = min_snps_per_block
		self.min_ks_p_value = min_ks_p_value
		self.min_fetal_fraction = min_fetal_fraction
		
		self.df = None
		
		self.type4a_df = None
		self.type4a_variant_count = None
		self.type4a_g_value = None
		self.type4a_d_value = None
		self.type4a_df_fwd = None
		self.type4a_df_fwd_haplotype_blocks = None
		self.type4a_df_rev = None
		self.type4a_df_rev_haplotype_blocks = None
		self.type4a_mean_snps_per_block = None
		
		self.type4b_df = None
		self.type4b_variant_count = None
		self.type4b_g_value = None
		self.type4b_d_value = None
		self.type4b_df_fwd = None
		self.type4b_df_fwd_haplotype_blocks = None
		self.type4b_df_rev = None
		self.type4b_df_rev_haplotype_blocks = None
		self.type4b_mean_snps_per_block = None
		
		self.type3_df = None
		self.type3_variant_count = None
		self.type3_df_fwd = None
		self.type3_df_fwd_haplotype_blocks = None
		self.type3_df_rev = None
		self.type3_df_rev_haplotype_blocks = None
		self.type3_mean_snps_per_block = None
		
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
		
		self.df = self.df[self.df['CHROM'] != 'X']
		
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
		self.paternal_sample = self.paternal_sample.replace('-', '_')
		self.proband_sample = self.proband_sample.replace('-', '_')
		self.plasma_sample = self.plasma_sample.replace('-', '_')

		self.maternal_sample = 'sample_' + self.maternal_sample
		self.paternal_sample = 'sample_' + self.paternal_sample
		self.proband_sample = 'sample_' + self.proband_sample
		self.plasma_sample = 'sample_' + self.plasma_sample
		
	
	def get_average_snp_depth(self):
		"""
		Get the average SNP depth across all samples.
		"""
		
		df = self.df
	
		snp_depths =  np.append(df[self.maternal_sample + '_DP'], [ \
			df[self.paternal_sample + '_DP'], \
			df[self.proband_sample + '_DP'], \
			df[self.plasma_sample + '_DP']])
		
		self.mean_snp_depth = snp_depths.mean()
	
		return self.mean_snp_depth
	
	def filter_on_depth(self):
		"""
		Drop any SNPs with a depth below the minimum in any of the samples.
		"""
		
		df = self.df
		
		df = df[(df[self.maternal_sample + '_DP'] > self.min_genotype_depth) &
		   (df[self.paternal_sample + '_DP'] > self.min_genotype_depth) &
		   (df[self.proband_sample + '_DP'] > self.min_genotype_depth) &
		   (df[self.plasma_sample + '_DP'] > self.min_genotype_depth)]
		
		df = df[(df[self.maternal_sample + '_GT'] != './.') &
		   (df[self.paternal_sample + '_GT'] != './.') &
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
		self.df.drop('Index', inplace=True,axis=1)
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
		
		for sample in [self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample]:
		
			self.df[sample + '_call'] = self.df.apply(NIPDTools.genotype_sample, axis=1, sample_name = sample)
			
	def get_alelle_depths(self):
		"""
		Get the allele depths for the ref and alt alleles for the ref and alts
		"""
		
		for sample in [self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample]:
			
			self.df[sample + '_AD_ref'] = self.df.apply(NIPDTools.get_allele_depth_ref, axis=1, sample_name = sample)
			self.df[sample + '_AD_alt'] = self.df.apply(NIPDTools.get_allele_depth_alt, axis=1, sample_name = sample)
			 
	def assign_snp_types(self):
		"""
		Assign each SNP as type from lo et al 2010.
		"""
		self.df['type'] = self.df.apply(NIPDTools.get_snp_type, axis=1, args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample,))
		
	def get_fetal_fraction(self):
		"""
		Calculate the fetal fraction on type 1 SNPs
		"""
		
		type1_snps = self.df[self.df['type'] == 'type_1'].copy(deep=True)
		
		type1_snps['autosomal_ff'] = type1_snps.apply(NIPDTools.calculate_autosomal_fetal_fraction, axis=1, args=(self.maternal_sample, self.paternal_sample, self.plasma_sample,))
		
		fetal_fraction = type1_snps['autosomal_ff'].mean()
		
		self.fetal_fraction = fetal_fraction
		
		if self.fetal_fraction < self.min_fetal_fraction:
			
			raise ValueError(f'The fetal fraction was below the minimum ({self.min_fetal_fraction})')
		
		return fetal_fraction
	
	def sprt_calculate_d_type4a(self):
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
		
		q0 = 0.5 

		# q1 is the allele imbalance if the null hypothesis is incorrect i.e. the fetus has the mutant allele.

		q1 = 0.5 + (self.fetal_fraction / 2)

		d_value = (1 - q1) / (1 - q0)
		
		self.type4a_d_value = d_value

		return d_value
	
	def sprt_calculate_g_type4a(self):
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
		

		q0 = 0.5

		# q1 is the allele imbalance if the null hypothesis is incorrect i.e. the fetus has the mutant allele.

		q1 = 0.5 + (self.fetal_fraction/2)

		g_value = ( q1*(1 - q0) ) / ( q0*(1 - q1) )
		
		self.type4a_g_value = g_value

		return g_value
	
	
	
	def sprt_calculate_d_type4b(self):
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

		q1 = 0.5 

		d_value = (1 - q1) / (1 - q0)
		
		self.type4b_d_value = d_value

		return d_value

	def sprt_calculate_g_type4b(self):
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
		
		q0 = 0.5 - (self.fetal_fraction / 2)

		# q1 is the allele imbalance if the null hypothesis is incorrect i.e. the fetus has the mutant allele.

		q1 = 0.5 

		g_value = ( q1*(1 - q0) ) / ( q0*(1 - q1) )
		
		
		self.type4b_g_value = g_value

		return g_value 
	
	
	def get_type4a_mean_snps_per_haplotype_block(self):
		"""
		Get the number of SNPs per Haplotype Block
		"""
		
		if self.type4a_df_fwd_haplotype_blocks != None and self.type4a_df_rev_haplotype_blocks != None:
			
			all_blocks = [x[3] for x in self.type4a_df_fwd_haplotype_blocks] + [x[3] for x in self.type4a_df_rev_haplotype_blocks]
			
			mean_snps_per_block = np.mean(all_blocks)
			
			return mean_snps_per_block
			
		else:
			
			raise ValueError('There are no typea haplotype blocks to process')
			
	def get_type4a_haplotype_block_counts(self):
		"""
		Return a count of the haplotype blocks
		"""
		
		if self.type4a_df_fwd_haplotype_blocks != None and self.type4a_df_rev_haplotype_blocks != None:
		
			haplotype_blocks = [x[4] for x in self.type4a_df_fwd_haplotype_blocks] + [x[4] for x in self.type4a_df_rev_haplotype_blocks]
		
			return collections.Counter(haplotype_blocks)
		
		else:
			
			 raise ValueError('There are no typea haplotype blocks to process') 
			
		
	def get_type4b_mean_snps_per_haplotype_block(self):
		"""
		Get the number of SNPs per Haplotype Block
		"""
		
		if self.type4b_df_fwd_haplotype_blocks != None and self.type4b_df_rev_haplotype_blocks != None:
			
			all_blocks = [x[3] for x in self.type4b_df_fwd_haplotype_blocks] + [x[3] for x in self.type4b_df_rev_haplotype_blocks]
			
			mean_snps_per_block = np.mean(all_blocks)
			
			return mean_snps_per_block
			
		else:
			
			raise ValueError('There are no typea haplotype blocks to process')
			
			
	def get_type4b_haplotype_block_counts(self):
		"""
		Return a count of the haplotype blocks
		"""
		
		if self.type4b_df_fwd_haplotype_blocks != None and self.type4b_df_rev_haplotype_blocks != None:
		
			haplotype_blocks = [x[4] for x in self.type4b_df_fwd_haplotype_blocks] + [x[4] for x in self.type4b_df_rev_haplotype_blocks]
		
			return collections.Counter(haplotype_blocks)
		
		else:
			
			 raise ValueError('There are no typea haplotype blocks to process') 

	def get_type3_mean_snps_per_haplotype_block(self):
		"""
		Get the number of SNPs per Haplotype Block
		"""
		
		if self.type3_df_fwd_haplotype_blocks != None and self.type3_df_rev_haplotype_blocks != None:
			
			all_blocks = [x[3] for x in self.type3_df_fwd_haplotype_blocks] + [x[3] for x in self.type3_df_rev_haplotype_blocks]
			
			mean_snps_per_block = np.mean(all_blocks)
			
			return mean_snps_per_block
			
		else:
			
			raise ValueError('There are no typea haplotype blocks to process')
	
	def get_type3_haplotype_block_counts(self):
		"""
		Return a count of the haplotype blocks
		"""
		
		if self.type3_df_fwd_haplotype_blocks != None and self.type3_df_rev_haplotype_blocks != None:
		
			haplotype_blocks = [x[4] for x in self.type3_df_fwd_haplotype_blocks] + [x[4] for x in self.type3_df_rev_haplotype_blocks]
		
			return collections.Counter(haplotype_blocks)
		
		else:
			
			 raise ValueError('There are no typea haplotype blocks to process') 
				
	def analyse_type4a_snps(self):
		"""
		Perfrom SPRT analysis on Type4 a SNPs
		
		"""
		
		self.type4a_df = self.df[self.df['type'] == 'type_4a'].copy(deep=True)
		
		self.type4a_variant_count = self.type4a_df.shape[0]
		
		self.type4a_df['hapA_count'] = self.type4a_df.apply(NIPDTools.get_hapA_counts_type_4a, axis=1, args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample,))
		self.type4a_df['hapB_count'] = self.type4a_df.apply(NIPDTools.get_hapB_counts_type_4a, axis=1, args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample,))
		
		self.sprt_calculate_d_type4a()
		self.sprt_calculate_g_type4a()
		
		self.type4a_df_fwd = NIPDTools.sprt_perform_test(self.type4a_df.copy(deep=True), self.type4a_d_value, self.type4a_g_value, self.min_snps_per_block)
		
		self.type4a_df_fwd_haplotype_blocks  = NIPDTools.call_haplotype_blocks(self.type4a_df_fwd)
		
		
		self.type4a_df_rev = NIPDTools.sprt_perform_test(self.type4a_df.sort_values(by='POS', ascending=False).copy(deep=True), self.type4a_d_value, self.type4a_g_value, self.min_snps_per_block)
		
		self.type4a_df_rev_haplotype_blocks  = NIPDTools.call_haplotype_blocks(self.type4a_df_rev)
		
		
	def analyse_type4b_snps(self):
		"""
		Perform SPRT analysis on Type 4 b SNPs
		
		"""
		
		self.type4b_df = self.df[self.df['type'] == 'type_4b'].copy(deep=True)
		
		self.type4b_variant_count = self.type4b_df.shape[0]
		
		self.type4b_df['hapA_count'] = self.type4b_df.apply(NIPDTools.get_hapA_counts_type_4b, axis=1, args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample,))
		self.type4b_df['hapB_count'] = self.type4b_df.apply(NIPDTools.get_hapB_counts_type_4b, axis=1, args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample,))
		
		self.sprt_calculate_d_type4b()
		self.sprt_calculate_g_type4b()
		
		self.type4b_df_fwd = NIPDTools.sprt_perform_test(self.type4b_df.copy(deep=True), self.type4b_d_value, self.type4b_g_value, self.min_snps_per_block)
		
		self.type4b_df_fwd_haplotype_blocks = NIPDTools.call_haplotype_blocks(self.type4b_df_fwd)
		
		
		self.type4b_df_rev = NIPDTools.sprt_perform_test(self.type4b_df.sort_values(by='POS', ascending=False).copy(deep=True), self.type4b_d_value, self.type4b_g_value, self.min_snps_per_block)
		
		self.type4b_df_rev_haplotype_blocks  = NIPDTools.call_haplotype_blocks(self.type4b_df_rev)
		
		
	def analyse_type3_snps(self):
		"""
		Perform KS testing on type 3 SNPs
		"""
		self.type3_df = self.df[self.df['type'] == 'type_3'].copy(deep=True)
		
		self.type3_variant_count = self.type3_df.shape[0]
		
		self.type3_df['HapA_count'] = self.type3_df.apply(NIPDTools.get_hapa_pat_count, axis=1,args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample))
		self.type3_df['HapB_count'] = self.type3_df.apply(NIPDTools.get_hapb_pat_count, axis=1,args=(self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample))
		
		self.type3_df_fwd = NIPDTools.perform_ks_test(self.type3_df.copy(deep=True), self.min_ks_p_value, self.min_snps_per_block)
		self.type3_df_fwd_haplotype_blocks = NIPDTools.call_haplotype_blocks(self.type3_df_fwd)
		
		self.type3_df_rev = NIPDTools.perform_ks_test(self.type3_df.sort_values(by='POS', ascending=False).copy(deep=True), self.min_ks_p_value, self.min_snps_per_block)
		self.type3_df_rev_haplotype_blocks = NIPDTools.call_haplotype_blocks(self.type3_df_rev)
		
		
	def run_analysis(self):
		
		self.read_family_csv()
		self.fix_column_names()
		self.get_average_snp_depth()
		self.filter_on_depth()
		self.filter_on_distance()
		self.assign_analysis_gene()
		self.genotype_samples()
		self.get_alelle_depths()
		self.assign_snp_types()
		self.get_fetal_fraction()
		self.analyse_type4a_snps()
		self.analyse_type4b_snps()
		self.analyse_type3_snps()
		
	def downsample(self, downsample_rate):
		"""
		For testing - downsample the counts in the dataframe
		"""
		
		for sample in [self.maternal_sample, self.paternal_sample, self.proband_sample, self.plasma_sample]:
			
			self.df[sample + '_AD_ref'] = np.round(self.df[sample + '_AD_ref'].astype(int) * downsample_rate)
			self.df[sample + '_AD_alt'] = np.round(self.df[sample + '_AD_alt'].astype(int) * downsample_rate)
			self.df[sample + '_DP'] = np.round(self.df[sample + '_DP'].astype(int) * downsample_rate)
		
			self.df[sample + '_AD_ref'] = self.df[sample + '_AD_ref'].astype(int)
			self.df[sample + '_AD_alt'] = self.df[sample + '_AD_alt'].astype(int)
			self.df[sample + '_DP'] = self.df[sample + '_DP'].astype(int)
			
	def run_down_sample_analysis(self, downsample_rate):
		"""
		For testing - run an analysis with downsampling enabled
		"""
		
		self.read_family_csv()
		self.fix_column_names()
		self.filter_on_depth()
		self.filter_on_distance()
		self.assign_analysis_gene()
		self.genotype_samples()
		self.get_alelle_depths()
		self.downsample(downsample_rate)
		self.get_average_snp_depth()
		self.assign_snp_types()
		self.get_fetal_fraction()
		self.analyse_type4a_snps()
		self.analyse_type4b_snps()
		self.analyse_type3_snps()    