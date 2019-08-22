import pandas as pd
import numpy as np
from scipy import stats
import statistics

np.seterr(divide='ignore', invalid='ignore')

class NIPDTools():
	"""
	A class to store functions for analysing NIPD data.
	
	These don't really fit as methods to a specific class.
	
	"""
	
	def get_genotypes_as_list(gt):
		"""
		Split a genotype like 1/0 or 1/1 into a list.
		
		For example: 1/1 to [1,1]
		
		"""
		
		if '/' in gt:

			gt = gt.split('/')

		elif '|' in gt:

			gt = gt.split('|')

		else:

			raise ValueError('Found a wierd genotype')

		return gt
	
	def genotype_sample(df, sample_name):
		"""
		Return whether a sample is HOM or HET

		Add a NA if sample is ./. or 1/. or something else wierd
		"""

		gt = df[sample_name+'_GT']

		gt = NIPDTools.get_genotypes_as_list(gt)

		if (gt[0] not in 'ATGC') or (gt[1] not in 'ATGC'):
			
			print (gt, sample_name, df['POS'])

			raise ValueError('Genotype is not valid - Not ATGC')

		if gt[0] == gt[1]:

			return 'HOM'
		else:
			return 'HET'
		
	def get_allele_depth_ref(df, sample_name):
		"""
		For a given sample split the AD column into AD for ref and alts
		"""
	
		ad = df[sample_name + '_AD']

		ad = ad.split(',')

		return ad[0]
	
	def get_allele_depth_alt(df, sample_name):
		"""
		For a given sample split the AD column into AD for ref and alts
		"""
		ad = df[sample_name + '_AD']

		ad = ad.split(',')

		return ad[1]
	
	def get_snp_type(df, maternal_sample, paternal_sample, proband_sample, plasma_sample):
		"""
		Get the SNP type. See Lo et al - https://www.ncbi.nlm.nih.gov/pubmed/21148127
		
		Type 1: Dad = hom and Mum = hom for different alleles.
		Type 2: Dad = hom and Mum = hom for same allele.
		Type 3: Dad = het and Mum = hom
		Type 4: Dad = hom and Mum = het
		Type 5: Dad = het and Mum = het
		
		Type 4 SNPs are then broken down into two further categories:
		
			Type 4a: The paternal SNP allele matches the maternal SNP allele on the affected haplotype.
			Type 4b: The paternal SNP allele matched the maternal SNP allele on the normal haplotype.
			
		"""
		
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# Type 1 SNPs
		if (df[maternal_sample + '_call'] == 'HOM') and (df[paternal_sample + '_call'] == 'HOM') and (df[maternal_sample +'_GT'] != df[paternal_sample +'_GT'] ):

			return 'type_1'
		
		# Type 2 SNPs
		elif (df[maternal_sample + '_call'] == 'HOM') and (df[paternal_sample + '_call'] == 'HOM') and (df[maternal_sample +'_GT'] == df[paternal_sample +'_GT']):

			return 'type_2'
		
		# Type 3 SNPs
		elif (df[maternal_sample + '_call'] == 'HOM') and (df[paternal_sample + '_call'] == 'HET'):

			return 'type_3'
		
		# Type 4 SNPs
		elif df[maternal_sample + '_call'] == 'HET' and (df[paternal_sample + '_call'] == 'HOM'):
			
			
			hap_from_mum = None
			
			# If proband is het
			if proband_gts[0] != proband_gts[1]:
				
				# then the haplotype is the one not in dad
				hap_from_mum = set(mat_gts) - set(pat_gts)

				hap_from_mum = list(hap_from_mum)[0]

			else:
				
				# otherwise it is the same at the paternal allele.
				hap_from_mum = pat_gts[0]

			# If the paternal allele is the same as those on the affected proband then type a
			if hap_from_mum == pat_gts[0]:

				return 'type_4a'

			# otherwise it is a type b
			else:

				return 'type_4b'
		
		# Type 5 SNPs
		elif df[maternal_sample + '_call'] == 'HET' and (df[paternal_sample + '_call'] == 'HET'):

			return 'type_5'
		
		else:
			
			return 'NA'
		
		
	def calculate_autosomal_fetal_fraction(df, maternal_sample, paternal_sample, plasma_sample):
		"""
		Calculates the fetal fraction based on snps which the mum and \
		dad are both homozygous for different alleles.
		
		Input should be a dataframe containing autosomal type 1 snps

		See the Supplementary information in the below link:

		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5386415/

		"""

		# Get the genoypes of relevant samples
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		allele_from_dad = pat_gts[0]

		# Is the allele from dad the ref or alt?

		dad_ref_or_alt = ''

		if allele_from_dad == df['REF']:

			dad_ref_or_alt = 'REF'

		elif allele_from_dad == df['ALT']:

			dad_ref_or_alt = 'ALT'

		allele_from_mum = mat_gts[0]

		# Is the allele from dad the ref or alt?

		mum_ref_or_alt = ''

		if allele_from_mum == df['REF']:

			mum_ref_or_alt = 'REF'

		elif allele_from_mum== df['ALT']:

			mum_ref_or_alt = 'ALT'        


		# How many of the dad's alleles do we have in the plasma sample?               
		if dad_ref_or_alt == 'REF':

			read_count_from_dad = df[plasma_sample + '_AD_ref']

		elif dad_ref_or_alt == 'ALT':

			read_count_from_dad = df[plasma_sample + '_AD_alt']

		# How many of mums alleles do we have in the plasma sample?          
		if mum_ref_or_alt == 'REF':

			read_count_from_mum = df[plasma_sample + '_AD_ref']

		elif mum_ref_or_alt == 'ALT':

			read_count_from_mum = df[plasma_sample + '_AD_alt']

		return (2*int(read_count_from_dad)) / (int(read_count_from_dad) + int(read_count_from_mum))
	
	def get_hapA_counts_type_4a(df, maternal_sample, paternal_sample, proband_sample, plasma_sample ):
		"""
		Count the number of reads in the plasma sample which are the same as the allele that the \
		affected proband got from mum.

		For type 4a samples this will be the allele which is the same as the paternal allele \

		i.e. the allele which the paternal sample is hom for.

		"""
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# allele which dad is hom for
		hapa_allele = pat_gts[0]


		#is that allele the ref or alt?
		affected_haplotype = ''

		if hapa_allele == df['REF']:

			affected_haplotype = 'REF'

		elif hapa_allele == df['ALT']:

			affected_haplotype = 'ALT'

		 # if the affected haplotype is the REF then get the ref AD for the plasma sample
		if affected_haplotype == 'REF':

			return df[plasma_sample+'_AD_ref']

		elif affected_haplotype == 'ALT':

			return df[plasma_sample+'_AD_alt']
		
		
	def get_hapB_counts_type_4a(df, maternal_sample, paternal_sample, proband_sample, plasma_sample ):
		"""
		Count the number of reads in the plasma sample which are the same as the allele that the \
		affected proband got from mum.

		For type 4a samples this will be the allele which is the same as the paternal allele \

		i.e. the allele which the paternal sample is hom for.

		"""
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# allele which dad is hom for
		hapa_allele = pat_gts[0]


		#is that allele the ref or alt?
		affected_haplotype = ''

		if hapa_allele == df['REF']:

			affected_haplotype = 'REF'

		elif hapa_allele == df['ALT']:

			affected_haplotype = 'ALT'

		 # if the affected haplotype is the REF then get the ref AD for the plasma sample
		if affected_haplotype == 'REF':

			return df[plasma_sample+'_AD_alt']

		elif affected_haplotype == 'ALT':

			return df[plasma_sample+'_AD_ref'] 
		
		
		
	def get_hapA_counts_type_4b(df, maternal_sample, paternal_sample, proband_sample, plasma_sample):
		"""
		Count the number of reads in the plasma sample which are the same as the allele that the \
		affected proband got from mum.

		i.e. the allele which the paternal sample is hom for.

		"""
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# allele which is in mum and not dad
		hapa_allele = set(mat_gts) - set(pat_gts)

		hapa_allele = list(hapa_allele)[0]

		#is that allele the ref or alt?
		affected_haplotype = ''

		if hapa_allele == df['REF']:

			affected_haplotype = 'REF'

		elif hapa_allele == df['ALT']:

			affected_haplotype = 'ALT'

		 # if the affected haplotype is the REF then get the ref AD for the plasma sample
		if affected_haplotype == 'REF':

			return df[plasma_sample+'_AD_ref']

		elif affected_haplotype == 'ALT':

			return df[plasma_sample+'_AD_alt']     
		
		
	def get_hapB_counts_type_4b(df, maternal_sample, paternal_sample, proband_sample, plasma_sample ):
		"""
		Count the number of reads in the plasma sample which are the same as the allele that the \
		affected proband got from mum.

		i.e. the allele which the paternal sample is hom for.

		"""
		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# allele which is in mum and not dad
		hapa_allele = set(mat_gts) - set(pat_gts)

		hapa_allele = list(hapa_allele)[0]

		#is that allele the ref or alt?
		affected_haplotype = ''

		if hapa_allele == df['REF']:

			affected_haplotype = 'REF'

		elif hapa_allele == df['ALT']:

			affected_haplotype = 'ALT'

		 # if the affected haplotype is the REF then get the ref AD for the plasma sample
		if affected_haplotype == 'REF':

			return df[plasma_sample+'_AD_alt']

		elif affected_haplotype == 'ALT':

			return df[plasma_sample+'_AD_ref']         
		
		
			
	def sprt_calculate_upper_boundry(total_read_count, d_value, g_value):
		"""
		Calculate the upper value for the SPRT classification.

		"""

		return (( np.log(1200) / total_read_count ) - np.log(d_value)) / np.log(g_value)
	
	def sprt_calculate_lower_boundry(total_read_count, d_value, g_value):
		"""
		Calculate the lower value for the SPRT classification.

		"""

		return (( np.log(float(1) / float(1200) )/ total_read_count ) - np.log(d_value)) / np.log(g_value)
	
	
	def sprt_perform_test(df, d_value, g_value, min_snps_per_block):
		"""
		Function to call the status of each snp as hapA, hapB or Unclassified.

		Input:

		df = The processed dataframe containing the X-linked SNP data.
		d_value = The value calculated by sprt_calculate_d()
		g_value = The value calculated by sprt_calculate_g()

		Output:

		df = A new dataframe containing seven additional rows:

		cumulative_sum_hapa
		cumulative_sum_hapb
		cumulative_ratio
		cumulative_total
		upper_boundry
		lower_boundry
		status


		"""

		#df = df.copy() 

		cumulative_sum_hapa = 0
		cumulative_sum_hapb = 0
		snp_count = 0

		# Loop through each snp
		for snp in df.itertuples():

			# Calculate the cumlative sum
			cumulative_sum_hapa = int(cumulative_sum_hapa) + int(snp.hapA_count)
			cumulative_sum_hapb = int(cumulative_sum_hapb) + int(snp.hapB_count)

			# Calculate the cumlative allele balance between hapA and hapB
			cumulative_ratio = cumulative_sum_hapa / (cumulative_sum_hapa + cumulative_sum_hapb)

			cumulative_total = cumulative_sum_hapa + cumulative_sum_hapb

			# Calculate upper and lower SPRT boundaries
			upper_boundry = NIPDTools.sprt_calculate_upper_boundry(cumulative_total, d_value, g_value )

			lower_boundry = NIPDTools.sprt_calculate_lower_boundry(cumulative_total, d_value, g_value)

			# call blocks either hapA, hapB or Unclassified

			if (cumulative_ratio > upper_boundry) and snp_count >= min_snps_per_block:

				status = 'hapA'
				cumulative_sum_hapa = 0
				cumulative_sum_hapb = 0
				snp_count = 0


			elif cumulative_ratio < lower_boundry and snp_count >= min_snps_per_block:

				status = 'hapB'
				cumulative_sum_hapa = 0
				cumulative_sum_hapb = 0
				snp_count = 0

			else:

				status = 'Unclassified'

			snp_count = snp_count + 1

			#update dataframe with calculated values
			df.at[snp.Index, 'cumulative_sum_hapa'] = cumulative_sum_hapa
			df.at[snp.Index, 'cumulative_sum_hapb'] = cumulative_sum_hapb
			df.at[snp.Index, 'cumulative_ratio'] = cumulative_ratio
			df.at[snp.Index, 'cumulative_total'] = cumulative_total
			df.at[snp.Index, 'upper_boundry'] = upper_boundry
			df.at[snp.Index, 'lower_boundry'] = lower_boundry
			df.at[snp.Index, 'status'] = status

		return df
	
	
	def call_haplotype_blocks(df):
		"""
		Call the haplotype blocks along the chromosome.
		
		Takes a dataframe with a status column - each row should be either \
		
		HapA, HapB, or Unclassified.
		
		"""
	
		block_df = []

		prev_status = 'Unclassified'
		snp_count = 0
		block_id = 0
		first = True

		for snp in df.itertuples():

			if first == True:

				prev_pos = snp.POS
				first = False

			status = snp.status

			if status != prev_status:

				block_df.append([block_id, prev_pos, snp.POS, snp_count, status, snp.CHROM, snp.GENE])

				prev_pos = snp.POS
				snp_count = 0
				block_id = block_id + 1


			snp_count = snp_count + 1

		return block_df
	
	def get_hapb_pat_count(df, maternal_sample, paternal_sample, proband_sample, plasma_sample):
		"""
		We define the hapb pat count for SNPS where:

		a) mum is hom and dad is het
		b) the affected proband is also hom (for same allele as mum)

		For these SNPS we get the plasma count for the allele which is unique to dad.

		This would be the allele in dad which is NOT in the affected proband. 

		"""

		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		# If the affected proband is hom for same allele as mum
		if proband_gts == mat_gts:

			# get count of normal pat allelle which is the one not in the mum

			not_in_mum =  set(pat_gts) - set(mat_gts)
			not_in_mum = list(not_in_mum)[0]


			normal_haplotype = ''

			if not_in_mum == df['REF']:

				normal_haplotype = 'REF'

			elif not_in_mum == df['ALT']:

				normal_haplotype = 'ALT'

			# if the affected haplotype is the REF then get the ref AD for the plasma sample
			if normal_haplotype == 'REF':

				return df[plasma_sample+'_AD_ref']

			elif normal_haplotype == 'ALT':

				return df[plasma_sample+'_AD_alt']

		else:

			return None
		
	def get_hapa_pat_count(df, maternal_sample, paternal_sample, proband_sample, plasma_sample):

		"""
		We define the hapa pat count for SNPS where:

		a) mum is hom and dad is het
		b) the affected proband is het (for same allele as mum)

		For these SNPS we get the plasma count for the allele which is unique to dad.

		This would be the allele in dad which IS in the affected proband. 

		"""

		mat_gts = NIPDTools.get_genotypes_as_list(df[maternal_sample + '_GT'])
		pat_gts = NIPDTools.get_genotypes_as_list(df[paternal_sample + '_GT'])
		proband_gts = NIPDTools.get_genotypes_as_list(df[proband_sample + '_GT'])
		plasma_gts = NIPDTools.get_genotypes_as_list(df[plasma_sample + '_GT'])

		if proband_gts != mat_gts:

			# get count of normal pat allelle which is the one not in the mum

			not_in_mum =  set(pat_gts) - set(mat_gts)
			not_in_mum = list(not_in_mum)[0]

			# get the count of not_in_mum allele

			normal_haplotype = ''

			if not_in_mum == df['REF']:

				normal_haplotype = 'REF'

			elif not_in_mum == df['ALT']:

				normal_haplotype = 'ALT'

			  # if the affected haplotype is the REF then get the ref AD for the plasma sample
			if normal_haplotype == 'REF':

				return df[plasma_sample+'_AD_ref']

			elif normal_haplotype == 'ALT':

				return df[plasma_sample+'_AD_alt']


		else:

			return None
		
		
	def perform_ks_test(df, min_p_value, min_snps_per_block):
		"""
		Perfrom the Kolmogorov-Smirnov (KS) test on the counts from:

		a) the normal pat haplotype (HapB)
		b) the mutant pat haplotype (HapA)

		"""

		normal_pat_array = []
		mutant_pat_array = []

		snp_count = 1

		for row in df.itertuples():

			if row.HapB_count != None and pd.isna(row.HapB_count) == False :
				
				normal_pat_array.append(int(row.HapB_count))

			if row.HapA_count != None and pd.isna(row.HapA_count) == False:

				mutant_pat_array.append(int(row.HapA_count))

			result = stats.ks_2samp(normal_pat_array, mutant_pat_array )

			if result.pvalue < min_p_value and snp_count > min_snps_per_block:

				if statistics.median(normal_pat_array) >  statistics.median(mutant_pat_array):

					status = 'hapB'

				elif  statistics.median(mutant_pat_array) >  statistics.median(normal_pat_array):

					status = 'hapA'

				else:

					status = 'Unclassified'

				normal_pat_array = []
				mutant_pat_array = []
				snp_count = 1

			else:

				status = 'Unclassified'

				snp_count = snp_count + 1

			df.at[row.Index, 'status'] = status
			df.at[row.Index, 'pvalue'] = result.pvalue

		return df


	def get_hapA_counts_xlinked(df, proband_sample, plasma_sample):
	    """
	    count the number of reads in the plasma sample which match the affected haplotype
	    
	    e.g. how many of the reads in the plasma match the genotype of the affected son?
	    
	    """
	    
	    # Get the genotype of the proband 
	    gt_proband = NIPDTools.get_genotypes_as_list(df[proband_sample+'_GT'])

	    # proband should be hom 
	    assert gt_proband[0] == gt_proband[1]

	    gt_proband = gt_proband[0]
	    
	    # is the affected haplotype ref or alt?
	    
	    affected_haplotype = ''
	    
	    if gt_proband == df['REF']:
	        
	        affected_haplotype = 'REF'
	        
	    elif gt_proband == df['ALT']:
	        
	        affected_haplotype = 'ALT'
	        
	    # if the affected haplotype is the REF then get the ref AD for the plasma sample
	        
	    if affected_haplotype == 'REF':
	        
	        return df[plasma_sample+'_AD_ref']
	    
	    elif affected_haplotype == 'ALT':
	        
	        return df[plasma_sample+'_AD_alt']



	def get_hapB_counts_xlinked(df, proband_sample, plasma_sample):
	    """
	    count the number of reads in the plasma sample which don't match the affected haplotype
	    
	    """
	                         
	    # Get the genotype of the proband 
	    gt_proband = NIPDTools.get_genotypes_as_list(df[proband_sample+'_GT'])
	    
	    # proband should be hom 
	    assert gt_proband[0] == gt_proband[1]

	    gt_proband = gt_proband[0]
	    
	    # is the affected haplotype ref or alt?
	    
	    affected_haplotype = ''
	    
	    if gt_proband == df['REF']:
	        
	        affected_haplotype = 'REF'
	        
	    elif gt_proband == df['ALT']:
	        
	        affected_haplotype = 'ALT'
	        
	    # if the affected haplotype is the REF then get the ref AD for the plasma sample
	        
	    if affected_haplotype == 'REF':
	        
	        return df[plasma_sample+'_AD_alt']
	    
	    elif affected_haplotype == 'ALT':
	        
	        return df[plasma_sample+'_AD_ref']
