"""
Pipeline for single gene Non Invasive Prenatal Diagnosis

See readme.md for details

"""
from pathlib import Path 


#-----------------------------------------------------------------------------------------------------------------#
# Configuration variables
#-----------------------------------------------------------------------------------------------------------------#

# Which YAMl config to use
config_location = "development_local.yaml"

configfile: config_location

# Worksheet ID
worksheet = config["seqID"]

# Get all family IDs from config file
families = config["families"].keys()

# How many lanes do we have?
folder, name, lanes = glob_wildcards("{folder}/{sample_number}_{lane}_R1_001.fastq.gz")
lanes =  list((set(lanes)))

# What are our samples and sample numbers?
folder, sample_names, sample_numbers = glob_wildcards("{folder}/{sample_name}_{sample_number}_L001_R1_001.fastq.gz")

# Get other data from config file
chromosomes = config["chromosomes"]
panel = config["panel"]

#-----------------------------------------------------------------------------------------------------------------#
# Utility Functions For Getting Files
#-----------------------------------------------------------------------------------------------------------------#


def get_fastqc(wildcards):
	"""	
	Function to return the fastqc file input into multiqc
	https://groups.google.com/forum/#!topic/snakemake/xxdADOSK7mY

	"""
	file_list = []
	for lane in lanes:

		for sample_name, sample_number in zip(sample_names,sample_numbers ):

			file_list.append("output/qc_reports/fastqc/" + sample_name + "_" + sample_number + "_" + lane + "_R1_001.qfilter_fastqc.zip" )
			file_list.append("output/qc_reports/fastqc/" + sample_name + "_" + sample_number + "_" + lane + "_R2_001.qfilter_fastqc.zip" )

	return file_list

#-----------------------------------------------------------------------------------------------------------------#
# Main pipeline
#-----------------------------------------------------------------------------------------------------------------#

# All function pulls all rules together
rule all:
	input:
		expand("output/family_csvs/{worksheet}_all_chr_qfiltered_anno_selected_{FAMID}.csv", worksheet = worksheet, FAMID = families),
		expand("output/vep_family_csvs/{worksheet}_all_chr_qfiltered_anno_norm_vep_{FAMID}.csv", worksheet = worksheet, FAMID = families),
		expand("output/qc_reports/multiqc/{worksheet}.html", worksheet = worksheet),
		expand("output/qfiltered_jointvcf_anno_norm_vep/{worksheet}_all_chr_qfiltered_anno_norm_vep.vcf", worksheet=worksheet)


#-----------------------------------------------------------------------------------------------------------------#
# Preprocessing and Read Level Quality Control
#-----------------------------------------------------------------------------------------------------------------#

# Run the fastp program to generate a read quality report and trim reads
rule fastp:
	input:
		fwd = "{sample_name}/{sample_name}_{sample_number}_{lane}_R1_001.fastq.gz",
		rev = "{sample_name}/{sample_name}_{sample_number}_{lane}_R2_001.fastq.gz"
	output:
		html = "output/qc_reports/fastp/{sample_name}_{sample_number}_{lane}_fastp.html",
		json = "output/qc_reports/fastp/{sample_name}_{sample_number}_{lane}_fastp.json",
		fwd = temp("output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz"),
		rev = temp("output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz")
	threads:
		config["fastp_threads"]
	shell:
		"fastp -i {input.fwd} "
		"-I {input.rev} "
		"-o {output.fwd} "
		"-O {output.rev} "
		"-h {output.html} "
		"-j {output.json} "
		"--length_required 35 "
		"--disable_quality_filtering "
		"--detect_adapter_for_pe "
		"-w {threads}"


# Run fastqc
rule fastqc:
	input:
		fwd = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz",
		rev = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz"
	output:
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R1_001.qfilter_fastqc.html",
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R1_001.qfilter_fastqc.zip",	
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R2_001.qfilter_fastqc.html",
		"output/qc_reports/fastqc/{sample_name}_{sample_number}_{lane}_R2_001.qfilter_fastqc.zip"
	threads:
		config["fastqc_threads"]
	params:
		temp_dir = config["fastqc_temp_dir"]
	shell:
		"fastqc "
		"--threads {threads} "
		"--dir {params.temp_dir} "
		"--outdir output/qc_reports/fastqc "
		"{input.fwd} "
		"{input.rev}"

#-----------------------------------------------------------------------------------------------------------------#
# Alignment
#-----------------------------------------------------------------------------------------------------------------#

# Align reads with bwa mem
rule bwa_align:
	input:
		fwd = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R1_001.qfilter.fastq.gz",
		rev = "output/qfiltered_reads/{sample_name}_{sample_number}_{lane}_R2_001.qfilter.fastq.gz"
	output:
		temp("output/alignments/{sample_name}_{sample_number}_{lane}.bam")
	threads:
		config["bwa_threads"]
	params:
		ref = config["bwa_reference"],
		worksheet = config["seqID"],
		centre = config["centre"],
		samtools_temp_dir = config["samtools_temp_dir"]
	shell:
		"bwa mem "
		"-t {threads} "
		"-M "
		"-R '@RG\\tID:{params.worksheet}.{wildcards.lane}\\tCN:{params.centre}\\tSM:{wildcards.sample_name}\\tLB:{params.worksheet}\\tPL:ILLUMINA' "
		"{params.ref} {input.fwd} {input.rev} | "
		"samtools view -Sb - | "
		"samtools sort -T {params.samtools_temp_dir} -O bam > {output}"

# Index bam file
rule index_original_bam:
	input:
		"output/alignments/{sample_name}_{sample_number}_{lane}.bam"
	output:
		temp("output/alignments/{sample_name}_{sample_number}_{lane}.bam.bai")
	shell:
		"samtools index {input}"

# Merge the bams and mark duplicates
rule merge_and_remove_duplicates:
	input:
		bams = expand("output/alignments/{{sample_name}}_{{sample_number}}_{lane}.bam", lane=lanes),
		bam_indexes = expand("output/alignments/{{sample_name}}_{{sample_number}}_{lane}.bam.bai", lane=lanes),
	output:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		metrics = "output/qc_reports/mark_duplicates/{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt"
	params:
		temp = config["picard_temp_dir"],
		merge_duplicates_max_records = config["merge_duplicates_max_records"],
		files = lambda wildcards, input: " I=".join(input.bams),
		java_options = config["picard_memory_options"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard {params.java_options} MarkDuplicates I={params.files} "
		"O={output.bam} "
		"METRICS_FILE={output.metrics} "
		"CREATE_INDEX=true "
		"MAX_RECORDS_IN_RAM={params.merge_duplicates_max_records} "
		"VALIDATION_STRINGENCY=SILENT "
		"TMP_DIR={params.temp} "


#-----------------------------------------------------------------------------------------------------------------#
# Post Alignment QC
#-----------------------------------------------------------------------------------------------------------------#

# Create an interval file from the BED file for use in Picard tools such as CollectHsMetrics
rule create_interval_file:
	input:
		bed = ancient(config["roi_bed_file"]),
	output:
		bed = temp("output/config/" + Path(config["roi_bed_file"]).name.split(".")[0] + ".interval_list"),
	params:
		sequence_dict = config["reference_sequence_dict"],
		java_home = config["java_home"]
	shell:
		"""
		export JAVA_HOME={params.java_home}
		picard BedToIntervalList I={input.bed} O={output.bed} SD={params.sequence_dict}
		"""

# Collect some insert size metrics using Picard
rule collect_insert_size_metrics:
	input:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		txt="output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.txt",
		pdf="output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.pdf"
	params:
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectInsertSizeMetrics I={input.bam} O={output.txt} HISTOGRAM_FILE={output.pdf}"


# Collect some HS metrics using picard
rule collect_hs_metrics:
	input:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		intervals_capture = "output/config/" + Path(config["roi_bed_file"]).name.split(".")[0] + ".interval_list",
	output:
		"output/qc_reports/hs_metrics/{sample_name}_{sample_number}_HsMetrics.txt"
	params:
		ref = config["reference"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectHsMetrics I={input.bam} O={output} R={params.ref} "
		"BAIT_INTERVALS={input.intervals_capture} TARGET_INTERVALS={input.intervals_capture}"

# Collect alignment summary metrics using picard
rule collect_alignment_metrics:
	input:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		"output/qc_reports/alignment_metrics/{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt"
	params:
		ref = config["reference"],
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard CollectAlignmentSummaryMetrics I={input.bam} O={output} R={params.ref}"

# Calculate per base coverage with sambamba
rule get_per_base_coverage:
	input:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		bed = "output/config/sorted_beds/{panel}_sorted.bed".format(panel=panel)
	output:
		"output/qc_reports/depth/{sample_name}_{sample_number}_per_base.coverage"
	shell:
		"sambamba depth base "
		"-L {input.bed} --min-coverage=0 "
		"{input.bam} > {output} "

# Calculate exon coverage with sambamba
rule get_exon_coverage:
	input:
		bam = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		"output/qc_reports/depth/{sample_name}_{sample_number}_exon.coverage"
	params:
		bed = config["exon_bed_file"]
	shell:
		"sambamba depth region "
		"-L {params.bed} -T 30 -T 100 -T 150 -T 200 -T 500 -T 750 -T 1000 "
		"{input.bam} > {output} "

# Get the coverage over each of the target snps
rule get_snp_coverage:
	input:
		"output/qc_reports/depth/{sample_name}_{sample_number}_per_base.coverage"
	output:
		snp_coverage = "output/qc_reports/depth/{sample_name}_{sample_number}_snps.coverage",
		per_base_bed = temp("output/qc_reports/depth/{sample_name}_{sample_number}_per_base.bed")
	params:
		snp_bed = config['snp_bed']
	shell:
		"""
		tail {input} -n+2 | awk 'BEGIN {{ OFS =\"\\t\" }} {{print $1,$2,$2+1,$3}}' > {output.per_base_bed}
		bedtools intersect -a {params.snp_bed} -b {output.per_base_bed} -loj > {output.snp_coverage}
		"""

# Relatedness Testing
rule relatedness_test:
	input:
		"output/qfiltered_jointvcf_anno/{worksheet}_all_chr_qfiltered_anno.vcf"
	output:
		"output/qc_reports/relatedness/{worksheet}.relatedness2"
	shell:
		"vcftools --relatedness2 "
		"--out output/qc_reports/relatedness/{wildcards.worksheet} "
		"--vcf {input} "
	
# Multiqc to compile all qc data into one file
rule multiqc:
	input:
		insert_size_metrics =  expand("output/qc_reports/insert_size_metrics/{sample_name}_{sample_number}_InsertSizeMetrics.txt", zip, sample_name=sample_names, sample_number=sample_numbers),
		hs_metrics_metrics = expand("output/qc_reports/hs_metrics/{sample_name}_{sample_number}_HsMetrics.txt", zip, sample_name=sample_names, sample_number=sample_numbers),
		alignment_metrics = expand("output/qc_reports/alignment_metrics/{sample_name}_{sample_number}_AlignmentSummaryMetrics.txt", zip, sample_name=sample_names, sample_number=sample_numbers),
		mark_duplicate_metrics = expand("output/qc_reports/mark_duplicates/{sample_name}_{sample_number}_MarkDuplicatesMetrics.txt", zip, sample_name=sample_names, sample_number=sample_numbers),
		base_depth = expand("output/qc_reports/depth/{sample_name}_{sample_number}_per_base.coverage", zip, sample_name=sample_names, sample_number=sample_numbers ),
		exon_depth = expand("output/qc_reports/depth/{sample_name}_{sample_number}_exon.coverage", zip, sample_name=sample_names, sample_number=sample_numbers ),
		snp_depth = expand("output/qc_reports/depth/{sample_name}_{sample_number}_snps.coverage", zip, sample_name=sample_names, sample_number=sample_numbers ),
		fastqc = get_fastqc,
		relatedness = expand("output/qc_reports/relatedness/{worksheet}.relatedness2", worksheet=worksheet)
	output:
		html = "output/qc_reports/multiqc/" + worksheet + ".html",
		data = directory("output/qc_reports/multiqc/" + worksheet + "_data")
	params:
		worksheet = worksheet
	shell:
		"multiqc --filename {params.worksheet} --exclude fastp --outdir output/qc_reports/multiqc/ output/qc_reports"

#-----------------------------------------------------------------------------------------------------------------#
# SNP and Small Indel Calling with GATK Haplotype Caller
#-----------------------------------------------------------------------------------------------------------------#

# Sort ROI bed for splitting by bedextract
rule sort_capture_bed:
	input:
		ancient(config["padded_roi_bed_file"])
	output:
		temp("output/config/sorted_beds/{{panel}}_sorted.bed".format(panel=panel))
	shell:
		"sort-bed {input} > {output}"

# Split the bed by chromosome for input into create_gvcfs
rule split_bed_by_chromosome:
	input:
		"output/config/sorted_beds/{panel}_sorted.bed".format(panel=panel)
	output:
		expand("output/config/split_capture_bed/{chr}.bed", chr=chromosomes)
	params:
		chromosomes = chromosomes
	shell:
		"for chr in {params.chromosomes}; do bedextract $chr {input} > output/config/split_capture_bed/$chr.bed; done"

# Create GVCF using Haplotype Caller for each sample chromosome combination
rule create_gvcfs:
	input:
		bam_file = "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index= "output/merged_bams/{sample_name}_{sample_number}_merged_nodups.bai",
		bed = "output/config/split_capture_bed/{chr}.bed"
	output:
		gvcf = temp("output/gvcfs/{sample_name}_{sample_number}_chr{chr}.g.vcf"),
		gvcf_index = temp("output/gvcfs/{sample_name}_{sample_number}_chr{chr}.g.vcf.idx")
	params:
		ref = config["reference"],
		padding = config['interval_padding_haplotype_caller'],
		java_options = config['gatk_hc_java_options']
	shell:
		"gatk --java-options '{params.java_options}' HaplotypeCaller -R {params.ref} "
		"-I {input.bam_file} "
		"--emit-ref-confidence GVCF "
		"-O {output.gvcf} "
		"-L {input.bed} "
		"--interval-padding {params.padding}"


# Consolidate all samples into a genomics db for joint genotyping
rule create_genomics_db:
	input:
		gvcfs = expand("output/gvcfs/{sample_name}_{sample_number}_chr{{chr}}.g.vcf" , zip, sample_name=sample_names, sample_number=sample_numbers),
		gvcf_indexes = expand("output/gvcfs/{sample_name}_{sample_number}_chr{{chr}}.g.vcf.idx" , zip, sample_name=sample_names, sample_number=sample_numbers),
		bed = "output/config/split_capture_bed/{chr}.bed"
	output:
		temp(directory("output/genomicdbs/{worksheet}_chr{chr}"))
	params:
		files = lambda wildcards, input: " -V ".join(input.gvcfs),
		java_options = config["gatk_genomics_db_java_options"],
		padding = config['interval_padding_haplotype_caller']
	shell:
		"gatk --java-options '{params.java_options}' "
		" GenomicsDBImport -V {params.files} "
		"--genomicsdb-workspace-path {output} "
		"-L {input.bed} "
		"--interval-padding {params.padding}"


# Genotype the gvcfs and produce a joint vcf
rule genotype_gvcfs:
	input:
		db = directory("output/genomicdbs/{worksheet}_chr{chr}"),
		bed = "output/config/split_capture_bed/{chr}.bed"
	output:
		vcf = temp("output/jointvcf_per_chr/{worksheet}_chr{chr}.vcf"),
		index = temp("output/jointvcf_per_chr/{worksheet}_chr{chr}.vcf.idx"),
	params:
		ref = config["reference"],
		java_options = config['gatk_hc_java_options'],
		padding = config['interval_padding_haplotype_caller']
	shell:
		"gatk --java-options '{params.java_options}'  GenotypeGVCFs -R {params.ref} "
		"-V gendb://{input.db} "
		"-G StandardAnnotation "
		"-O {output.vcf} "
		"-L {input.bed} "
		"--interval-padding {params.padding} "

# Combine the chromsome vcfs into one final vcf with all samples and all chromosomes
rule collect_vcfs:
	input:
		vcf = expand("output/jointvcf_per_chr/{{worksheet}}_chr{chr}.vcf", chr= chromosomes),
		index = expand("output/jointvcf_per_chr/{{worksheet}}_chr{chr}.vcf.idx", chr= chromosomes),
	output:
		vcf = temp("output/jointvcf/{worksheet}_all_chr.vcf"),
		index = temp("output/jointvcf/{worksheet}_all_chr.vcf.idx")
	params:
		files = lambda wildcards, input: " I=".join(input.vcf),
		java_home = config["java_home"]
	shell:
		"export JAVA_HOME={params.java_home}; picard GatherVcfs "
		"I={params.files} "
		"O={output.vcf}"

#-----------------------------------------------------------------------------------------------------------------#
# SNP and Small Indel Calling for de novo calling - TODO - Platypus or Vardict caller? Depth and PL information? recalling?
#-----------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------------------#
# Filter and Annotate Variants
#-----------------------------------------------------------------------------------------------------------------#

# Use hard filtering on quality attributes
# See https://gatkforums.broadinstitute.org/gatk/discussion/6925
# for more information on the values chosen here.
rule hard_filter_vcf_gatk:
	input:
		vcf = "output/jointvcf/{worksheet}_all_chr.vcf",
		index = "output/jointvcf/{worksheet}_all_chr.vcf.idx",
	output:
		vcf = temp("output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf"),
		index = temp("output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf.idx")
	params:
		ref = config["reference"],
		min_QD = config['min_QD'],
		max_FS = config['max_FS'],
		max_SOR = config['max_SOR'],
		min_MQ = config['min_MQ'],
		min_MQRankSum = config['min_MQRankSum'],
		min_ReadPosRankSum = config['min_ReadPosRankSum']
	shell:
		"gatk VariantFiltration -R {params.ref} -O {output.vcf} "
		"--variant {input.vcf} "
		"--filter-expression 'QD < {params.min_QD}' --filter-name 'LOW_QD' "
		"--filter-expression 'FS > {params.max_FS}' --filter-name 'HIGH_FS' "
		"--filter-expression 'SOR > {params.max_SOR}' --filter-name 'HIGH_SOR' "
		"--filter-expression 'MQ < {params.min_MQ}' --filter-name 'LOW_MQ' "
		"--filter-expression 'MQRankSum < {params.min_MQRankSum}' --filter-name 'LOW_MQRankSum' "
		"--filter-expression 'ReadPosRankSum < {params.min_ReadPosRankSum}' --filter-name 'LOW_ReadPosRankSum' "

# Compress and index vcf so we can annotate with gene name
rule compress_and_index_vcf_gatk:
	input:
		"output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf"
	output:
		vcf = temp("output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf.gz"),
		index = temp("output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf.gz.tbi")
	shell:
		"bgzip {input} && tabix {output.vcf}"


# Compress and index bedfile so we can use bcftools to annotate with gene
rule compress_and_index_bed_file_gatk:
	input:
		ancient(config["gene_bed_file"])
	output:
		bed = "output/config/" + Path(config["gene_bed_file"]).name.split(".")[0] + ".bed.gz",
		index = "output/config/" + Path(config["gene_bed_file"]).name.split(".")[0] + ".bed.gz.tbi"
	shell:
		"bgzip --stdout {input} > {output.bed} && tabix {output.bed}"


# Annotate with gene bed file.
rule annotate_vcf_with_gene_gatk:
	input:
		vcf = "output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf.gz",
		vcf_index = "output/qfiltered_jointvcf/{worksheet}_all_chr_qfiltered.vcf.gz",
		bed = "output/config/" + Path(config["gene_bed_file"]).name.split(".")[0] + ".bed.gz",
		bed_index = "output/config/" + Path(config["gene_bed_file"]).name.split(".")[0] + ".bed.gz.tbi"
	output:
		vcf = temp("output/qfiltered_jointvcf_anno/{worksheet}_all_chr_qfiltered_anno.vcf")
	shell:
		"bcftools annotate -a {input.bed} "
		"-c CHROM,FROM,TO,GENE "
		"-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">') "
		"{input.vcf} -o {output.vcf} "


# Use vt to split multiallelics and normalise variants
rule decompose_and_normalise_gatk:
	input:
		"output/qfiltered_jointvcf_anno/{worksheet}_all_chr_qfiltered_anno.vcf"
	output:
		temp("output/qfiltered_jointvcf_anno_norm/{worksheet}_all_chr_qfiltered_anno_norm.vcf")
	params:
		ref = config["reference"]
	shell:
		"cat {input} | "
		"vt decompose -s - | "
		"vt normalize -r {params.ref} - > {output}"

# Annotate the GATK vcf using VEP
rule annotate_vep_gatk:
	input:
		"output/qfiltered_jointvcf_anno_norm/{worksheet}_all_chr_qfiltered_anno_norm.vcf"
	output:
		vcf = "output/qfiltered_jointvcf_anno_norm_vep/{worksheet}_all_chr_qfiltered_anno_norm_vep.vcf",
		summary = temp("output/qfiltered_jointvcf_anno_norm_vep/{worksheet}_all_chr_qfiltered_anno_norm_vep.vcf_summary.html"),
		warnings = temp("output/qfiltered_jointvcf_anno_norm_vep/{worksheet}_all_chr_qfiltered_anno_norm_vep.vcf_warnings.txt")
	params:
		vep_cache = config["vep_cache_location"],
		ref = config["reference"],
		gnomad_genomes = config["gnomad_genomes"],
		gnomad_exomes = config["gnomad_exomes"],
	threads:
		config["vep_threads"]
	shell:
		"vep --verbose "
		"--format vcf "
		"--everything "
		"--fork {threads} "
		"--species homo_sapiens "
		"--assembly GRCh37  "
		"--input_file {input}  "
		"--output_file {output.vcf} "
		"--force_overwrite "
		"--cache "
		"--dir  {params.vep_cache} "
		"--fasta {params.ref} "
		"--offline "
		"--cache_version 94 "
		"--no_escape "
		"--shift_hgvs 1 "
		"--vcf "
		"--refseq "
		"--flag_pick "
		"--pick_order biotype,canonical,appris,tsl,ccds,rank,length "
		"--exclude_predicted "
		"--custom {params.gnomad_genomes},gnomADg,vcf,exact,0,AF_POPMAX "
		"--custom {params.gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX "

# Create VEP Family CSVs
rule convert_gatk_vep_vcf_to_csv_per_family:
	input:
		"output/qfiltered_jointvcf_anno_norm_vep/{worksheet}_all_chr_qfiltered_anno_norm_vep.vcf"
	output:
		expand("output/vep_family_csvs/{{worksheet}}_all_chr_qfiltered_anno_norm_vep_{FAMID}.csv", FAMID=config['families'].keys())
	params:
		config = config_location
	shell:
		"python utils/pipeline_scripts/make_gatk_vep_csv.py "
		"--input {input} "
		"--config {params.config} "
		"--output_dir output/vep_family_csvs/ "
		"--output_prefix {wildcards.worksheet}_all_chr_qfiltered_anno_norm_vep"


# Select only Biallelic SNPs which pass filtering i.e. exclude indels and multialleleics and fails
rule select_relevant_variants_gatk:
	input:
		vcf = "output/qfiltered_jointvcf_anno/{worksheet}_all_chr_qfiltered_anno.vcf",
	output:
		vcf = temp("output/qfiltered_jointvcf_anno_selected/{worksheet}_all_chr_qfiltered_anno_selected.vcf"),
		index = temp("output/qfiltered_jointvcf_anno_selected/{worksheet}_all_chr_qfiltered_anno_selected.vcf.idx")
	params:
		ref = config["reference"]
	shell:
		"gatk SelectVariants "
		"-R {params.ref} "
		"-V {input} "
		"-O {output.vcf} "
		"-select-type SNP "
		"--restrict-alleles-to BIALLELIC "
		"--exclude-filtered true "

# Split each vcf by family as specified in the config file
rule split_vcf_by_family_gatk:
	input:
		"output/qfiltered_jointvcf_anno_selected/{worksheet}_all_chr_qfiltered_anno_selected.vcf"
	output:
		temp(expand("output/family_vcfs/{{worksheet}}_all_chr_qfiltered_anno_selected_{FAMID}.vcf", FAMID=config['families'].keys())),
		temp(expand("output/family_vcfs/{{worksheet}}_all_chr_qfiltered_anno_selected_{FAMID}.vcf.idx", FAMID=config['families'].keys()))
	params:
		ref = config["reference"],
		config = config_location
	shell:
		"python utils/pipeline_scripts/split_joint_vcf_by_family.py "
		"--input {input} "
		"--config {params.config} "
		"--ref {params.ref} "
		"--output_dir output/family_vcfs "
		"--output_prefix {wildcards.worksheet}_all_chr_qfiltered_anno_selected"

# Create family CSVs from VCFs ready for SPRT analysis
rule create_family_csv_gatk:
	input:
		"output/family_vcfs/{worksheet}_all_chr_qfiltered_anno_selected_{FAMID}.vcf"
	output:
		"output/family_csvs/{worksheet}_all_chr_qfiltered_anno_selected_{FAMID}.csv"
	shell:
		"gatk VariantsToTable "
		"-V {input} "
		"-O {output} "
		"-F CHROM "
		"-F POS "
		"-F REF "
		"-F ALT "
		"-F ID "
		"-F GENE "
		"-GF GT "
		"-GF AD "
		"-GF DP "
		"-GF GQ"


