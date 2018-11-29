"""
Pipeline for single gene Non Invasive Prenatal Diagnosis

"""
from pathlib import Path 

#-----------------------------------------------------------------------------------------------------------------#
# Configuration variables
#-----------------------------------------------------------------------------------------------------------------#

config_location = "config/development.yaml"
configfile: config_location
initial_files = ["HG002_son_DMD","HG004_mum_DMD","HG004_plasma_DMD"]
sample_numbers =   ["S1","S2","S3"]
worksheet = 'test'

#-----------------------------------------------------------------------------------------------------------------#
# Main pipeline
#-----------------------------------------------------------------------------------------------------------------#



# Run the fastp program to generate a read quality report and trim reads
rule fastp_lane1:
	input:
		fwd = "input/{sample_name}/{sample_name}_{sample_number}_L001_R1_001.fastq",
		rev = "input/{sample_name}/{sample_name}_{sample_number}_L001_R2_001.fastq"
	output:
		html = "output/{sample_name}/qc_reports/fastp/{sample_name}_{sample_number}_fastp.html",
		json = "output/{sample_name}/qc_reports/fastp/{sample_name}_{sample_number}_fastp.json",
		fwd = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R1_001.qfilter.fastq",
		rev = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R2_001.qfilter.fastq"
	threads:
		config['fastp_threads']
	shell:
		"fastp -i {input.fwd} -I {input.rev} "
		"-o {output.fwd} -O {output.rev} "
		"-h {output.html} -j {output.json} "
		"-w {threads} "
		"-3"


# Check for inter species contamination
rule fastq_screen_lane1:
	input:
		fwd = "input/{sample_name}/{sample_name}_{sample_number}_L001_R1_001.fastq",
		rev = "input/{sample_name}/{sample_name}_{sample_number}_L001_R2_001.fastq"
	output:
		"output/{sample_name}/qc_reports/fastq_screen/{sample_name}_{sample_number}_L001_R1_001_screen.html",
		"output/{sample_name}/qc_reports/fastq_screen/{sample_name}_{sample_number}_L001_R1_001_screen.txt",
		"output/{sample_name}/qc_reports/fastq_screen/{sample_name}_{sample_number}_L001_R2_001_screen.html",
		"output/{sample_name}/qc_reports/fastq_screen/{sample_name}_{sample_number}_L001_R2_001_screen.txt"	
	threads:
		config['fastq_screen_threads']
	params:
		fastq_screen_config = config["fastq_screen_config"]
	shell:
		"fastq_screen --aligner bwa --threads {threads} "
		"--outdir output/{wildcards.sample_name}/qc_reports/fastq_screen/ "
		"--conf {params.fastq_screen_config} "
		"{input.fwd} {input.rev} "

# Run fastqc
rule fastqc_lane1:
	input:
		fwd = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R1_001.qfilter.fastq",
		rev = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R2_001.qfilter.fastq"
	output:
		"output/{sample_name}/qc_reports/fastqc/{sample_name}_{sample_number}_L001_R1_001.qfilter_fastqc.html",
		"output/{sample_name}/qc_reports/fastqc/{sample_name}_{sample_number}_L001_R1_001.qfilter_fastqc.zip",	
		"output/{sample_name}/qc_reports/fastqc/{sample_name}_{sample_number}_L001_R2_001.qfilter_fastqc.html",
		"output/{sample_name}/qc_reports/fastqc/{sample_name}_{sample_number}_L001_R2_001.qfilter_fastqc.zip"
	threads:
		config['fastqc_threads']
	params:
		temp_dir = config['fastqc_temp_dir']
	shell:
		"fastqc --threads {threads} --dir {params.temp_dir} "
		 "--outdir output/{wildcards.sample_name}/qc_reports/fastqc "
		 "{input.fwd} {input.rev} "


#-----------------------------------------------------------------------------------------------------------------#
# Alignment
#-----------------------------------------------------------------------------------------------------------------#

# Align reads with bwa mem
rule bwa_align_lane1:
	input:
		fwd = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R1_001.qfilter.fastq",
		rev = "output/{sample_name}/qfiltered_reads/{sample_name}_{sample_number}_L001_R2_001.qfilter.fastq"
	output:
		"output/{sample_name}/alignments/{sample_name}_{sample_number}_L001.bam"
	threads:
		config['bwa_threads']
	params:
		ref = config["bwa_reference"],
		worksheet = config["worksheet"],
		flowcell = config["flowcell"],
		centre = config["centre"],
	shell:
		"bwa mem -t {threads} -M -k 18 "
		"-R '@RG\\tID:{params.worksheet}.{params.flowcell}.L001\\tCN:{params.centre}\\tSM:{wildcards.sample_name}\\tLB:{params.worksheet}\\tPL:ILLUMINA' "
		"{params.ref} {input.fwd} {input.rev} | "
		"samtools view -Sb - | "
		"samtools sort -T {wildcards.sample_name}.temp -O bam > {output}"

# Index bam file
rule index_original_bam_lane1:
	input:
		"output/{sample_name}/alignments/{sample_name}_{sample_number}_L001.bam"
	output:
		"output/{sample_name}/alignments/{sample_name}_{sample_number}_L001.bam.bai"
	shell:
		"samtools index {input}"


# Mark duplicates using Picard
rule merge_and_mark_duplicates:
	input:
		lane1_bam = "output/{sample_name}/alignments/{sample_name}_{sample_number}_L001.bam",
		lane1_bam_index = "output/{sample_name}/alignments/{sample}_{sample_number}_L001.bam.bai"
	output:
		bam = "output/{sample_name}/merged/{sample}_{sample_number}_merged_nodups.bam",
		index = "output/{sample_name}/merged/{sample}_{sample_number}_merged_nodups.bai",
		metrics = "output/{sample_name}/qc_reports/{sample}_{sample_number}_MarkDuplicatesMetrics.txt"
	params:
		temp = config["picard_temp_dir"]
	shell:
		"picard MarkDuplicates I={input.lane1_bam} "
		"O={output.bam} "
        "METRICS_FILE={output.metrics} "
		"CREATE_INDEX=true "
        "MAX_RECORDS_IN_RAM=2000000 "
        "VALIDATION_STRINGENCY=SILENT "
        "TMP_DIR={params.temp} "

# Create an interval file from the BED file for use in Picard tools such as CollectHsMetrics
rule create_interval_file:
	input:
		config["capture_bed_file"]
	output:
		"output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list"
	params:
		sequence_dict = config["reference_sequence_dict"]
	shell:
		"picard BedToIntervalList I={input} O={output} SD={params.sequence_dict}" 


# Collect some insert size metrics using Picard
rule collect_insert_size_metrics:
	input:
		bam = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bam",
		index = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		txt="output/{sample_name}/qc_reports/insert_size_metrics/{sample_name}_{sample_number}.txt",
		pdf="output/{sample_name}/qc_reports/insert_size_metrics/{sample_name}_{sample_number}.pdf"	
	shell:
		"picard CollectInsertSizeMetrics I={input.bam} O={output.txt} HISTOGRAM_FILE={output.pdf}"


# Collect some HS metrics using picard
rule collect_hs_metrics:
	input:
		bam = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bam",
		index = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bai",
		intervals = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".interval_list"
	output:
		"output/{sample_name}/qc_reports/hs_metrics/{sample_name}_{sample_number}.txt"
	params:
		ref = config['reference']
	shell:
		"picard CollectHsMetrics I={input.bam} O={output} R={params.ref} "
		"BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals}"

# Collect alignment summary metrics using picard
rule collect_alignment_metrics:
	input:
		bam = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bam",
		index = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		"output/{sample_name}/qc_reports/alignment_metrics/{sample_name}_{sample_number}.txt"
	params:
		ref = config['reference']
	shell:
		"picard CollectAlignmentSummaryMetrics I={input.bam} O={output} R={params.ref}"



#-----------------------------------------------------------------------------------------------------------------#
# Variant Calling
#-----------------------------------------------------------------------------------------------------------------#


# Create GVCF for each sample
rule create_gvcfs:
	input:
		bam_file = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bam",
		bam_index = "output/{sample_name}/merged/{sample_name}_{sample_number}_merged_nodups.bai"
	output:
		gvcf_file = "output/{sample_name}/gvcfs/{sample_name}_{sample_number}.g.vcf"
	params:
		ref = config["reference"],
		dbsnp_vcf = config["dbsnp_vcf"],
		bed_file = config["capture_bed_file"],
		padding = config["interval_padding"]
	shell:
		"gatk HaplotypeCaller -R {params.ref} "
		"-I {input.bam_file} "
		"--emit-ref-confidence GVCF "
		"--dbsnp {params.dbsnp_vcf} "
		"-O {output} "
		"-L {params.bed_file} --interval-padding {params.padding}"


# Consolidate all samples into a genomics db for joint genotyping
rule create_genomics_db:
	input:
		gvcfs = expand("output/{sample}/gvcfs/{sample}_{sample_number}.g.vcf",
		 zip,
		 sample=initial_files,
		 sample_number=sample_numbers),
	output:
		directory("output/genomicsdb")
	params:
		files = lambda wildcards, input: " -V ".join(input)
	shell:
		"gatk GenomicsDBImport -V {params.files} --genomicsdb-workspace-path {output} "
		"-L X "

# Genotype the gvcfs and produce a joint vcf
rule genotype_gvcfs:
	input:
		directory("output/genomicsdb")
	output:
		"output/jointvcf/{worksheet}.vcf"
	params:
		ref = config["reference"],
		bed_file = config['capture_bed_file'],
		padding = config['interval_padding']
	shell:
		"gatk GenotypeGVCFs -R {params.ref} "
		"-V gendb://{input} "
		"-G StandardAnnotation "
		"-O {output} "
		"-L {params.bed_file} --interval-padding {params.padding}"


# Use hard filtering on quality attributes
# See https://gatkforums.broadinstitute.org/gatk/discussion/6925
# for more information on the values chosen here.
rule hard_filter_vcf:
	input:
		"output/jointvcf/{worksheet}.vcf"
	output:
		"output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf"
	params:
		ref = config["reference"],
		min_QD = config['min_QD'],
		max_FS = config['max_FS'],
		max_SOR = config['max_SOR'],
		min_MQ = config['min_MQ'],
		min_MQRankSum = config['min_MQRankSum'],
		min_ReadPosRankSum = config['min_ReadPosRankSum']
	shell:
		"gatk VariantFiltration -R {params.ref} -O {output} "
		"--variant {input}  --filter-expression 'QD < {params.min_QD}' --filter-name 'LOW_QD' "
		"--filter-expression 'FS > {params.max_FS}' --filter-name 'HIGH_FS' "
		"--filter-expression 'SOR > {params.max_SOR}' --filter-name 'HIGH_SOR' "
		"--filter-expression 'MQ < {params.min_MQ}' --filter-name 'LOW_MQ' "
		"--filter-expression 'MQRankSum < {params.min_MQRankSum}' --filter-name 'LOW_MQRankSum' "
		"--filter-expression 'ReadPosRankSum < {params.min_ReadPosRankSum}' --filter-name 'LOW_ReadPosRankSum' "

# compress and index vcf so we can annotate with gene name
rule compress_and_index_vcf:
	input:
		"output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf"
	output:
		vcf = "output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf.gz",
		index = "output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf.gz.tbi"
	shell:
		"bgzip {input} && tabix {output.vcf}"

# compress and index bedfile so we can use bcftools to annotate with gene
rule compress_and_index_bed_file:
	input:
		config["capture_bed_file"]
	output:
		bed = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".bed.gz",
		index = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".bed.gz.tbi"
	shell:
		"bgzip --stdout {input} > {output.bed} && tabix {output.bed}"


# Annotate with gene bed file.
rule annotate_vcf_with_gene:
	input:
		vcf = "output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf.gz",
		vcf_index = "output/qfiltered_jointvcf/{worksheet}_qfiltered.vcf.gz.tbi",
		bed = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".bed.gz",
		bed_index = "output/config/" + Path(config["capture_bed_file"]).name.split(".")[0] + ".bed.gz.tbi"
	output:
		vcf = "output/qfiltered_jointvcf_gene/{worksheet}_qfiltered_gene.vcf",
	shell:
		"bcftools annotate -a {input.bed} "
		"-c CHROM,FROM,TO,GENE "
		"-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">') "
		"{input.vcf} -o {output.vcf} "


# Select only Biallelic SNPs which pass filtering i.e. exclude indels and multialleleics and fails
rule select_relevant_variants:
	input:
		"output/qfiltered_jointvcf_gene/{worksheet}_qfiltered_gene.vcf",
	output:
		"output/qfiltered_jointvcf_gene_selected/{worksheet}_qfiltered_gene_selected.vcf"
	params:
		ref = config["reference"]
	shell:
		"gatk SelectVariants "
		"-R {params.ref} "
		"-V {input} "
		"-O {output} "
		"-select-type SNP "
		"--restrict-alleles-to BIALLELIC "
		"--exclude-filtered true "

# Split each vcf by family as specified in the config file
rule split_vcf_by_family:
	input:
		"output/qfiltered_jointvcf_gene_selected/{worksheet}_qfiltered_gene_selected.vcf"
	output:
		expand("output/family_vcfs/{{worksheet}}_qfiltered_selected_{FAMID}.vcf", FAMID=config['families'].keys())
	params:
		ref = config["reference"],
		config = config_location
	shell:
		"python utils/pipeline_scripts/split_joint_vcf_by_family.py "
		"--input {input} "
		"--config {params.config} "
		"--ref {params.ref} "
		"--output_dir output/family_vcfs "

#Create family CSVs from VCFs
rule create_family_csv:
	input:
		"output/family_vcfs/{worksheet}_qfiltered_selected_{FAMID}.vcf"
	output:
		"output/family_csvs/{worksheet}_qfiltered_selected_{FAMID}.csv"
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
		"-GF DP"


rule final:
	input:
		"output/family_csvs/{worksheet}_qfiltered_selected_FAM001.csv".format(worksheet=worksheet),
		"output/{sample_name}/qc_reports/fastq_screen/{sample_name}_{sample_number}_L001_R1_001_screen.html",
		"output/{sample_name}/qc_reports/fastqc/{sample_name}_{sample_number}_L001_R1_001.qfilter_fastqc.html",
		"output/{sample_name}/qc_reports/insert_size_metrics/{sample_name}_{sample_number}.txt",
		#"output/{sample_name}/qc_reports/hs_metrics/{sample_name}_{sample_number}.txt",
		#"output/{sample_name}/qc_reports/alignment_metrics/{sample_name}_{sample_number}.txt"
	output:
		"output/finished/{sample_name}_{sample_number}_finished.txt"
	shell:
		"touch {output}"
