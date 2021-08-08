############################################################################################
# 2021-08-02
# Xiaoqi Luo
# WES pipeline
############################################################################################

# snakemake -s ChipSeq_snakemake.py -p -j 1 --config workspace=/mnt/d/lxq/Training/Chip-Seq \
# genomes=/mnt/d/lxq/Training/WES/GATK/hg38/bwa_index \
# blacklist=/mnt/d/lxq/Training/reference \
# GATK=/mnt/d/lxq/Training/WES/GATK/gatk-4.1.7.0/gatk 

files=os.listdir(config['workspace']+'/rawfastq/')
SRR = []
for f in files[0:3:2]:
    SRR.append(re.split(r'_', f)[0])

rule all:
    input:
        expand(
            config['workspace'] + '/call_peaks/{sample}_peaks.narrowPeak',
            sample=SRR
		)
  
rule quality:
	output:
		fq1=config['workspace'] + '/qc/{sample}_1_trim.fastq.gz',
		fq2=config['workspace'] + '/qc/{sample}_2_trim.fastq.gz',
		json=config['workspace'] + '/qc/{sample}_fastp.json',
		html=config['workspace'] + '/qc/{sample}_fastp.html'
	input:
		fq1=config['workspace'] + '/rawfastq/{sample}_1.fastq.gz',
		fq2=config['workspace'] + '/rawfastq/{sample}_2.fastq.gz'
	log:
		config['workspace'] + '/logs/{sample}_trim.log'
	shell:
		'fastp -i {input.fq1} -I {input.fq2} '
		'-o {output.fq1} -O {output.fq2} '
		'-g -q 5 -u 50 -n 15 '
		'-j {output.json} -h {output.html} 2>{log}'

rule mapping:
	output:
		config['workspace'] + '/mapped/{sample}.bam'
	input:
		fq1=config['workspace'] + '/qc/{sample}_1_trim.fastq.gz',
		fq2=config['workspace'] + '/qc/{sample}_2_trim.fastq.gz'
	log:
		config['workspace'] + '/logs/{sample}_bwa_mem.log'
	params:
		index=config['genomes']+'/gatk_hg38', #genomes: D:\lxq\Training\WES\GATK\hg38\bwa_index
	shell:
		'bwa mem {params.index} -v 1 -T 30 -h 5'
		'{input.fq1} {input.fq2} '
		'| samtools sort -o {output} -'

rule rmblacklist:
    output:
        config['workspace'] + '/mapped/{sample}_rmblacklist.bam'
    input:
        bam=config['workspace'] + '/mapped/{sample}.bam'
    params:
        blist=config['blacklist'] + '/wgEncodeHg19ConsensusSignalArtifactRegions.bed'
    shell:
        'bedtools intersect -v -abam {input.bam} -b {params.blist} > {output}'

rule markdup:
	output:
		bam=config['workspace'] + '/mapped/{sample}_marked.bam',
		metrics=config['workspace'] + '/mapped/{sample}_marked.metrics.txt'
	input:
		config['workspace'] + '/mapped/{sample}.bam'
	log:
		config['workspace'] + '/logs/{sample}_markdup.log'
	params:
		extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
		GATK=config['GATK'] 
	shell:
		'{params.GATK} --java-options {params.extra} MarkDuplicates '
		'-I {input} -O {output.bam} -M {output.metrics}'

#The output bam file was then filtered to remove unmapped reads and reads with Mapping Quality less than 5
rule filter_unmap:
    output:
        config['workspace'] + '/mapped/{sample}_rmblacklist_marked_intermediate.bam'
    input:
        config['workspace'] + '/mapped/{sample}_rmblacklist_marked.bam'
    shell:
        'samtools view -b -F 4 -q 5 {input} > {output}'

#The intermediate output bam file was then filtered to remove PCR or optical duplicate reads
rule filter_PCR:
    output:
        config['workspace'] + '/mapped/{sample}_output.bam'
    input:
        config['workspace'] + '/mapped/{sample}_rmblacklist_marked_intermediate.bam'
    shell:
        'samtools view -b -F 1024 {input} > {output}'

rule callpeak:
    output:
        config['workspace'] + '/call_peaks/{sample}_peaks.narrowPeak'
    input:
        config['workspace'] + '/mapped/{sample}_output.bam'
    params:
        resultdir=config['workspace'] + '/call_peaks',
        name='{sample}'
    shell:
        'macs2 callpeak -p 1e-5 --keep-dup all -B -n {params.name} -t {input} --outdir {params.resultdir}'


# macs2 callpeak -p 1e-3 --keep-dup all -B -n SRR13242842 -t SRR13242842_output.bam --outdir ./call_peak
