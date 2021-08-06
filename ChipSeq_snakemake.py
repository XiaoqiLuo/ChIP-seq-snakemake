############################################################################################
# 2021-08-02
# Xiaoqi Luo
# WES pipeline
############################################################################################

# snakemake -s ChipSeq_snakemake.py -p -j 1 --config workspace=/mnt/d/lxq/Training/Chip-Seq \
# genomes=/mnt/d/lxq/Training/WES/GATK/hg38/bwa_index \
# blacklist=/mnt/d/lxq/Training/reference
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
		index=config['genomes']+'gatk_hg38', #genomes: D:\lxq\Training\WES\GATK\hg38\bwa_index
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
        'bedtools intersect -v -abam {input.bam} -b {params.blist} > {outputa}'

rule buildindex:
    output:
        config['workspace'] + '/mapped/{sample}_rmblacklist.bam.bai'
    input:
        config['workspace'] + '/mapped/{sample}_rmblacklist.bam'
    shell:
        'samtools index {input} {output}'

rule  chromsize:
    output:
        config['workspace'] + '/mapped/{sample}_chromosome_size.bed'
    input:
        config['workspace'] + '/mapped/{sample}_rmblacklist.bam'
    shell:
        'samtools idxstats {input} | grep -v \'_\' | grep -v \'-\' | grep -v chrM'  
        '| awk \'{{print $1, 0, $2}}\' > {output}'

rule cleanbam:
    output:
        config['workspace'] + '/mapped/{sample}_clean.sorted.bam'
    input:
        bam=config['workspace'] + '/mapped/{sample}_rmblacklist.bam',
        bed=config['workspace'] + '/mapped/{sample}_chromosome_size.bed'
    shell:
        'samtools view -b -f 2 -F 1024'
        ' -L {input.bed} {input.bam} > {output}'

rule cleanindex:
    output:
        config['workspace'] + '/mapped/{sample}_clean.sorted.bam.bai'
    input:
        config['workspace'] + '/mapped/{sample}_clean.sorted.bam'
    shell:
        'samtools index {input} {output}'
        
rule callpeak:
    output:
        config['workspace'] + '/call_peaks/{sample}_peaks.narrowPeak'
    input:
        config['workspace'] + '/mapped/{sample}_clean.sorted.bam'
    params:
        resultdir=config['workspace'] + '/call_peaks',
        name='{sample}'
    shell:
        'macs2 callpeak -p 1e-2 -f BAMPE --keep-dup all -B -n {params.name} -t {input} --outdir {params.resultdir}'
