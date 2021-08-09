# Chip-Seq-snakemake

## Conda Enviornment
```
conda install -c bioconda -c conda-forge bowtie sra-tools samtools fastp bowtie2 bwa bowtie bedtools snakemake picard
```

## GATK Download
```
wget  https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
unzip gatk-4.0.6.0.zip
cd /GATK/gatk-4.0.6.0
./gatk --help
```

## Files needed for analysis
### blacklist
```
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz
```

## How To Execute
### File Tree (fastq files are saved in rawfastq file)
>workspace<br>
>>ChipSeq_snakemake.py<br>
>>rawfastq<bc>
>>>SRR539135_1.fastq.gz<br>
>>>SRR539135_2.fastq.gz<br>
>>>SRR13242842_1.fastq.gz<br>
>>>SRR13242842_2.fastq.gz<br>

```
snakemake -s ChipSeq_snakemake.py -p -j 1 --config workspace=path-to-workspace \
genomes=path-to-bwa-index \
blacklist=path-to-blacklist-file \
GATK=path-to-gatk 
```
 
### example
```
snakemake -s ChipSeq_snakemake.py -p -j 1 --config workspace=/mnt/d/lxq/Training/Chip-Seq \
genomes=/mnt/d/lxq/Training/WES/GATK/hg38/bwa_index \
blacklist=/mnt/d/lxq/Training/reference \
GATK=/mnt/d/lxq/Training/WES/GATK/gatk-4.1.7.0/gatk 
```
## result
The pipeline generates 17 files for one sample, which ***call_peaks*** file is the result file needed for downstream analysis.
File tree shown below:<br>
![image](https://github.com/XiaoqiLuo/Chip-Seq-snakemake/blob/main/Chipseq-workspace.png)
