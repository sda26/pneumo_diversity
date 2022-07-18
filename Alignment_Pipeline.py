#!/usr/bin/env python3
#

import os
import sys
import subprocess
import glob
import shutil

print ('The current working directory is: ' + os.getcwd())

# reading manifest file and running bowtie2
# fil= open("manifestfile.csv")

# index = Index reference genome for alignment
subprocess.call("bwa/bwa index reference/pE539_amplicon.fasta", shell=True)

# Remove Illumina adapters and trim reads to remove low-quality bases
subprocess.call('java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE reads/SDA06/SDA06_R1_001.fastq.gz reads/SDA06/SDA06_R2_001.fastq.gz reads/Trimmed_F_P.fastq reads/Trimmed_F_UP.fastq reads/Trimmed_R_P.fastq reads/Trimmed_R_UP.fastq ILLUMINACLIP:Trimmomatic-0.39/adapters/illumina_adapter.fa:2:30:10 SLIDINGWINDOW:3:20 LEADING:15 TRAILING:15 MINLEN:75', shell=True)

# mem = align reads to reference genome and generate SAM file, -O = gap penalty, so 100 restricts gaps.
# gzip = bgzip compression. -4 = intermediate quality compression (-1 is worst/fastest, -9 is slowest/best)
subprocess.call('bwa/bwa mem -O 100 -A 10 -B 2 reference/pE539_amplicon.fasta reads/Trimmed_F_P.fastq reads/Trimmed_R_P.fastq > output/SDA06/joined_iga.sam', shell = True)

# The R script for calling haplotypes from the .sam file can be incorporated into this pipeline, but I would advise running it manually at first, until confident.
#subprocess.call('Rscript Haplotype_Calling.R', shell = True)

#subprocess.call('Rscript Haplotype_Analysis.R', shell = True)