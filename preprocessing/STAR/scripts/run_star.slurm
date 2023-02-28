#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=STAR_align
#SBATCH --output run_star.out
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200gb


sample_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/samples_merged.txt"
fastq_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs/merged"
output_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/preprocessing/STAR"
#add STAR to path
export PATH=$PATH:/usr/bin/STAR/source

#build index - only need to run once
# STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/users/jared/STAR_output \
#   --genomeFastaFiles /reference/genomes/mouse/gencode/vM27/GRCm39.primary_assembly.genome.fa \
#   --sjdbGTFfile /reference/genomes/mouse/gencode/vM27/gencode.vM27.annotation.gtf \
# --limitGenomeGenerateRAM 100000000000

cd $fastq_path || exit

#setup directories
parallel --tmpdir /scratch/users/jared -j 5 -a ${sample_path} \
  mkdir -p ${output_path}/{}
  
#run mapping
parallel --tmpdir /scratch/users/jared -j 4 -a ${sample_path} \
  STAR --runThreadN 4 --genomeDir /data/users/jared/STAR_output \
  --outFileNamePrefix=${output_path}/{} \
  --readFilesIn {}_R1.fq {}_R2.fq

#     #run mapping
# STAR --runThreadN 10 --genomeDir /data/users/jared/STAR_output \
# --outFileNamePrefix=/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/preprocessing/STAR/ \
# --readFilesCommand gunzip --readFilesIn /mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs/17mo-F-1_S6_L001_R1_001.fastq.gz \
#/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs/17mo-F-1_S6_L001_R2_001.fastq.gz
