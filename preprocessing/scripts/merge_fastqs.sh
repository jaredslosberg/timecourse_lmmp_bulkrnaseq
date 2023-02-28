#!/bin/bash

fastq_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs"
sample_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/samples_merged.txt"

output_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs/merged"
mkdir -p $output_path

cd $fastq_path || exit

echo ${output_path}

parallel --tmpdir /scratch/users/jared -j 9 -a ${sample_path} \
   cat {}*R1_001.fastq.gz ">" ${output_path}/{}_R1.fq.gz
   
parallel --tmpdir /scratch/users/jared -j 9 -a ${sample_path} \
   cat {}*R2_001.fastq.gz ">" ${output_path}/{}_R2.fq.gz


