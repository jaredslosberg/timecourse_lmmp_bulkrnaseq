#!/bin/bash

fastq_path="/mnt/morbo/Data/Public/kutlu_agudelo_dio_2021/raw/fastqs"
srr_path="/mnt/morbo/Data/Public/kutlu_agudelo_dio_2021/raw/SRP340286_run_accession_list.txt"

index_path="/reference/indexes/kallisto_indexes/mouse/mm_vM27_index.idx"
output_path="/mnt/morbo/Data/Public/kutlu_agudelo_dio_2021/kallisto_out_pseudobam"

cd $fastq_path || exit

source activate kallistobus

kallisto version

#kallisto quant -i $index_path -o $output_path --single --pseudobam -l 120 
parallel -j 2 -a ${srr_path} \
    kallisto quant --pseudobam -i ${index_path} -o ${output_path}/{} -t 4 {}_1.fastq {}_2.fastq


