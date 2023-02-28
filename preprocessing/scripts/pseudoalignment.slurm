#!/bin/bash
#SBATCH --account=jared
#SBATCH --job-name=Pseudoalign
#SBATCH --output kallisto_bus_pseudoalign.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400gb


##---------------------------------------------------------------------
## SLURM job script to run kallisto pseudoalign
##---------------------------------------------------------------------

date

source  ~/anaconda3/etc/profile.d/conda.sh

fastq_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/fastqs"
sample_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/samples.txt"

index_path="/reference/indexes/kallisto_indexes/mouse/mm_vM27_index.idx"
output_path="/mnt/morbo/Data/Users/jslosberg/aged_bulkrna/kallisto_out/"

cd $fastq_path || exit

conda activate kallistobus

kallisto version

#kallisto quant -i $index_path -o $output_path --single --pseudobam -l 120 
parallel --tmpdir /scratch/users/jared -j 4 -a ${sample_path} \
    kallisto quant -i ${index_path} -o ${output_path}/{} -t 4 {}_R1_001.fastq.gz {}_R2_001.fastq.gz


