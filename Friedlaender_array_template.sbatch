#!/bin/bash
#SBATCH -J F[BATCH][SHORTREF]QC
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=12
#Make sure to specify enough memory for loading the kraken2 index
#For PlusPFP, that's more than 90 GB
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -e logs/Friedlaender_[BATCH]_[REF]_QC_%j.stderr
#SBATCH -o logs/Friedlaender_[BATCH]_[REF]_QC_%j.stdout
#SBATCH --mail-user [EMAIL]
#SBATCH --mail-type ALL

#SBATCH -a 1-[BATCHSIZE]

#Match this to --cpus-per-task:
cores=12

PIPELINEDIR="/home/pfr8/bin/QCPipeline"

#Metadata file replaces the need to determine the sample files above:
metadata="Friedlaender_[BATCH]_[REF]_QC_metadata.tsv"

#Special options include:
#no_markdup: Skip Picard MarkDuplicates
#BBMerge: Use BBMerge instead of SeqPrep for merging OLPE read pairs
#subset_#: Set FastQ-Screen's --subset option to #
#adapter1_[ACGT]+: Set adapter1 for trimming/merging
#adapter2_[ACGT]+: Set adapter2 for trimming/merging
special="[SPECIAL]"

pushd /gpfs/gibbs/pi/tucci/pfr8/Friedlaender/[BATCH]/

/usr/bin/time -v ${PIPELINEDIR}/array_call.sh ${SLURM_ARRAY_TASK_ID} ${metadata} ${cores} ${special}

popd
