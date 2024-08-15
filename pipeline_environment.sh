#Path to PseudoreferencePipeline (https://github.com/YourePrettyGood/PseudoreferencePipeline):
#Used for mapping and marking duplicates
export PSEUDOREF="/home/pfr8/bin/PseudoreferencePipeline"

#All of these variables must contain the full path to the executables
#This pipeline depends on:
#Samtools (compatible with versions 1.0+, may have issues with older)
#Picard (can obtain via `git clone https://github.com/broadinstitute/picard --recursive`)

export BBMERGE="/home/pfr8/bin/bbmap/bbmerge.sh"
export SEQPREP="/home/pfr8/bin/SeqPrep/SeqPrep"
export FASTQC="/ysm-gpfs/apps/software/FastQC/0.11.8-Java-1.8/fastqc"
export FASTQSCREEN="/home/pfr8/bin/FastQ-Screen/fastq_screen"
export KRAKEN="/home/pfr8/bin/kraken2/kraken2"
export BRACKEN="/home/pfr8/bin/Bracken/bracken"
export BWA="/home/pfr8/bin/bwa/bwa"
export SAMTOOLS="/home/pfr8/bin/samtools/bin/samtools"
export PICARD="/home/pfr8/bin/picard/picard-2.23.3.jar"
export MUTSERVE="/home/pfr8/bin/mutserve/mutserve"
export HAPLOGREP="/home/pfr8/bin/haplogrep-cmd/haplogrep"
#For Rscript, load the module:
module load R/4.1.0-foss-2020b

#Configuration files and databases:
#FastQ-Screen configuration file:
#This configuration is generic for human samples:
#export SCREENCONF="/home/pfr8/bin/FastQ-Screen/fastq_screen_human_contam.conf"
#Use this configuration for Plasmodium-focused screening:
export SCREENCONF="/gpfs/gibbs/pi/tucci/pfr8/bin/QCPipeline/fastq_screen_human_plasmodium.conf"
#Kraken2 pre-compiled database:
#For Standard (built by PFR):
#export KRAKEN2DB="/gpfs/gibbs/pi/tucci/pfr8/Kraken2_DB/"
#For Standard+Protozoa+Fungi+Plants (precompiled with Bracken databases):
#2020/12/02 PlusPFP database actually missed Fungi:
#export KRAKEN2DB="/gpfs/gibbs/pi/tucci/pfr8/Kraken2_DB/PlusPFP/"
#Instead, the 2021/01/28 PlusPFP database is fixed:
export KRAKEN2DB="/gpfs/gibbs/pi/tucci/pfr8/Kraken2_DB/PlusPFP_20210128/"
#mtDNA reference FASTA for mutserve:
export RCRS="/home/pfr8/bin/mutserve/rCRS.fasta"
