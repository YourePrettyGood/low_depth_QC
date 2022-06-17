# Low-depth QC pipeline

This repository is a cluster-oriented pipeline for QC of low-depth
sequencing data, mainly written for human genomic data. The pipeline format
follows the general format of the [PseudoreferencePipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline),
but I will eventually make a Nextflow version of it. Template SBATCH scripts
are provided for the SLURM case, and two setup scripts are also provided
to simplify the process of FASTQ organization and SBATCH script generation.

## What does the pipeline actually do?

This pipeline assesses the quality of your input libraries in a number of
different ways:

1. FastQC (looking for technical issues in the reads, like flowcell failures)
1. Library characteristics (insert size distributions for merged and unmerged reads)
1. Mapping rates to a selected reference (e.g. hs37d5, GRCh38DH, etc.)
1. Contamination detection (FastQ-Screen and kraken2+Bracken)
1. mtDNA and Y haplogroup estimation (mutserve+haplogrep and Yleaf)
1. Genetic sex estimation (both Y-restricted estimators like Skoglund's R_Y and sex chromosome ploidy estimators)

## Dependencies:

1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (tested with version 0.11.8)
1. [FastQ-Screen](https://github.com/StevenWingett/FastQ-Screen) (tested with commit 3168d23)
1. Perl's GD::Graph module (tested with version 1.54 using Perl 5.28.0)
1. [SeqPrep](https://github.com/jstjohn/SeqPrep) for adapter trimming and read pair merging (tested with commit 575507b)
1. [kraken2](https://github.com/DerrickWood/kraken2) for contamination screening and the [PlusPFP database from here](https://benlangmead.github.io/aws-indexes/k2) (tested with commit 00ceda0)
1. NCBI BLAST+ (dependency of kraken2, tested with version 2.8.1)
1. [Bracken](https://github.com/jenniferlu717/Bracken) for contamination abundance estimates (tested with commit f420218, which is version 2.5.3)
1. [BWA-MEM](https://github.com/lh3/bwa) for read mapping (tested with commit 13b5637, aka 0.7.17-r1198)
1. [SAMtools](https://github.com/samtools/samtools) for BAM processing and QC (tested with commits af811a6 and cc4e1a6)
1. [HTSlib](https://github.com/samtools/htslib) as a dependency of SAMtools
1. [mutserve](https://github.com/seppinho/mutserve) for mtDNA variant calling (tested 2.0.0-rc12)
1. [haplogrep](https://github.com/seppinho/haplogrep-cmd) for mtDNA haplogroup classification (tested 2.4.0)
1. [Yleaf](https://github.com/genid/Yleaf) for Y haplogroup inference (tested commit 7c33ca0)
1. [Picard](https://github.com/broadinstitute/picard) for various steps (tested 2.24.0)
1. [R](https://r-project.org) for some plotting
1. tidyverse packages for R (`install.packages(c("tidyverse"))`) for the R script

## Setup of reference files:

It's important to set up the reference and associated files appropriately
before starting the pipeline. Not only do you need to set up the reference
for the mapping rate and insert size sections, but you also need to index
any references needed by FastQ-Screen.

TODO: Fill in more detail here

## How to run the pipeline

I've provided a script `prep_metadata.sh` that creates the metadata
file used by the pipeline to run each sample in parallel.

TODO: Fill in more detail here

### Config files



## Results/outputs

Currently, the pipeline outputs the following files and directories
for each sample:

`logs/`: Logs for each major step

`[Sample ID]_FastQC/`: HTML and ZIP files of the reports from FastQC

`[Sample ID]_fastq_screen/`: HTML and TXT files of the reports from FastQ-Screen

TODO: Fill in more info here

## Extra scripts

