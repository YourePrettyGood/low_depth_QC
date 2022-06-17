#!/bin/bash

#This script uses a metadata file to run a pipeline job
# when submitted as bash or GNU parallel job

#The arguments are:
#1) Task ID (line of the metadata file to use)
#2) metadata file (TSV comprised of prefix, ref, read file(s))
#4) Optional override of # cores used
#5) special parameters (e.g. no_markdup, etc.)

#PFR notes to avoid bugs:
#outputdir variable *MUST* add the slash to avoid attempting to
# write to root directory when prefix is local
#prefix for indexAlignDictMarkDup_v2.sh should be the same as
# the input prefix, so prepended ${outputdir} is included
#The whole IFS=$'\a' and substituted tr "\t" "\a" enforces
# tab separation in the metadata file, but maybe we should
# be more flexible with that and just let bash do it's
# default whitespace splitting
#One error I ran into was accidentally compiling the dependencies
# with different versions of GCC (and libgcc, etc.) due to stale
# modules during compilation. So if you get errors saying stuff
# like "GLIBCXX_3.4.21 not found", that's the likely cause.

task_id=$1
metadata=$2
if [[ ! -e "${metadata}" ]]; then
   echo "The metadata file does not exist! Did you make a typo? The file you specified is: ${metadata}"
   exit 2;
fi
if [[ ! -z "$3" ]]; then
   cores="$3"
fi
special=$4

whichsample=1
while IFS=$'\a' read -r -a metadatafields
   do
   if [[ ${whichsample} -eq ${task_id} ]]; then
      if [[ ${#metadatafields[@]} -lt "3" ]]; then
         echo "No read files supplied! Is your metadata file actually tab-separated (not space-separated)?"
         exit 8;
      fi
      prefix="${metadatafields[0]}"
      ref="${metadatafields[1]}"
      if [[ ! -e "${ref}" ]]; then
         echo "Reference ${ref} does not exist! Did you make a typo?"
         exit 5;
      fi
      if [[ ! -z "${metadatafields[3]}" ]]; then
         read1="${metadatafields[2]}"
         read2="${metadatafields[3]}"
#         reads="${metadatafields[2]} ${metadatafields[3]}"
         if [[ ! -e "${metadatafields[3]}" ]]; then
            echo "Read 2 file ${metadatafields[3]} does not exist! Did you make a typo?"
            exit 7;
         fi
      else
         read1="${metadatafields[2]}"
#         reads="${metadatafields[2]}"
         echo "Only 1 read file provided for sample ${prefix}, this pipeline expects paired-end data."
         echo "Please check your metadata file, as I read this for line ${task_id}:"
         echo "Sample Prefix = ${prefix}"
         echo "Reference Genome = ${ref}"
         echo "Read files = ${read1}"
         exit 2;
      fi
      if [[ ! -e "${metadatafields[2]}" ]]; then
         echo "Read 1 file ${metadatafields[2]} does not exist! Did you make a typo?"
         exit 6;
      fi
   fi
   (( whichsample++ ))
done < <(tr "\t" "\a" < ${metadata})
if [[ -z "${prefix}" ]]; then
   echo "Unable to find sample ${task_id} in metadata file. Skipping."
   exit 4
fi

#Determine the directory of the array_call.sh script so we can source
# the pipeline_environment.sh script for useful paths
scriptdir=`dirname $0`

#Now source pipeline_environment.sh
source ${scriptdir}/pipeline_environment.sh

#Check for the existence of the programs necessary for this pipeline:
#bbmerge.sh
if [[ ! -e "${BBMERGE}" ]]; then
   echo "Could not find bbmerge.sh at BBMERGE=${BBMERGE}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#SeqPrep
if [[ ! -x "$(command -v ${SEQPREP})" ]]; then
   echo "Could not find SeqPrep at SEQPREP=${SEQPREP}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#fastqc
if [[ ! -e "${FASTQC}" ]]; then
   echo "Could not find fastqc at FASTQC=${FASTQC}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#fastq_screen
if [[ ! -e "${FASTQSCREEN}" ]]; then
   echo "Could not find fastq_screen at FASTQSCREEN=${FASTQSCREEN}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#kraken2
if [[ ! -e "${KRAKEN}" ]]; then
   echo "Could not find kraken2 at KRAKEN=${KRAKEN}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#bracken
if [[ ! -e "${BRACKEN}" ]]; then
   echo "Could not find bracken at BRACKEN=${BRACKEN}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#bwa
if [[ ! -x "$(command -v ${BWA})" ]]; then
   echo "Could not find bwa at BWA=${BWA}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#samtools
if [[ ! -x "$(command -v ${SAMTOOLS})" ]]; then
   echo "Could not find samtools at SAMTOOLS=${SAMTOOLS}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#Picard
if [[ ! -e "${PICARD}" ]]; then
   echo "Could not find picard.jar at PICARD=${PICARD}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#mutserve
if [[ ! -x "$(command -v ${MUTSERVE})" ]]; then
   echo "Could not find mutserve at MUTSERVE=${MUTSERVE}, please fix the path in pipeline_environment.sh"
   exit 3
fi
if [[ ! -e "${RCRS}" ]]; then
   echo "Could not find mtDNA reference FASTA for mutserve at RCRS=${RCRS}, please fix the path in pipeline_environment.sh"
   exit 3
fi
#haplogrep
if [[ ! -x "$(command -v ${HAPLOGREP})" ]]; then
   echo "Could not find haplogrep at HAPLOGREP=${HAPLOGREP}, please fix the path in pipeline_environment.sh"
   exit 3
fi
##Rscript


#Make a logs directory if it doesn't already exist
mkdir -p logs

#Parse out the output directory from the sample prefix:
outputdir=""
if [[ ${prefix} =~ \/ ]]; then #If the prefix has a path
   outputdir="`dirname ${prefix}`/"
   prefix=`basename ${prefix}`
fi

#Get a short name for the reference
#Note: This is fairly custom for human, but can be general if
# you consistently name the reference FASTA as [short name]_[rest]
# or [short name].[extension]
reffn=`basename ${ref}`
refid=${reffn%%[_.]*}
#Also determine the name of the mtDNA contig based on this short name:
#This is used mainly for mutserve/haplogrep and sexing.
if [[ "${refid}" == "hs37d5" ]]; then
   #hs37d5 names chromosomes with their number/letter, so mtDNA is MT
   mt="MT"
   #And the rest are sex chromosomes and autosomes with similar lengths:
   x="X"
   y="Y"
   autox="7"
   autoy="19"
elif [[ "${refid}" == "GRCh38" ]]; then
   #GRCh38 names chromosomes as chr([0-9]{1,2}|[XYM]), so mtDNA is chrM
   mt="chrM"
   #And the rest are sex chromosomes and autosomes with similar lengths:
   x="chrX"
   y="chrY"
   autox="chr7"
   autoy="chr19"
elif [[ "${refid}" == "Kamilah" ]]; then
   #Kamilah is the current NCBI RefSeq Gorilla gorilla gorilla reference
   #So all of the chromosome IDs are NC_* accessions
   #Note that Kamilah is female, so there isn't a Y chromosome in the assembly
   mt="NC_011120.1"
   x="NC_044625.1"
   y=""
   autox="NC_044609.1"
   autoy=""
else
   echo "Reference ${refid} not recognized, was the name ${reffn} parsed incorrectly?"
   exit 4
fi

#Java Memory for things like Picard and BBMerge:
#For initial testing, use 4 GB RAM in the heap
javamem="4g"
javamemre='(mem)_([0-9]+[MGmg]?)'
if [[ ${special} =~ $javamemre ]]; then
   javamem="${BASH_REMATCH[2]}"
fi
echo "Running java steps with -Xmx${javamem} for sample ${sample}"

#Selection of merging program:
if [[ "${special}" =~ "BBMerge" ]]; then
   merger="BBmerge"
else
   merger="SeqPrep"
fi

#QC summary file and header:
qc_summary="${outputdir}${prefix}_${refid}_${merger}_QC_summary.tsv"
printf "Sample\tRef\tMerger" > ${qc_summary} #metadata
printf "\tRawReadPairs\tRawBases" >> ${qc_summary} #awk raw
printf "\tMergedReads\tMergedBases" >> ${qc_summary} #awk M
printf "\tUnmergedReadPairs\tUnmergedBases" >> ${qc_summary} #awk UM
printf "\tClassifiedMerged\tUnclassifiedMerged\tPctClassifiedMerged" >> ${qc_summary} #kraken2 M
printf "\tHumanAbundanceMerged" >> ${qc_summary} #Bracken M
printf "\tClassifiedUnmerged\tUnclassifiedUnmerged\tPctClassifiedUnmerged" >> ${qc_summary} #kraken2 UM
printf "\tHumanAbundanceUnmerged" >> ${qc_summary} #Bracken UM
if [[ "${refid}" == "hs37d5" || "${refid}" == "GRCh38" ]]; then
   printf "\tmtHaplogroupMerged\tmtHaplogroupQualityMerged" >> ${qc_summary} #Haplogrep M
   printf "\tmtHaplogroupUnmerged\tmtHaplogroupQualityUnmerged" >> ${qc_summary} #Haplogrep UM
fi
printf "\tmtDepthMerged\tmtDepthUnmerged" >> ${qc_summary} #mtDNA depth
if [[ ! -z "${y}" ]]; then
   printf "\tYDepthMerged\tAutoYDepthMerged\tYNormalizedDepthMerged" >> ${qc_summary} #Y normalized depth M, could be third way to estimate sex
   printf "\tYDepthUnmerged\tAutoYDepthUnmerged\tYNormalizedDepthUnmerged" >> ${qc_summary} #Y normalized depth UM, as above
   printf "\tYMappedMerged\tTotalMappedMerged\tYFractionMerged\tYEstimatedSexMerged" >> ${qc_summary} #Serena's method to estimate sex from percent of total reads mapping to Y (merged)
   printf "\tYMappedUnmerged\tTotalMappedUnmerged\tYFractionUnmerged\tYEstimatedSexUnmerged" >> ${qc_summary} #As above, but for unmerged
fi
if [[ ! -z "${x}" ]]; then
   printf "\tXDepthMerged\tAutoXDepthMerged\tXNormalizedDepthMerged" >> ${qc_summary} #PFR's method to estimate sex from normalized depth of X (see e.g. Mongue and Walters on Danaus plexippus Z (merged)
   printf "\tXDepthUnmerged\tAutoXDepthUnmerged\tXNormalizedDepthUnmerged" >> ${qc_summary} #As above, but for unmerged
   printf "\tXEstimatedSexMerged\tXEstimatedSexUnmerged" >> ${qc_summary}
fi
printf "\n" >> ${qc_summary}

#Steps in the pipeline:
#1a) Merging
#1b) FastQC on orig reads
#2a) Contaminant screening (FastQ-Screen, kraken2)
#2b) Mapping (orig, merged, unmerged)
#3a) Mapping rates
#3b) Insert length distribution
#4) mtDNA haplogroup
#5) A few average depths, including for sexing

#Merging:
#Illumina standard adapter sequences for trimming:
illumina_r1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
illumina_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
adapter1re='(adapter1)_([ACGTacgt]+)'
if [[ ${special} =~ $adapter1re ]]; then
   illumina_r1="${BASH_REMATCH[2]}"
fi
adapter2re='(adapter2)_([ACGTacgt]+)'
if [[ ${special} =~ $adapter2re ]]; then
   illumina_r2="${BASH_REMATCH[2]}"
fi
echo "Trimming and merging with adapter sequences:"
if [[ "${merger}" == "BBMerge" || "${special}" =~ "adapter[12]_" ]]; then
   echo "R1: ${illumina_r1}"
   echo "R2: ${illumina_r2}"
else
   echo "SeqPrep defaults"
fi

if [[ "${merger}" == "BBMerge" ]]; then
   #BBMerge is VERY fast and takes very little memory
   #Run BBMerge
   #Separate input FASTQs into R1 and R2:
   r1="${read1}"
   r2="${read2}"
   #Set up the output filenames:
   ihist="${outputdir}${prefix}_${merger}_insert_size_histogram.txt"
   unmergedfq1="${outputdir}${prefix}_${merger}_R1.fastq.gz"
   unmergedfq2="${outputdir}${prefix}_${merger}_R2.fastq.gz"
   mergedfq="${outputdir}${prefix}_${merger}_merged.fastq.gz"
   #Set up the options:
   opts="-Xmx${javamem} usejni=t adapter1=${illumina_r1} adapter2=${illumina_r2} t=${cores}"
   #Run the command:
   echo "Running bbmerge.sh to trim and merge ${prefix}"
   echo "BBMerge options: ${opts}"
   echo "Logs at logs/bbmerge_${prefix}.std{err,out}"
   /usr/bin/time -v ${BBMERGE} ${opts} in1=${r1} in2=${r2} outm=${mergedfq} outu1=${unmergedfq1} outu2=${unmergedfq2} ihist=${ihist} 2> logs/bbmerge_${prefix}.stderr > logs/bbmerge_${prefix}.stdout
else
   #SeqPrep is slower than BBMerge and takes more memory
   #Run SeqPrep (default)
   #Separate input FASTQs into R1 and R2:
   r1="${read1}"
   r2="${read2}"
   #Set up the output filenames:
   discardedfq1="${outputdir}${prefix}_${merger}_discarded_R1.fastq.gz"
   discardedfq2="${outputdir}${prefix}_${merger}_discarded_R2.fastq.gz"
   unmergedfq1="${outputdir}${prefix}_${merger}_R1.fastq.gz"
   unmergedfq2="${outputdir}${prefix}_${merger}_R2.fastq.gz"
   mergedfq="${outputdir}${prefix}_${merger}_merged.fastq.gz"
   #Set up the argument string:
   opts="-f ${r1} -r ${r2} -1 ${unmergedfq1} -2 ${unmergedfq2} -3 ${discardedfq1} -4 ${discardedfq2} -s ${mergedfq}"
   #Tack on adapter sequences unless not specified in SPECIAL:
   if [[ ${special} =~ ${adapter1re} || ${special} =~ ${adapter2re} ]]; then
      opts="${opts} -A ${illumina_r1} -B ${illumina_r2}"
   fi
   #Run the command:
   echo "Running SeqPrep to trim and merge ${prefix}"
   echo "SeqPrep options: ${opts}"
   echo "Logs at logs/seqprep_${prefix}.std{err,out}"
   /usr/bin/time -v ${SEQPREP} ${opts} 2> logs/seqprep_${prefix}.stderr > logs/seqprep_${prefix}.stdout
fi

#QC summary TSV beginning after header:
printf "${prefix}\t${refid}\t${merger}" >> ${qc_summary}
#Simple summary info for raw data and merging results:
gzip -dc ${r1} | awk 'NR%4==2{sum+=length($0);count++;}END{printf "\t%i\t%i", count, 2*sum;}' >> ${qc_summary}
gzip -dc ${mergedfq} | awk 'NR%4==2{sum+=length($0);count++;}END{printf "\t%i\t%i", count, sum;}' >> ${qc_summary}
gzip -dc ${unmergedfq1} | awk 'NR%4==2{sum+=length($0);count++;}END{printf "\t%i\t%i", count, 2*sum;}' >> ${qc_summary}

#FastQC:
#FastQC takes a while, but takes little memory
fastqcdir="${outputdir}${prefix}_FastQC"
mkdir -p ${fastqcdir}
echo "Running FastQC on ${prefix} original reads ${read1} and ${read2} and merging results ${mergedfq}, ${unmergedfq1}, and ${unmergedfq2} from merger ${merger}"
echo "Logs at logs/fastqc_${prefix}_${merger}.std{err,out}"
/usr/bin/time -v ${FASTQC} -t ${cores} -o ${fastqcdir} ${read1} ${read2} ${mergedfq} ${unmergedfq1} ${unmergedfq2} 2> logs/fastqc_${prefix}_${merger}.stderr > logs/fastqc_${prefix}_${merger}.stdout

#Contaminant screening (FastQ-Screen, kraken2):
#FastQ-Screen:
#FastQ-Screen is slower than FastQC since it uses the entire dataset (not default),
# and uses more memory because it actually maps the reads
#Map with BWA-MEM, use the entire dataset, use settings in config file ${SCREENCONF},
# force overwrite if --outdir exists, and outputs go to ${fqsdir}:
#Default is --subset 0, but add a special option to set it:
fqssubset=0
fqssubsetre='(subset)_([0-9]+)'
if [[ ${special} =~ $fqssubsetre ]]; then
   fqssubset="${BASH_REMATCH[2]}"
fi
fqsdir="${outputdir}${prefix}_${merger}_fastq_screen"
opts="--aligner bwa --threads ${cores} --subset ${fqssubset} --conf ${SCREENCONF} --force --outdir ${fqsdir}"
echo "Running FastQ-Screen on ${prefix} original reads ${read1} ${read2} and merging results ${mergedfq}, ${unmergedfq1}, and ${unmergedfq1} using configuration ${SCREENCONF}"
echo "Using BWA-MEM, ${cores} cores, --subset ${fqssubset} of reads (0 means all), and outputting to ${fqsdir}"
echo "Logs at logs/fastq_screen_${prefix}_${merger}.std{err,out}"
/usr/bin/time -v ${FASTQSCREEN} ${opts} ${read1} ${read2} ${mergedfq} ${unmergedfq1} ${unmergedfq2} 2> logs/fastq_screen_${prefix}_${merger}.stderr > logs/fastq_screen_${prefix}_${merger}.stdout

#Kraken2 against PlusPFP:
#kraken2 takes a TON of memory, more than 90 GB for the StandardPlusPFP database
#However, kraken2 is quite fast, taking a couple minutes for about 1M read pairs
#bracken is super fast and takes minimal memory

#Hard-code the options for the moment:
#Set a suffix for output files that indicates the database:
suffix="PlusPFP"

#Set the parameters for Bracken abundance estimation:
#Read length (150 bp):
readlen=150
#Taxonomic level for classification (S=species, G=genus, etc.):
taxlevel="S"
#Threshold for classification (minimum reads to allow down-propagation
# of reads from higher/broader taxa):
minreads=10

#Set the options, including multi-threading, database path,
# use taxon names, and decompress input with gzip:
k2opts="--threads ${cores} --db ${KRAKEN2DB} --use-names --gzip-compressed"
bopts="-d ${KRAKEN2DB} -r ${readlen} -l ${taxlevel} -t ${minreads}"

#Now run kraken2 and bracken:
#Merged:
k2report="${outputdir}${prefix}_${merger}_M_${suffix}.k2report"
echo "Kraken2 on ${merger} merged reads from sample ${prefix} against the Standard${suffix} database at ${KRAKEN2DB}:"
echo "kraken2 options: ${k2opts}"
/usr/bin/time -v ${KRAKEN} ${k2opts} --report ${k2report} ${mergedfq} | gzip -9 > ${outputdir}${prefix}_k2_${merger}_M_${suffix}.tsv.gz
echo "Bracken abundance estimation on ${merger} merged reads from sample ${prefix} against the Standard${suffix} database:"
echo "bracken options: ${bopts}"
/usr/bin/time -v ${BRACKEN} ${bopts} -i ${k2report} -o ${outputdir}${prefix}_k2_${merger}_M_${suffix}.bracken
#Unmerged:
k2report="${outputdir}${prefix}_${merger}_UM_${suffix}.k2report"
echo "Kraken2 on ${merger} unmerged reads from sample ${prefix} against the Standard${suffix} database at ${KRAKEN2DB}:"
echo "kraken2 options: ${k2opts} --paired"
/usr/bin/time -v ${KRAKEN} ${k2opts} --report ${k2report} --paired ${unmergedfq1} ${unmergedfq2} | gzip -9 > ${outputdir}${prefix}_k2_${merger}_UM_${suffix}.tsv.gz
echo "Bracken abundance estimation on ${merger} unmerged reads from sample ${prefix} against the Standard${suffix} database:"
echo "bracken options: ${bopts}"
/usr/bin/time -v ${BRACKEN} ${bopts} -i ${k2report} -o ${outputdir}${prefix}_k2_${merger}_UM_${suffix}.bracken

#Add the kraken2 and Bracken result summaries to the QC summary file:
#Note that the Bracken abundance is a percentage *of classified reads*,
# not *of all reads*.
#For now, we're hard-coding the taxid for Homo sapiens here:
taxid="9606"
#Kraken classified and unclassified counts, and classification percentage:
awk '$5=="0"{unclassified=$2;}$5=="1"{classified=$2;}END{printf "\t%i\t%i\t%f", classified, unclassified, classified*100/(classified+unclassified);}' ${outputdir}${prefix}_${merger}_M_${suffix}.k2report >> ${qc_summary}
#Bracken percent abundance for the designated taxid:
awk -v "taxid=${taxid}" 'BEGIN{FS="\t";}$2==taxid{printf "\t%f", 100*$7;}' ${outputdir}${prefix}_k2_${merger}_M_${suffix}.bracken >> ${qc_summary}
#Kraken classified and unclassified counts, and classification percentage:
awk '$5=="0"{unclassified=$2;}$5=="1"{classified=$2;}END{printf "\t%i\t%i\t%f", classified, unclassified, classified*100/(classified+unclassified);}' ${outputdir}${prefix}_${merger}_UM_${suffix}.k2report >> ${qc_summary}
#Bracken percent abundance for the designated taxid:
awk -v "taxid=${taxid}" 'BEGIN{FS="\t";}$2==taxid{printf "\t%f", 100*$7;}' ${outputdir}${prefix}_k2_${merger}_UM_${suffix}.bracken >> ${qc_summary}

#Mapping:
#BWA-MEM uses some memory, but not too much. Sorting takes more memory, but can
# be limited. The whole process is probably the slowest step in this pipeline
echo "Mapping ${prefix} merged reads to ${refid} with BWA-MEM and special options ${special}"
echo "Logs at logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_M.std{err,out}"
/usr/bin/time -v ${PSEUDOREF}/indexAlignDictMarkDup_v2.sh ${outputdir}${prefix}_${refid}_${merger}_M ${ref} ${mergedfq} ${cores} ${special} logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_M.stderr > logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_M.stdout
echo "Mapping ${prefix} unmerged reads to ${refid} with BWA-MEM and special options ${special}"
echo "Logs at logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_UM.std{err,out}"
/usr/bin/time -v ${PSEUDOREF}/indexAlignDictMarkDup_v2.sh ${outputdir}${prefix}_${refid}_${merger}_UM ${ref} ${unmergedfq1} ${unmergedfq2} ${cores} ${special} logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_UM.stderr > logs/PseudoreferencePipeline_MAP_${prefix}_${refid}_${merger}_UM.stdout
   #Params: PREFIX REF READS CORES SPECIAL
#   cmd="${PSEUDOREF}/indexAlignDictMarkDup_v2.sh ${prefix} ${ref} ${read1} ${read2} ${cores} ${special}"

#Indel Realignment if needed:
#echo "Realigning merged reads around indels for ${prefix} on ${refid} with GATK IndelRealigner and special options ${special}"
#echo "Logs at logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_M.std{err,out}"
#/usr/bin/time -v ${PSEUDOREF}/IndelRealignment.sh ${outputdir}${prefix}_${refid}_${merger}_M ${ref} ${special} 2> logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_M.stderr > logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_M.stdout
#echo "Realigning unmerged reads around indels for ${prefix} on ${refid} with GATK IndelRealigner and special options ${special}"
#echo "Logs at logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_UM.std{err,out}"
#/usr/bin/time -v ${PSEUDOREF}/IndelRealignment.sh ${outputdir}${prefix}_${refid}_${merger}_UM ${ref} ${special} 2> logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_UM.stderr > logs/PseudoreferencePipeline_IR_${prefix}_${refid}_${merger}_UM.stdout

   #Params: PREFIX REF SPECIAL
#   cmd="${PSEUDOREF}/IndelRealignment.sh ${prefix} ${ref} ${special}"

#$cmd

#Mapping rates (merged, unmerged):
#Calculating mapping rates is pretty fast and takes very little memory
echo "Calculating mapping rates of ${prefix} to reference ${refid} for merged and unmerged reads based on merger ${merger}"
${SAMTOOLS} view -F0xB00 ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | \
   ${scriptdir}/mapping_rate.awk -v "indiv=${prefix}" -v "ref=${refid}" -v "merger=${merger}" -v "source=merged" > ${outputdir}${prefix}_${refid}_${merger}_mapping_rates.tsv
${SAMTOOLS} view -F0xB00 ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | \
   ${scriptdir}/mapping_rate.awk -v "indiv=${prefix}" -v "ref=${refid}" -v "merger=${merger}" -v "source=unmerged" -v "noheader=1" >> ${outputdir}${prefix}_${refid}_${merger}_mapping_rates.tsv

#Add the mapping rate summaries to the QC summary file:
awk 'BEGIN{FS="\t";}$4=="M"{printf "\t%s", $7;}$4=="UM"{printf "\t%s\t%s", $7, $9;}' ${outputdir}${prefix}_${refid}_${merger}_mapping_rates.tsv >> ${qc_summary}

#Insert size distributions (merged+mapped, merged+unmapped, and unmerged+mapped):
#Insert size distributions are slow to generate, but don't require too much
# memory, since we use GNU coreutils sort
echo "Computing insert size distributions for ${prefix} mapped to reference ${refid} based on merger ${merger}"
#Use read length stored in column 10 for insert size of merged reads
#Exclude unmapped reads for merged+mapped (-F0x4)
#Exclude secondary and supplementary alignments and vendor QC fails (-F0xB00)
${SAMTOOLS} view -F0xB04 ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | \
   awk 'BEGIN{FS="\t";OFS=FS;}{print length($10);}' | \
   sort -k1,1n | \
   uniq -c | \
   sort -k2,2n | \
   awk -v "indiv=${prefix}" -v "ref=${refid}" -v "merger=${merger}" -v "datatype=merged" -v "state=mapped" 'BEGIN{OFS="\t";}{print indiv, ref, merger, datatype, state, $1, $2;}' > ${outputdir}${prefix}_UMproperlypaired_read_length_ePMFs.tsv
#Use read length as above
#Include unmapped reads for merged+unmapped (-f0x4)
#Exclude secondary and supplementary alignments and vendor QC fails (-F0xB00)
${SAMTOOLS} view -f0x4 -F0xB00 ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | \
   awk 'BEGIN{FS="\t";OFS=FS;}{print length($10);}' | \
   sort -k1,1n | \
   uniq -c | \
   sort -k2,2n | \
   awk -v "indiv=${prefix}" -v "ref=${refid}" -v "merger=${merger}" -v "datatype=merged" -v "state=unmapped" 'BEGIN{OFS="\t";}{print indiv, ref, merger, datatype, state, $1, $2;}' >> ${outputdir}${prefix}_UMproperlypaired_read_length_ePMFs.tsv
#Use template length as insert size (column 9), but only include positive
# lengths, since every pair will have two template lengths
#Include properly paired reads for unmerged+mapped (-f0x2)
#Exclude unmapped reads for unmerged+mapped (-F0x4)
#Exclude secondary and supplementary alignments and vendor QC fails (-F0xB00)
${SAMTOOLS} view -F0xB04 -f0x2 ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | \
   awk 'BEGIN{FS="\t";OFS=FS;}$9>=0{print $9;}' | \
   sort -k1,1n | \
   uniq -c | \
   sort -k2,2n | \
   awk -v "indiv=${prefix}" -v "ref=${refid}" -v "merger=${merger}" -v "datatype=unmerged" -v "state=mapped" 'BEGIN{OFS="\t";}{print indiv, ref, merger, datatype, state, $1, $2;}' >> ${outputdir}${prefix}_UMproperlypaired_read_length_ePMFs.tsv

if [[ "${refid}" == "hs37d5" || "${refid}" == "GRCh38" ]]; then
   #mtDNA haplogroup determination:
   #mutserve is pretty fast, haplogrep is faster, both take minimal memory
   opts="--reference ${RCRS} --threads ${cores}"
   echo "Running mutserve to call mtDNA variants from ${merger} merged reads for sample ${prefix} mapped to reference ${refid}"
   /usr/bin/time -v ${MUTSERVE} call ${opts} --output ${outputdir}${prefix}_${refid}_${merger}_M_noMD_mutserve.vcf.gz ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam 2> logs/mutserve_${prefix}_${refid}_${merger}_M_noMD.stderr > logs/mutserve_${prefix}_${refid}_${merger}_M_noMD.stdout
   echo "Running mutserve to call mtDNA variants from ${merger} unmerged reads for sample ${prefix} mapped to reference ${refid}"
   /usr/bin/time -v ${MUTSERVE} call ${opts} --output ${outputdir}${prefix}_${refid}_${merger}_UM_noMD_mutserve.vcf.gz ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam 2> logs/mutserve_${prefix}_${refid}_${merger}_UM_noMD.stderr > logs/mutserve_${prefix}_${refid}_${merger}_UM_noMD.stdout

   opts="--format vcf"
   echo "Running haplogrep to assign mtDNA haplogroup for ${merger} merged reads for sample ${prefix} mapped to reference ${refid}"
   /usr/bin/time -v ${HAPLOGREP} classify ${opts} --in ${outputdir}${prefix}_${refid}_${merger}_M_noMD_mutserve.vcf.gz --out ${outputdir}${prefix}_${refid}_${merger}_M_noMD_haplogrep.tsv 2> logs/haplogrep_${prefix}_${refid}_${merger}_M_noMD.stderr > logs/haplogrep_${prefix}_${refid}_${merger}_M_noMD.stdout
   echo "Running haplogrep to assign mtDNA haplogroup for ${merger} unmerged reads for sample ${prefix} mapped to reference ${refid}"
   /usr/bin/time -v ${HAPLOGREP} classify ${opts} --in ${outputdir}${prefix}_${refid}_${merger}_UM_noMD_mutserve.vcf.gz --out ${outputdir}${prefix}_${refid}_${merger}_UM_noMD_haplogrep.tsv 2> logs/haplogrep_${prefix}_${refid}_${merger}_UM_noMD.stderr > logs/haplogrep_${prefix}_${refid}_${merger}_UM_noMD.stdout
   #Output the mtDNA haplogroup to the QC summary file:
   awk 'BEGIN{FS="\t";}NR==2{gsub(/"/, ""); printf "\t%s\t%s", $2, $4;}' ${outputdir}${prefix}_${refid}_${merger}_M_noMD_haplogrep.tsv >> ${qc_summary}
   awk 'BEGIN{FS="\t";}NR==2{gsub(/"/, ""); printf "\t%s\t%s", $2, $4;}' ${outputdir}${prefix}_${refid}_${merger}_UM_noMD_haplogrep.tsv >> ${qc_summary}
else
   #Run mtArchitect?
   echo "Running mtArchitect not yet implemented"
fi

#Calculate some average read depths for merged and unmerged reads:
#These depth calculations are reasonably fast, especially due to using -r,
# though their runtime is dependent on chromosome length. Memory is basically
# constant and low
#First is mtDNA depth:
#printf "mtDNA (${mt}) average depth (merged) = "
${SAMTOOLS} depth -aa -r "${mt}" ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{printf "\t%f", sum/count;}' >> ${qc_summary}
#printf "mtDNA (${mt}) average depth (unmerged) = "
${SAMTOOLS} depth -aa -r "${mt}" ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{printf "\t%f", sum/count;}' >> ${qc_summary}

if [[ ! -z "${y}" ]]; then
   #Now do Y depth:
#   printf "Y (${y}) average depth (merged) = "
   ymdp=`${SAMTOOLS} depth -aa -r "${y}" ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`
#   printf "Y (${y}) average depth (unmerged) = "
   yumdp=`${SAMTOOLS} depth -aa -r "${y}" ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`

   #Compare to chr19 depth:
#   printf "Autosome with Y-like length (${autoy}) average depth (merged) = "
   autoymdp=`${SAMTOOLS} depth -aa -r "${autoy}" ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`
#   printf "Autosome with Y-like length (${autoy}) average depth (unmerged) = "
   autoyumdp=`${SAMTOOLS} depth -aa -r "${autoy}" ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`

   #Output depths and ratios to QC summary file:
   yratiom=$(awk -v "ymdp=${ymdp}" -v "autoymdp=${autoymdp}" 'BEGIN{print ymdp/autoymdp;}')
   yratioum=$(awk -v "yumdp=${yumdp}" -v "autoyumdp=${autoyumdp}" 'BEGIN{print yumdp/autoyumdp;}')
   printf "\t${ymdp}\t${autoymdp}\t${yratiom}\t${yumdp}\t${autoyumdp}\t${yratioum}" >> ${qc_summary}

   #Use Serena's method based on fraction of total reads mapping to Y
   # > about 0.2% => male to estimate sex:
   ymapthresh="0.002"
   ${SAMTOOLS} view -F0xB04 ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk -v "y=${y}" -v "thresh=${ymapthresh}" 'BEGIN{FS="\t";}!/^@/{count[$3]++;total++;}END{if (count[y]/total > thresh) {sex="M";} else {sex="F";}; printf "\t%i\t%i\t%f\t%s", count[y], total, count[y]/total, sex;}' >> ${qc_summary}
   ${SAMTOOLS} view -F0xB04 ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk -v "y=${y}" -v "thresh=${ymapthresh}" 'BEGIN{FS="\t";}!/^@/{count[$3]++;total++;}END{if (count[y]/total > thresh) {sex="M";} else {sex="F";}; printf "\t%i\t%i\t%f\t%s", count[y], total, count[y]/total, sex;}' >> ${qc_summary}
fi

if [[ ! -z "${x}" ]]; then
   #Now do X depth:
#   printf "X (${x}) average depth (merged) = "
   xmdp=`${SAMTOOLS} depth -aa -r "${x}" ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`
#   printf "X (${x}) average depth (unmerged) = "
   xumdp=`${SAMTOOLS} depth -aa -r "${x}" ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`

   #Compare to chr7 depth:
#   printf "Autosome with X-like length (${autox}) average depth (merged) = "
   autoxmdp=`${SAMTOOLS} depth -aa -r "${autox}" ${outputdir}${prefix}_${refid}_${merger}_M_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`
#   printf "Autosome with X-like length (${autox}) average depth (unmerged) = "
   autoxumdp=`${SAMTOOLS} depth -aa -r "${autox}" ${outputdir}${prefix}_${refid}_${merger}_UM_sorted.bam | awk 'BEGIN{FS="\t";}{sum+=$3;count+=1;}END{print sum/count;}'`

   #Output depths and ratios to QC summary file:
   xratiom=$(awk -v "xmdp=${xmdp}" -v "autoxmdp=${autoxmdp}" 'BEGIN{print xmdp/autoxmdp;}')
   xratioum=$(awk -v "xumdp=${xumdp}" -v "autoxumdp=${autoxumdp}" 'BEGIN{print xumdp/autoxumdp;}')
   printf "\t${xmdp}\t${autoxmdp}\t${xratiom}\t${xumdp}\t${autoxumdp}\t${xratioum}" >> ${qc_summary}
   #Estimate sex from x to autosome ratio > about 0.7 => female to estimate sex:
   xratiothresh="0.7"
   awk -v "xratio=${xratiom}" -v "thresh=${xratiothresh}" 'BEGIN{if (xratio > thresh) {sex="F";} else {sex="M";}; printf "\t%s", sex;}' >> ${qc_summary}
   awk -v "xratio=${xratioum}" -v "thresh=${xratiothresh}" 'BEGIN{if (xratio > thresh) {sex="F";} else {sex="M";}; printf "\t%s", sex;}' >> ${qc_summary}
#   if [[ "${xratiom}" -gt "${xratiothresh}" ]]; then
#      printf "\tF" >> ${qc_summary}
#   else
#      printf "\tM" >> ${qc_summary}
#   fi
#   if [[ "${xratioum}" -gt "${xratiothresh}" ]]; then
#      printf "\tF" >> ${qc_summary}
#   else
#      printf "\tM" >> ${qc_summary}
#   fi
fi

echo "The QC summary file can be found at: ${qc_summary}"
