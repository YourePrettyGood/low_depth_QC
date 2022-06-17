#!/bin/bash

if [[ "$#" -lt "4" ]]; then
   echo "prep_sbatch.sh usage:"
   echo "prep_sbatch.sh [Reference ID] [Batch ID] [Batch size] [Yale e-mail] [Special options]"
   echo "Required arguments:"
   echo " Reference ID can be hs37d5, GRCh38, or Kamilah"
   echo "  hs37d5 is the previous version of the Human reference genome"
   echo "  GRCh38 is the latest version of the Human reference genome"
   echo "  Kamilah is the latest version of the Gorilla reference genome"
   echo " Batch ID can be any short string"
   echo "  For instance, Pilot6, or Baining13"
   echo " Batch size must be a positive integer that is the number of samples"
   echo "  in the current batch."
   echo "  For instance, for the Baining13 batch it was 13"
   echo " Yale e-mail is your full Yale e-mail address"
   echo "Optional arguments:"
   echo " Special options can be left blank/omitted, but is meant for finer"
   echo "  control over the QC pipeline."
   exit 1;
fi

#Determine the short ID for the reference (used in the SLURM job name):
ref=$1
if [[ "${ref}" == "hs37d5" ]]; then
   shortref="hg19";
elif [[ "${ref}" == "GRCh38" ]]; then
   shortref="hg38";
elif [[ "${ref}" == "Kamilah" ]]; then
   shortref="GGO";
else
   echo "Reference ${ref} not recognized, quitting." > "/dev/stderr";
   exit 2;
fi

batch=$2

#Validate that batchsize is a positive number:
batchsize=$3
if [[ "${batchsize}" -lt "1" ]]; then
   echo "Invalid batch size less than 1: ${batchsize}. Quitting." > "/dev/stderr";
   exit 3;
fi

#Validate the e-mail address:
email=$4
if [[ ! "${email}" =~ ^[A-Za-z0-9._+-]+[@]yale\.edu$ ]]; then
   echo "Invalid Yale e-mail specified: ${email}. Quitting." > "/dev/stderr";
   exit 4;
fi

#No validation on special options:
special=$5

#Determine the directory of the prep_sbatch.sh script so we have the proper
# path for the array_template.sbatch file no matter where prep_sbatch.sh is
# called:
scriptdir=`dirname $0`
sed -e "s/\[BATCH\]/${batch}/g" -e "s/\[BATCHSIZE\]/${batchsize}/g" -e "s/\[EMAIL\]/${email}/g" -e "s/\[REF\]/${ref}/g" -e "s/\[SHORTREF\]/${shortref}/g" -e "s/\[SPECIAL\]/${special}/g" < ${scriptdir}/array_template.sbatch
