#!/bin/bash

ref=$1

if [[ "${ref}" == "hs37d5" ]]; then
   refpath="/gpfs/gibbs/pi/tucci/pfr8/refs/1kGP/hs37d5/hs37d5.fa";
elif [[ "${ref}" == "GRCh38" ]]; then
   refpath="/gpfs/gibbs/pi/tucci/pfr8/refs/1kGP/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa";
elif [[ "${ref}" == "Kamilah" ]]; then
   refpath="/gpfs/gibbs/pi/tucci/pfr8/refs/Gorilla/Kamilah_GGO_v0/Kamilah_GGO_v0.fasta";
else
   echo "Reference ${ref} not recognized, quitting." > "/dev/stderr";
   exit 2;
fi

while read line;
   do
   prefix="${line}/${line}"
   r1=(${prefix}_*_R1*fastq.gz);
   r2=(${prefix}_*_R2*fastq.gz);
   printf "${prefix}\t${refpath}\t${r1}\t${r2}\n";
done
