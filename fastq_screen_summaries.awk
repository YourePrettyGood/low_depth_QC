#!/bin/awk -f
#This script reformats the *_screen.txt outputs from FastQ-Screen.
#We only output the numbers of reads and ignore all percentages,
# as these can be recomputed on-the-fly. We also add sample ID
# and read file info to each line.
#Output is in "tall" tabular format, with one line per statistic
# per input file.
#The original purpose is to concatenate QC summaries from multiple samples,
# but retain a header for easy import into R.
BEGIN{
   FS="\t";
   OFS=FS;
   #Print a header:
   print "SampleID", "ReadFile", "Bait", "Total", "Unmapped", "OneHitOneGenome", "ManyHitsOneGenome", "OneHitManyGenomes", "ManyHitsManyGenomes";
}
FNR==1{
   n_pathelems=split(FILENAME, fnpath, "/");
   n_fnelems=split(fnpath[n_pathelems], fnparts, "_");
   fn=fnpath[n_pathelems];
   id=fnparts[1];
}
!/^#/&&FNR>2&&NF>1{
   print id, fn, $1, $2, $3, $5, $7, $9, $11;
}
