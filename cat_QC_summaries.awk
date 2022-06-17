#!/bin/awk -f
#This script is super simplistic, but it basically only retains the first
# line from the first file supplied, and retains all other lines from
# all files. It's a simple way to only retain the header line of one
# file from multiple and concatenate them together.
#The original purpose is to concatenate QC summaries from multiple samples,
# but retain a header for easy import into R.
FNR==NR&&FNR==1{
   print;
}
FNR>1{
   print;
}
