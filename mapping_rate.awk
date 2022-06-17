#!/bin/awk -f
#This script partitions primary alignments and calculates the
# mapping rate, rate of proper pairs, mapped orphans (i.e.
# single-ended reads in the input), and singletons (i.e.
# reads from a pair where the other read is unmapped).
#The input SAM is expected to consist of only primary alignments,
# so it would typically be produced with:
# samtools view -F0xB00
# (0x100 is secondary, 0x200 is fail vendor QC, and 0x800 is supplementary)
#Input arguments:
# indiv:    ID of individual/sample (e.g. NWM212)
# ref:      ID for reference genome (e.g. hs37d5)
# merger:   Name of program used for PE merging (e.g. SeqPrep or BBMerge)
# source:   Type of input data (e.g. unmerged or merged)
#SAM can be piped into STDIN or passed as first positional argument
#Note that indiv, ref, merger, and source are passed straight through
# to the output file.
#Optional argument:
# noheader: Do not output a header line describing the columns
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(indiv) == 0) {
      print "Missing indiv argument, quitting." > "/dev/stderr";
      exit 2;
   };
   if (length(ref) == 0) {
      print "Missing ref argument, quitting." > "/dev/stderr";
      exit 3;
   };
   if (length(merger) == 0) {
      print "Missing merger argument, quitting." > "/dev/stderr";
      exit 4;
   };
   if (length(source) == 0) {
      print "Missing source argument, quitting." > "/dev/stderr";
      exit 5;
   };
   if (length(noheader) == 0) {
      print "Indiv", "Ref", "Merger", "Source", "NumMapped", "NumPrimaryAlns", "PctMapped", "NumProperlyPairedReads", "PctProperlyPaired", "NumSE", "PctSE", "NumSingletons", "PctSingletons";
   };
   #Primary alignments are any alignments in the SAM that do not have
   # the secondary (0x100), QC/vendor fail (0x200), or supplementary (0x800)
   # flags set, so they are the longest aligned segment of a given read
   # mapping to the best location. Primary alignments can also be unmapped,
   # in which case they correspond to the entire read.
   primaries=0;
   #Mapped here corresponds to mapped primary alignments, so only excluding
   # unmapped primary alignments. Unmapped alignments are indicated by the
   # 0x4 flag.
   mapped=0;
   #Proper here counts only once for each read pair that maps concordantly
   # within the expected size distribution. Concordant mapping means that
   # the reads point toward each other to form a fragment. Proper pairing
   # is indicated by the 0x2 flag, which only technically has meaning if
   # 0x1 is set. We're not meeting spec here, unfortunately.
   proper=0;
   #A singleton is a read from a pair that is mapped, but its mate is
   # unmapped. Thus the 0x1 flag is set, 0x4 is unset, and 0x8 is set.
   # The 0x8 flag means the mate is unmapped.
   singletons=0;
   #Here we use "orphan" to entail an unpaired read, i.e. single end. We use
   # it as a catch-all for processing merged and unmerged datasets, but it
   # would also capture cases in BAMs where one end of a pair is dropped
   # prior to mapping.
   orphans=0;
}
!/^@/{
   if (and($2, 0x004) == 0) {
      mapped+=1;
      if (and($2, 0x002) > 0) {
         proper+=1;
      } else if (and($2, 0x001) > 0) {
         if (and($2, 0x008) > 0) {
            singletons+=1;
         };
      } else {
         orphans+=1;
      };
   };
   primaries+=1;
}
END{
   if (primaries > 0) {
      print indiv, ref, merger, source, mapped, primaries, mapped*100/primaries, proper, proper*100/primaries, orphans, orphans*100/primaries, singletons, singletons*100/primaries;
   };
}
