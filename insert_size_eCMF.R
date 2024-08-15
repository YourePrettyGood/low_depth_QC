#!/usr/bin/env Rscript
#Read in arguments from the command line:
options <- commandArgs(trailingOnly=TRUE)

#First argument on the command line should be the individual's ID:
ID <- options[1]

#Load tidyverse for manipulating dataframes and plotting:
library(tidyverse)

#Set some variables for plotting:
width_cm <- 16.0 #Plot width in cm
height_cm <- 12.0 #Plot height in cm
dpi <- 500 #Resolution of the plot in dots per inch
line_size <- 0.5 #Line thickness
point_size <- 0.5 #Point radius
alpha <- 0.5 #Make the lines and points somewhat transparent
legend_pos <- c(0.85, 0.8) #Legend position, x and y from bottom left
max_lens <- c(1000) #Upper limits on insert length, one per plot
#The limits used in the report from Claret Bio/Ed Green would be c(400, 1000)
#The plots look pretty clean with just one maximum at 1000
#There are occasional gigantic TLENs which we don't want to plot, as they're
# either mapping or reference artifacts.
max_count <- 1000000 #This makes the plots consistent with each other
#The limit used in the report from Claret Bio/Ed Green is pretty reasonable
# for around < 10% mapping for ~1M total read pairs, but we're dealing with
# more like 100% mapping for ~20-30M total read pairs

#Load the insert size distribution TSV into a dataframe:
insert_lengths <- read.table(paste0(ID, "/",
                                    ID, "_UMproperlypaired_read_length_ePMFs.tsv"),
                             colClasses=c(rep("character", 5),
                                          rep("numeric", 2)),
                             col.names=c("Individual", "Reference", "Merger",
                                         "DataType", "State", "Count", "Length"),
                             header=FALSE,
                             sep="\t")

#Check the input to make sure we're only dealing with one reference and
# merger at a time:
stopifnot(length(unique(insert_lengths$Reference)) == 1,
          length(unique(insert_lengths$Merger)) == 1)

#Now we get to the plots:
for (max_len in max_lens) {
   for (indiv in unique(insert_lengths$Individual)) {
      for (ref in unique(insert_lengths$Reference)) {
         for (merger in unique(insert_lengths$Merger)) {
            insert_lengths %>% unite("Pipeline", DataType, State, sep=", ") %>%
                               filter(Individual == indiv,
                                      Reference == ref,
                                      Merger == merger,
                                      Length >= 0,
                                      Length <= max_len) %>%
                               ggplot(aes(x=Length,
                                          y=Count,
                                          colour=Pipeline,
                                          fill=Pipeline)) +
                                  geom_line(size=line_size, alpha=alpha) +
                                  geom_point(size=point_size, alpha=alpha) +
                                  scale_y_log10(limits=c(1, max_count)) +
                                  annotation_logticks(sides="l") +
                                  scale_x_continuous(limits=c(0, max_len)) +
                                  scale_colour_brewer(palette="Set2") +
                                  scale_fill_brewer(palette="Set2") +
                                  labs(x="Insert Length (bp)",
                                       y="Read (Pair) Count",
                                       title=paste("Insert length distribution for individual", indiv),
                                       subtitle=paste("on", ref, "using", merger)) +
                                  theme_bw() +
                                  theme(legend.position=legend_pos,
                                        legend.direction="vertical",
                                        legend.title=element_blank(),
                                        axis.title=element_text(face="bold"),
                                        plot.title=element_text(hjust=0.5),
                                        plot.subtitle=element_text(hjust=0.5))
         }
      }
      ggsave(paste0(indiv, "/",
             paste(indiv, ref, merger,
                   "max_len", max_len, "insert_lengths.pdf",
                   sep="_")),
             width=width_cm,
             height=height_cm,
             units="cm",
             dpi=dpi)
   }
}
