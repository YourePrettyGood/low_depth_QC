

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
max_lens <- c(400, 1000) #Upper limits on insert length, one per plot
#We use the same limits as the report from Claret Bio/Ed Green
#The plots look pretty clean with these two maxima
#There are occasional gigantic TLENs which we don't want to plot, as they're
# either mapping or reference artifacts.
max_count <- 10000 #This makes the plots consistent with each other
#We use the same limit as the report from Claret Bio/Ed Green
#This is pretty reasonable for around < 10% mapping for ~1M total read pairs

#Load the insert size distribution TSV into a dataframe:
insert_lengths <- read.table(paste(),
                             colClasses=c(),
                             col.names=c(),
                             header=FALSE)


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

ggsave(paste0(, paste(indiv, ref, merger, "max_len", max_len, "insert_lengths.pdf", sep="_")), width=width_cm, height=height_cm, dpi=dpi)
