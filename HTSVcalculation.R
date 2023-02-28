# R script to get information of mapping reads, to add virus information and to plot viral reads
# Publication: Bioinformatics pipeline for high-throughput sequencing detection of plant viruses
# These steps are optional in our pipeline
# Author: Livia Donaire (ldonaire@abiopep.com/ldonaire@cebas.csic.es)

## Install required R packages
# Install dplyr
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")

# Install ggplot2 and purrs
# The easiest way to get ggplot2 and purrr is to install the whole tidyverse:
install.packages("tidyverse")

## Load required R packages
library("dplyr")
library("ggplot2")
library("purrr")

## Set working directory ##
path <-"/mnt/Linux_data/MiMB_protocol" #Please, write here the complete part of your working directory
setwd(path)

# Required files (obtained in previous steps in the pipeline):
# Raw fastq.gz files in a folder called "Original_files"
# Filtered fastq.gz files in a folder called "trimming"
# Optional: Clean fastq.gz files in a folder called "bowtie2_plant" 
# sort.bam files of read alignments against viral genomes saved in a folder called "bwa"

############################# Part one ################################################
############################# Counting reads ##################################
## Code to count number of raw, filtered, cleaned reads and viral reads and to export results to a tab-delimited file ##

# Create variables to save results
raw <- integer()
clean <- integer()
filter<- integer()

# List files to count
myfileRaw <- list.files("Original_files", pattern=".fastq.gz$")
myfileClean <- list.files("trimming", pattern=".fastq.gz$") 
myfileFilt <- list.files("bowtie2_plant", pattern=".fastq.gz$") 

## Count reads in raw FASTQ files
for (i in 1:length(myfileRaw)){
  raw[i] <- as.numeric(system(sprintf("echo $(zcat Original_files/%s|wc -l)/4|bc", myfileRaw[i]), intern=TRUE)) 
}
# Sum PE counts
rawT <- data.frame("Raw_reads"=sapply(seq(1,length(raw),by=2),function(i) sum(raw[c(i,i+1)],na.rm=FALSE)))

## Count reads in cleaned FASTQ files
for (i in 1:length(myfileClean)){
  clean[i] <- as.numeric(system(sprintf("echo $(zcat trimming/%s|wc -l)/4|bc", myfileClean[i]), intern=TRUE)) 
}
# sum PE and single counts
cleanT <- data.frame("Cleaned_reads"=sapply(seq(1,length(clean),by=4),function(i) sum(clean[c(i,i+1,i+2,i+3)],na.rm=FALSE)))
# add percentages
cleanT$Percentage <- round((cleanT$Cleaned_reads/rawT$Raw_reads)*100,2)

## Count reads in filtered FASTQ files (optional)
for (i in 1:length(myfileFilt)){
  filter[i] <- as.numeric(system(sprintf("echo $(zcat bowtie2_plant/%s|wc -l)/4|bc", myfileFilt[i]), intern=TRUE)) 
}
# sum PE and single counts
filterT <- data.frame("Filtered_reads"=sapply(seq(1,length(filter),by=3),function(i) sum(filter[c(i,i+1,i+2)],na.rm=FALSE)))
# add percentages
filterT$Percentage <- round((filterT$Filtered_reads/cleanT$Cleaned_reads)*100,2)

## Count reads mapping to reference viruses using SeqKit
system("cat sample.names.txt | parallel \"seqkit bam -c bwa/{}_virus_reads.txt bwa/{}_bwaVirusFound.sort.bam\"")
# Load files in R
myCountNames<-paste("COUNT", mySamples, sep="_")
myCountFiles <- list.files(path="bwa", pattern="*virus_reads.txt", full.names = TRUE)
for (i in 1:length(myCountNames)) assign(myCountNames[i], read.csv(myCountFiles[i], sep="\t", header=TRUE, stringsAsFactors = FALSE,
                                                                   col.names=c("Accession","Counts","Second_counts","Super_counts")))
# Load all tables in lists
myCountList <-  mget(myCountNames)
myCountList <- lapply(myCountList, function(x) x <- x[rowSums(x[,c(2:4)])>0,]) # remove not mapping data
myCountList <- lapply(myCountList, function(x) x[,c(1:2)]) # select only primary counts
# Get total number of viral reads for each sample in a data frame
virusT <- data.frame("Viral_reads"=unlist(lapply(myCountList, function(x) x <- colSums(x["Counts"]))))
# add percentages
virusT$Percentage <- round((virusT$Viral_reads/filterT$Filtered_reads)*100,2)

## Join all in a table
myCountRes<- cbind(rawT,cleanT,filterT,virusT)
# in case no filtered reads are available use:
# myCountRes<- cbind(rawT,cleanT,virusT)
rownames(myCountRes) <- unique(sub('\\_1.fastq.gz|_2.fastq.gz',"",myfileRaw))

## Export the table to a file
write.csv(myCountRes,"Summary_of_mapping_reads.csv")

############################# Part two ##################################################
############################# Calculating basic parameters on viral reads ###############
############################# Add virus information #####################################

## Get sequencing depth using SAMtools depth from the BAM files generated before
system("cat sample.names.txt | parallel \"samtools depth bwa/{}_bwaVirusFound.sort.bam > bwa/{}_depth.txt\"")

## Get the length of the viruses using SAMtools faidx
system("samtools faidx virusDB/VirusFound.fasta")
myVirusLengths <- read.csv("virusDB/VirusFound.fasta.fai", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# Lengths are column 2
myVirusLengths <- myVirusLengths[,1:2]
colnames(myVirusLengths) <- c("Accession","Virus_length")

## Load files in R and make some calculations
myDepthNames<-paste("DEPTH", mySamples, sep="_")
myDepthFiles <- list.files(path="bwa", pattern="*depth.txt", full.names = TRUE)
for (i in 1:length(myDepthNames)) assign(myDepthNames[i], read.csv(myDepthFiles[i], sep="\t", header=FALSE, stringsAsFactors = FALSE,
                                                                   col.names=c("Accession","Position","Depth")))

# Load all tables in lists
myDepthList <-  mget(myDepthNames)

# Create a function to make the calculations (don't forget to load dplyr package)
myCoverage<-function(x,myVirusLenghts){ # where y is myDEPTH_list
  myMeandDepth <- aggregate(cbind(x$Depth)~x$Accession, x, mean)
  colnames(myMeandDepth) <- c("Accession","Average_depth")
  myPosCov <- aggregate(cbind(x$Depth)~x$Accession, x, length)
  colnames(myPosCov) <- c("Accession","Length_covered")
  myStats <- full_join(myMeandDepth,myPosCov,"Accession") %>% full_join(.,myVirusLengths,"Accession")
  myStats <- transform(myStats,Perc_covered=(myStats$Length_covered/myStats$Virus_length)*100)
  myStats <- return(myStats)
}

# Apply the function myCoverage to each item in myDepthNames
mySTAT <- list()
for (i in 1:length(myDepthNames)) {
  mySTAT <- lapply(myDepthList, myCoverage)
}

# remove rows with empty cells (in case there is no data for an specific virus in a sample)
mySTAT <- lapply(mySTAT, function(x) na.omit(x))

## add virus information 
# load the csv file downloaded from virus NCBI portal (instructions in the pipeline paper)
myVirusInfo <- read.csv("virusDB/pvDB.csv")

# create a list to merge with mySTAT list
myVirusInfoList <- lapply(1:length(mySamples), function(x) myVirusInfo)

# merge mySTAT with myCountList (created in part 1) (don't forget to load purrr package)
mySTATtmp <- map2(mySTAT , myCountList, left_join , by="Accession")
# merge with myVirusInfoList
mySTATfin <- map2(mySTATtmp , myVirusInfoList , left_join , by="Accession")
names(mySTATfin) <- mySamples
# save results into individual csv files
lapply(names(mySTATfin), function(x) {
  f <- mySTATfin[[x]]
  write.csv(f, file=paste0(x,"_final_results.csv"),row.names=FALSE)
}) 

############################# Part three ##################################################
############################# Create plots of read coverage ###############################

# Create a function to make the same plot for different samples (don't forget to load dplyr and ggplot2 packages)
myPlotFunc <- function(y) {   # where y is myDEPTH_list
  y <- left_join(y,myVirusInfo,"Accession") %>% left_join(.,myVirusLengths, "Accession") # Join myVirusInfo with myVirusLengths (created previously)
  y$Virus <- paste(y$Species,y$Segment,y$Accession,paste0("(",y$Virus_length," bp)"))   # Join Species, segment and accession in a new column called Virus
  y$Virus <- as.factor(y$Virus)           # transform Virus column to factor
  y$Virus <- factor(y$Virus,levels = rev(levels(y$Virus))) # order the levels
  Mycolor   <- "blue3"                # color for plotting reads
  ggplot(y, aes(x=Virus, y=Position))+    # plot Virus Vs Position
    geom_boxplot() + 
    geom_jitter(shape = 20,
                position = position_jitter(0.25),
                size = 1.5,
                colour = Mycolor) +
    coord_flip() +
    theme_bw() +
    labs(x = "", y = "Genome position") +
    theme(plot.title = element_text(size = 12),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 8, face = "bold.italic"),
          axis.text.x = element_text(size = 8, face = "bold"))
}

# Create the plots in a List
plots <- lapply(myDepthList, myPlotFunc)
names(plots) <- mySamples

# Add a title to each plot according to mySamples
for (i in seq_along(mySamples)) {
  plots[[i]] <- plots[[i]] + ggtitle(mySamples[i])
}

# Use print(plots[[1]]) to print individual plots in the Plots tab (right bottom panel) (change 1 by successive numbers according to the number of samples)

# Export plots in a pdf file per each sample
lapply(names(plots), function(x) {
  pdf(file=paste0(x,"_coverage_plot.pdf"),onefile=FALSE) # open graphical device
  print(plots[[x]])
})
graphics.off() # close all graphical devices

########################### END ###############################################################