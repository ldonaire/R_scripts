# R script for virus discovery using HTS 
# Publication: Bioinformatics pipeline for high-throughput sequencing detection of plant viruses
# Contains commands to run the full pipeline
# Author: Livia Donaire (ldonaire@abiopep.com/ldonaire@cebas.csic.es)


#######################################################################################################
########## Note: full path to different programs must be specified if they are not in PATH ############
#######################################################################################################

## Install required R packages
# Install dplyr
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")

# Install ggplot2 and purrr
install.packages("tidyverse")

## Load required R packages
library("dplyr")
library("ggplot2")
library("purrr")

## Set working directory ##
path <-"/mnt/Linux_data/MiMB_protocol" #Please, write here the complete path of your working directory
setwd(path)

## Build a nucleotide BLAST database ####
dir.create("virusDB")
# Instructions to download "pvDB.fasta" in the protocol paper
system("makeblastdb -in virusDB/pvDB.fasta -dbtype nucl -parse_seqids -out virusDB/pvDB")

## Build a plant host database to subtract plant reads ####
# Download Tomato Heinz1706 genomic sequences from SGN (ITAG2.4 release) (accessed Feb 2023)
dir.create("plantDB")
MyUrl="https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/Heinz1706/annotation/ITAG2.4_release/ITAG2.4_genomic.fasta"
download.file(MyUrl, "plantDB/ITAG2.4_genomic.fasta")
# Create a Bowtie 2 index
system("bowtie2-build plantDB/ITAG2.4_genomic.fasta plantDB/ITAG2.4_genomic")

## Get example data ####
dir.create("Original_files")
system("fasterq-dump SRR14066947 SRR14066948  --outdir Original_files")

## Create an R variable and an external file containing sample names ####
mySamples <- list.files("Original_files", pattern=".fastq.gz$")
mySamples<-unique(sub('\\_1.fastq.gz|_2.fastq.gz','', mySamples)) 
write(mySamples,"sample.names.txt")

## QC of raw reads ####
dir.create("FastQC_raw")
system("fastqc -t 4 Original_files/*.fastq.gz -o FastQC_raw")
# -t n (number of samples analyzed simultaneously)
# -o string (folder to save the results)

## Trimming of raw reads ####
dir.create("trimming")
# use Trimmomatic v0.39 with parallel for more than one sample
system("cat sample.names.txt | parallel \"java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 Original_files/{}_1.fastq.gz Original_files/{}_2.fastq.gz trimming/{}_1.PEtrim.fastq.gz trimming/{}_1.UnPEtrim.fastq.gz  trimming/{}_2.PEtrim.fastq.gz trimming/{}_2.UnPEtrim.fastq.gz ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:13 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36\"")
# Note that the full path to the "adapters" folder within the Trimmomatic directory should be add in the option ILLUMINACLIP
# Alternatively, make a copy of this folder or the specific file in your wd

## QC of cleaned reads ####
dir.create("FastQC_clean")
system("fastqc -t 8 trimming/*.fastq.gz -o FastQC_clean")

## Mapping cleaned reads to the plantDB using Bowtie 2 #####
dir.create("bowtie2_plant")
# mapping to the plant database
setwd("plantDB") # run from the directory where the plant database is saved
system("cat ../sample.names.txt | parallel \"bowtie2 --threads 12 -x ITAG2.4_genomic -1 ../trimming/{}_1.PEtrim.fastq.gz -2 ../trimming/{}_2.PEtrim.fastq.gz -U ../trimming/{}_1.UnPEtrim.fastq.gz,../trimming/{}_2.UnPEtrim.fastq.gz --un-gz ../bowtie2_plant/{}_noPl.fastq.gz --un-conc-gz ../bowtie2_plant/{}_noPlconc.fastq.gz -S ../{}_plant.sam\"")

## De novo assembly of filtered reads using Trinity ####
dir.create("Trinity")
# Join PE and single reads in an unique file per each sample
system("cat sample.names.txt | parallel \"zcat bowtie2_plant/{}_noPlconc.1.fastq.gz bowtie2_plant/{}_noPlconc.2.fastq.gz bowtie2_plant/{}_noPl.fastq.gz >bowtie2_plant/{}_allNoPl.fastq\"")
# Run Trinity for all samples:
for (i in 1:length(mySamples)){
  system(sprintf("Trinity --seqType fq --single bowtie2_plant/%s_allNoPl.fastq --run_as_paired --CPU 12 --max_memory 50G --output Trinity/trinity_%s", mySamples[i],mySamples[i]))
}

## Mapping contigs against a nt virus database using BLASTn ####
dir.create("blastn")
# blast output in tab format (6):
system("cat sample.names.txt | parallel \"blastn -query Trinity/trinity_{}.Trinity.fasta -db virusDB/plantVDB -outfmt '6 qseqid qlen sseqid slen qstart qend length pident evalue' -num_threads 12 -max_target_seqs 1 -word_size 7 > blastn/{}_VirNt.txt\"")

## load BLASTn outputs in R:
myblastnList <- list.files(path="blastn", pattern="*.txt", full.names = TRUE)

for (i in 1:length(mySamples)) assign (mySamples[i], read.delim(myblastnList[i], header=FALSE, stringsAsFactors = FALSE,
                                                                col.names=c("qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "length", "pident", "evalue")))

myBlastList <- mget(mySamples) # Load all tables in a list
myBlastFil <- myBlastList #duplicate the variable

## Filter BLASTn results 
for (i in 1:length(mySamples)) {
  myBlastFil[[i]] <- subset(myBlastFil[[i]], myBlastFil[[i]]['qlen']>500 & 
                              myBlastFil[[i]]['qlen']<14000 & 
                              myBlastFil[[i]]['length']>300 & 
                              myBlastFil[[i]]['evalue']<0.00001)
}
# export results into individual csv files
lapply(names(myBlastFil), function(x) {
  f <- myBlastFil[[x]]
  write.csv(f, file=paste0("blastn/",x,"_blastnFilt.csv"),row.names=FALSE)
})

## Get the sequences of the reference viruses found by BLASTn ####
# Get accessions of the identified viruses
virus <- character()
for (i in 1:length(mySamples)){
  virus[i] <- rbind(myBlastFil[[i]][3])
}
myVirus <- unique(unlist(virus)) # transform list to vector
myVirus <- gsub("^\\w+\\||\\|","",myVirus) # remove additional characters in the accession
# write virus accessions to a file
write(myVirus, "myVirusIDs.txt")
# extract the sequences from the pvDB.fasta file using SeqKit:
system("seqkit grep -f myVirusIDs.txt virusDB/pvDB.fasta > virusDB/VirusFound.fasta")

## Get the sequences of the viral contigs ####
# First add a specific identifier to each contig according to the sample from which they come from
myBlastFil2 <- myBlastFil # duplicate
for (i in 1:length(mySamples)){
  myBlastFil2[[i]][1] <- sapply(myBlastFil2[[i]][1],function(x) paste0(x, "_",mySamples[i]))
}
# Get the id of the contigs from BLASTn tables after filtering results
contigs <- character()
for (i in 1:length(mySamples)){
  contigs[i] <- rbind(myBlastFil2[[i]][1])
}
myContigs <- unique(unlist(contigs)) # transform list to vector
# write contig ids into a file
write(myContigs, "myContigsNames.txt")
# Change also the name in a copy of the contig FASTA files
for (i in 1:length(mySamples)) {
  system(sprintf("sed 's/ .*//' Trinity/trinity_%s.Trinity.fasta > Trinity/%s.temp.fasta",mySamples[i],mySamples[i]))
  system(sprintf("sed 's/>.*/&_%s/' Trinity/%s.temp.fasta > Trinity/%s.temp2.fasta" ,mySamples[i],mySamples[i],mySamples[i]))
}
# Save all contigs in a common file for all samples
system("cat Trinity/*.temp2.fasta > Trinity/all_contigs_Trinity.fasta")
# extract the sequences from the virus DB file using SeqKit
system("seqkit grep -f myContigsNames.txt Trinity/all_contigs_Trinity.fasta > VirusContigs.fasta")

## Create index files to run BWA using the viruses found as reference #### 
system("bwa index virusDB/VirusFound.fasta")

## Map reads against the identified viruses using BWA MEM algorithm ####
dir.create("bwa")
system("cat sample.names.txt | parallel \"bwa mem virusDB/VirusFound.fasta bowtie2_plant/{}_allNoPl.fastq > bwa/{}_bwaVirusFound.sam\"")

## Transform sam to bam, sort bam files and create indexed bam files ####
system("cat sample.names.txt | parallel \"samtools view -b -S -h bwa/{}_bwaVirusFound.sam > bwa/{}_bwaVirusFound.bam\"")
system("cat sample.names.txt | parallel \"samtools sort bwa/{}_bwaVirusFound.bam -o bwa/{}_bwaVirusFound.sort.bam\"")
system("cat sample.names.txt | parallel \"samtools index bwa/{}_bwaVirusFound.sort.bam\"")

## Counting reads ####
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

## Count reads mapping to reference viruses using SeqKit and the BAM files generated before
system("cat sample.names.txt | parallel \"seqkit bam -c bwa/{}_virus_reads.txt bwa/{}_bwaVirusFound.sort.bam\"")
# Load files in R
myCountNames<-paste("COUNT", mySamples, sep="_")
myCountFiles <- list.files(path="bwa", pattern="*virus_reads.txt", full.names = TRUE)
for (i in 1:length(myCountNames)) assign(myCountNames[i], read.csv(myCountFiles[i], sep="\t", header=TRUE, stringsAsFactors = FALSE,
                                                                   col.names=c("Accession","Counts","Second_counts","Super_counts")))
# Load in a list
myCountList <-  mget(myCountNames)
myCountList <- lapply(myCountList, function(x) x <- x[rowSums(x[,c(2:4)])>0,]) # remove not mapping data
myCountList <- lapply(myCountList, function(x) x[,c(1:2)]) # select only primary counts in column 2
# Get total number of viral reads for each sample in a data frame
virusT <- data.frame("Viral_reads"=unlist(lapply(myCountList, function(x) x <- colSums(x["Counts"]))))
# add percentages
virusT$Percentage <- round((virusT$Viral_reads/filterT$Filtered_reads)*100,2)

## Join all in a table
myCountRes<- cbind(rawT,cleanT,filterT,virusT)
# in case no filtered reads are available use:
# myCountRes<- cbind(rawT,cleanT,virusT)
rownames(myCountRes) <- unique(sub('\\_1.fastq.gz|_2.fastq.gz',"",myfileRaw))
# Export the table to a file
write.csv(myCountRes,"Summary_of_mapping_reads.csv")

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

# merge mySTAT with myCountList (created before) (don't forget to load purrr package)
mySTATtmp <- map2(mySTAT , myCountList, left_join , by="Accession")
# merge with myVirusInfoList
mySTATfin <- map2(mySTATtmp , myVirusInfoList , left_join , by="Accession")
names(mySTATfin) <- mySamples
# save results into individual csv files
lapply(names(mySTATfin), function(x) {
  f <- mySTATfin[[x]]
  write.csv(f, file=paste0(x,"_final_results.csv"),row.names=FALSE)
}) 

## Create plots of read coverage ####
# create a function to make the same plot for different samples (don't forget to load dplyr and ggplot2 packages)
myPlotFunc <- function(y) {   # where y is myDEPTH_list
  y <- left_join(y,myVirusInfo,"Accession") %>% left_join(.,myVirusLengths, "Accession") # Join myVirusInfo with myVirusLengths (created previously)
  y$Virus <- paste(y$Species,y$Segment,y$Accession,paste0("(",y$Virus_length," bp)"))   # Join Species, segment and accession in a new column called Virus
  y$Virus <- as.factor(y$Virus)           # transform Virus column to factor
  y$Virus <- factor(y$Virus,levels = rev(levels(y$Virus))) # order the levels alphabetically
  Mycolor   <- "blue3"                # color for plotting reads
    ggplot(y, aes(x=Virus, y=Position))+    # plot Virus Vs Position
      geom_boxplot()+
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
