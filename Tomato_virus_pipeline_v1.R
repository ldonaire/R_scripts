# Commands in R for Tomato virus pipeline (reference publication Maachi et al. 2021)
# wrote by Livia Donaire and Ayoub Maachi 
# version1 March 2021

# Example: raw reads of one sample called sampleA (PE, 150 bp)

########## Requiered programs ######

# Java OpenJDK Runtime Environment (build 1.8.0_152-release-1056-b12)
# fastqc version 0.11.9
# Trimmomatic 0.39
# SeqTK Version 1.2-r94 https://github.com/lh3/seqtk 
# bbmap, 38.82 https://sourceforge.net/projects/bbmap/
# bowtie2 version 2.4.1 
# Trinity version 2.8.4 https://github.com/trinityrnaseq/trinityrnaseq/releases 
# quast-5.0.0 https://github.com/ablab/quast 
# blast 2.10.0 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
# Tablet 1.19.09.03 https://ics.hutton.ac.uk/tablet/download-tablet/
# bwa Version: 0.7.17-r1188 
# samtools Version: 1.7 
# bcftools Version 1.6
# pyfasta 0.5.2 https://pypi.org/project/pyfasta/#files 
# cd-hit

########## Requiered R packages #####

# install.packages("dplyr")
# install.packages("data.table")
# install.packages("tidyr")
# install.packages("tidyverse")

########## Creating subsets of raw reads #########

# To create a random subset of your reads
system("seqtk sample -s 100 path_to_original_dir/raw_reads_sampleA_1.fq n > path_to_output_dir/sampleA_1_subset.fq")
system("seqtk sample -s 100 path_to_original_dir/raw_reads_sampleA_2.fq n > path_to_output_dir/sampleA_2_subset.fq")

# n = number of substracted reads
# NOTE: for PE reads use the same ramdom seed (-s option) to keep pairing

########## Quality control of raw reads #####

# using FastQC v0.11.9
system("fastqc -t n path_to_original_directory/sampleA*.fq")

# option -t (number of threads used to run the analysis)
# visualize the results by opening the html file in a web browser and decide the next step according to the quality of the files

########## Trimming raw reads ######

# using Trimmomatic 0.39
system("java -jar path_to_trimmomatics/trimmomatic-0.39.jar SE -phred33 path_to_original_directory/sampleA_1_subset.fq path_to_trimming_dir/sampleA_1.trim.fq ILLUMINACLIP:path_to_adapter_file/fastqc_adapters.fa:2:30:10  LEADING:n TRAILING:x HEADCROP:y MINLEN:70 AVGQUAL:30")
system("java -jar path_to_trimmomatics/trimmomatic-0.39.jar SE -phred33 path_to_original_directory/sampleA_2_subset.fq path_to_trimming_dir/sampleA_2.trim.fq ILLUMINACLIP:path_to_adapter_file/fastqc_adapters.fa:2:30:10  LEADING:n TRAILING:x HEADCROP:y MINLEN:70 AVGQUAL:30")

# n, x, y = integers indicating the number of nucleotides to crop according to the FastQC report.

# Trimmomatic options:
# SE: treat the samples as if they where SE and the pairs will be repaired above
# -phred33 - quality score of fastq files (Illumina since 2015) https://drive5.com/usearch/manual/quality_score.html
# ILLUMINACLIP: cut adapter and other illumina-specific sequences from the read.
# LEADING: cut bases off the start of a read, if below a threshold quality (Q30)
# TRAILING: Cut bases off the end of a read, if below a threshold quality (Q30)
# HEADCROP: Cut the specified number of bases from the start of the read (decide according to per base sequence content 
# analysis of FastQC)
# MINLEN: Dropthe read if it is below a specified length (depending of your read lenght)
# AVGQUAL: Drop the read if the average quality is below the specified level (Q30)


########## Repair paired end data (only for PE reads) ######

# using bbmap v38.82:
system (" /path_to_bbmap_dir/repair.sh in1=path_to_trimming_dir/sampleA_1.trim.fq in2=path_to_trimming_dir/sampleA_1.trim.fq out1=path_to_trimming_dir/sampleA_1.repaired.fq out2=path_to_trimming_dir/sampleA_2.repaired.fq outsingle=path_to_trimming_dir/sampleA.single.fq")


########## Filtering out reads mapping with the plant genome #######

###### Build bowtie2 database 
system("bowtie2-build path_to_plantDB_dir/plant_genome.fasta path_to_plantDB_dir/plant_genome")

# this command will generate 6 files with .bt2 extension in plantDB directory

###### Mapping reads to the plantDB using bowtie2 
system("bowtie2 -x path_to_plantDB_dir/plant_genome -1 path_to_trimming_dir/sampleA_1.repaired.fq -2 path_to_trimming_dir/sampleA_2.repaired.fq -U path_to_trimming_dir/sampleA.single.fq  --un-conc path_to_bowtie2_dir/sampleA_noPlant.fq -S sampleA_delete.sam")

# Note: the output file sampleA_delete.sam contains reads mapping to the plant genome and can be removed

########## De novo assembly of reads ########################################

# de novo assembly using Trinity v2.8.4
system("path_to_Trinity_dir/Trinity --seqType fq --single path_to_bowtie2_dir/sampleA_*.fq  --no_normalize_reads --run_as_paired --CPU n --max_memory x")

# the output is a folder called "trinity_out_dir"

########## Quality of assembly using QUAST #####
### using quast v5.0.0:
system("python path_to_quast/quast.py -o ./QUAST_output_dir/Trinity_assembly --rna-finding --eukaryote --plots-format pdf path_to_trinity_out_dir")

# see report in html or pdf

########## Mapping contigs to the virus sequences using blastn ###########################
###### Build Blast nt database
## Removing low complex sequences
system("dustmasker -in path_to_virus_seq/plant_virus_seq.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out path_to_virus_seq/plant_virus_seq.asnb")
## Creating the database
system("makeblastdb -in path_to_virus_seq/plant_virus_seq.fa -input_type fasta -dbtype nt -parse_seqids -mask_data path_to_virus_seq/plant_virus_seq.asnb -out path_to_virus_seq/plantVirusDB -title 'Plant viruses nucleotide database, v1'")

# fasta headers should contain information about "Species", "Genus", "Family" and "Segment"

##### Get DB info:
system("blastdbcmd -db path_to_virus_seq/plantVirusDB -info")

##### mapping using Blastn (output in tab format)
system("blastn -query path_to_trinity_out_dir/contigs.fa -db path_to_virus_seq/plantVirusDB -outfmt '6 qseqid qlen sseqid stitle slen qstart qend length pident evalue' -num_threads n -max_target_seqs 1 -max_hsps 1 > path_to_blastn_out_dir/sampleA_virus.txt")

##### Filtering blastn results 
# load blastn output into R:
sampleA_blast<- read.delim(sampleA_virus.txt, header=FALSE, stringsAsFactors = FALSE,
                                                            col.names=c("qseqid", "qlen", "sseqid", "stitle", "slen", "qstart", "qend", "length", "pident", "evalue"))
# separate column stitle by comma:
library(tidyr)
sampleA <- separate(SampleA_blast,"stitle",into=c("Species", "Genus", "Family","Segment"),sep=",")
# Filtering results:
myFilter <- function(x){
  x <- x[x$evalue<0.00001,] # evalue < 0.00001
  x <- x[x$qlen>1000 & x$qlen<14000,] # contig length bewteen 1000 and 14000
  x <- x[x$length>500,]
}
sampleAfil <- myFilter(sampleA)
write.csv(sampleAfil, "path_to_blastn_out_dir/sampleA_blastn_filtered.csv", row.names=FALSE)

########## Realingning reads against the reference viruses found #################

##### Get the sequences of the viruses found by blastn
# Get id of viruses with mapping reads from blastn tables loaded in step 4.2 and remove special character
myVirus <- unique(sampleAfil[3])
myVirus$sseqid <- gsub("ref","",myVirus$sseqid) 
myVirus$sseqid <- gsub("[/|]","",myVirus$sseqid) 

# write virus ids into a file
library("data.table")
fwrite(myVirus, "path_to_virus_seq/myVirusNames.txt",col.names=FALSE)

# in order to substract fasta sequences it is necessary to remove all the information in the fasta header after the first whitespace:
system("cut -d ' ' -f1 path_to_virus_seq/plant_virus_seq.fa  >path_to_virus_seq/plant_virus_seq.tmp.fa")

# extract the sequences from the virus DB in fasta format file using pyfasta v0.5.2:
system("pyfasta extract --header --fasta path_to_virus_seq/plant_virus_seq.tmp.fa --file path_to_virus_seq/myVirusNames.txt > path_to_virus_seq/Viruses_found.fa")

###### Create index files to run BWA 
system("bwa index path_to_virus_seq/Viruses_found.fa")

###### Map reads against new Virus found Database using bwa mem algorithm 
system ("bwa mem -t 12  path_to_virus_seq/Viruses_found path_to_bowtie2_dir/sampleA_noPlant.1.fq path_to_bowtie2_dir/sampleA_noPlant.1.fq >path_to_bwa_dir/sampleA_bwa_map.sam")

# Optional: visualize the results using Tablet and the sam file

########## Calculation of genome coverage and average read depth ##########

# Transform to bam, sort the file and create index using samtools v1.7
system("samtools view -b -S -h path_to_bwa_dir/sampleA_bwa_map.sam > path_to_bwa_dir/sampleA_bwa_map.bam")
system("samtools sort path_to_bwa_dir/sampleA_bwa_map.bam -o path_to_bwa_dir/sampleA_bwa_map.sort.bam")
system("samtools index path_to_bwa_dir/sampleA_bwa_map.sort.bam")

# Get sequencing depth using samtools depth 
system("samtools depth path_to_bwa_dir/sampleA_bwa_map.sort.bam > path_to_bwa_dir/sampleA_bwa_depth.txt")

# Get the lenght of the viruses found using samtools faidx
system("samtools faidx path_to_virus_seq/Viruses_found.fa")
myVirusLengths <- read.csv("path_to_virus_seq/Viruses_found.fa.fai", sep="\t", header=FALSE, stringsAsFactors = FALSE)
# Lenghts are column 2
myVirusLengths <- myVirusLengths[,1:2]
colnames(myVirusLengths) <- c("ID","Virus_length")

# load files in R to make the calculations
library(dplyr)
myDepthSampleA<-read.csv("path_to_bwa_dir/sampleA_bwa_depth.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("ID","Position","Reads"))

myCoverage<-function(x,myVirusLenghts){
  myMeandDepth <- aggregate(cbind(x$V3)~x$V1, x, mean)
  colnames(myMeandDepth) <- c("ID","Average_depth")
  myPosCov <- aggregate(cbind(x$V3)~x$V1, x, length)
  colnames(myPosCov) <- c("ID","Length_covered")
  myStats <- full_join(myMeandDepth,myPosCov,"ID") %>% full_join(.,myVirusLengths,"ID")
  myStats <- transform(myStats,Perc_covered=(myStats$Length_covered/myStats$Virus_length)*100)
}

# Apply the function myCoverage 
mySTAT_sampleA <- myCoverage(myDepthSampleA)

# Add virus information:
# requiered a csv file with virus information (ID, Species, Genus, Family, Segment)
VirusData <- read.csv("path_to_virus_seq/plant_virus_seq_info.csv", col.names=c("ID","Species","Genus","Family","Segment"))

# join with mySTAT:
mySTAT <- left_join(mySTAT_sampleA,VirusData,"ID")
write.csv(mySTAT, file="path_to_bwa_dir/sampleA_calculations.txt", row.names=FALSE)

########## Getting  consensus sequences of viruses ##########
# Using samtools and bcftools
system("samtools mpileup -uf path_to_virus_seq/Viruses_found.fa path_to_bwa_dir/sampleA_bwa_map.sort.bam | bcftools call -c | vcfutils.pl vcf2fq > path_to_cns_out_folder/virus_sampleA_cns.fastq")
system("seqtk seq -aQ64 -q20 -n N path_to_cns_out_folder/virus_sampleA_cns.fastq > path_to_cns_out_folder/virus_sampleA_cns.fasta")


########## Discovery of new viruses by BLASTX ################################
##### Build blast protein database 
## Collapsing protein sequences with aa identity higher than 98% 
system("cd-hit -i path_to_virus_seq/plant_virus_prot.fa -o path_to_virus_seq/plant_virus_prot_collapsed.fa -c 0.98 -n 5 -M 16000 -d 0 -T 8")
## Removing low complex sequences
system("dustmasker -in path_to_virus_seq/plant_virus_prot_collapsed.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out path_to_virus_seq/plant_virus_prot_collapsed.asnb")
## Removing low complex sequences
system("makeblastdb -in path_to_virus_seq/plant_virus_prot_collapsed.fa -input_type fasta -dbtype prot -parse_seqids -mask_data path_to_virus_seq/plant_virus_prot_collapsed.asnb -out path_to_virus_seq/plant_virus_prot_collapsed_DB -title 'Plant virus protein database, collapsed, version 1'")

##### Mapping contigs using blastx
system("blastx -query path_to_trinity_out_dir/contigs.fa -db path_to_virus_seq/plant_virus_prot_collapsed_DB -outfmt '6 qseqid qlen sseqid stitle slen qstart qend length pident evalue' -num_threads 12 -max_target_seqs 1 -max_hsps 1 > path_to_blastx_out_dir/sampleA_NewViruses.txt")
