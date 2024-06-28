#!/usr/bin/env Rscript
library(dplyr)
args <- commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please, add a file with one specie name / line and the path /path/Project/", call.=FALSE)
}
##### Libraries #############################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("BiocInstaller")
#BiocManager::install("S4Vectors")
#BiocManager::install("IRanges", version = "3.8")
#biocLite("DESeq2")

#BiocInstaller::biocLite("GenomicRanges")
#BiocManager::install("S4Vectors")
library(S4Vectors)
library(GenomicRanges)
#library(dplyr)

#############################################
#"/beegfs/home/bguinet/M2_script/short_file_species_name.txt"
path <- args[1] #Fichier avec tous les noms d'espèces pour boucler
#path<-"/beegfs/home/bguinet/M2_script/short_file_species_name.txt"
path2 <- args[2] #Fichier avec tous les noms d'espèces pour boucler
#Permet de récupérer chacune des lignes dans un vecteur names
names=c()
conn <- file(path,open="r")
lines <- readLines(conn)
cat ("\n")
cat("Files ready to by analysed: ")
cat ("\n")
for (i in 1:length(lines)){
  names[i] <- lines[i]
  print(lines[i])
}
cat ("\n")
cat(" Overlapping processing... ")
cat ("\n")


for (i in names){

  ##############################
  ###. SUMMARY BUSCO. ##########
  ##############################
  
  
#  tab_busco=read.table(paste0(path2,i,"/run_busco/run_BUSCO_v3/full_table_",i,"_BUSCO_v3_strand.tsv"), sep="\t",header=T,fill = TRUE,na.strings=c("", "NA"))
##  tab_busco$X<-NULL
 # colnames(tab_busco) <- c("query", "Status", "target","sstart",  "ssend", "Score", "Length","Strand")
  #reduce is to use target to split the ranges into groups, then calculate the range of each group, and simplify that as much as possible
  #So I created a GRanges with target as an additonal column
 # tab_busco<-tab_busco[!(tab_busco$Status=="Missing"),]
 # tab_busco<-tab_busco[!(tab_busco$Strand=="NA"),]
 # tab_busco[rowSums(is.na(tab_busco)) != ncol(tab_busco), ]
 # tab_busco<-tab_busco[complete.cases(tab_busco), ]
 # genes_busco <- with(tab_busco, GRanges(seqnames=target,ranges=IRanges(start=sstart,end=ssend), target=query,strand=Strand))
  # Allows to creat a Granges object and add the sseqif informations in order to sort the tab with it.\n)
 # tab_busco <- reduce(genes_busco) #Merge all overlapping sequences
#  tab_busco <- data.frame(tab_busco)
  #2Resume_scaffold_busco=count(tab_busco,tab_busco$seqnames)
  
 # write.table(tab_busco, file=paste0(path2,i,"/run_busco/result_blast_summary_H.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote=F)
 #2 write.table(Resume_scaffold_busco, file=paste0(path2,i,"/run_busco/resume_scaffold_busco_H.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote=F)


  ##############################
  ###. SUMMARY VIRAL. ##########
  ##############################
  
  mmseq_viral_tab <- read.table(paste0(path2,i,"/run_mmseqs2_V/result_mmseqs2_strand_V.m8"),header=T,sep='')
  
  #reduce is to use target to split the ranges into groups, then calculate the range of each group, and simplify that as much as possible
  #So I created a GRanges with target as an additonal column 
  genes_viral1 <- with(mmseq_viral_tab, GRanges(seqnames=query,ranges=IRanges(start=qstart,end=qend), target=target,strand=strand,len=qlen))
  # Allows to creat a Granges object and add the sseqif informations in order to sort the tab with it.\n")
  mmseq_viral_tab1 <- reduce(genes_viral1) #Merge all overlapping sequences
  mmseq_viral_tab1 <- data.frame(mmseq_viral_tab1)
  
  genes_viral <- with(mmseq_viral_tab1, GRanges(seqnames=seqnames,ranges=IRanges(start=start,end = end)))
  genes_viral$strand2 <- mmseq_viral_tab1$strand
  genes_viraldf <- data.frame(genes_viral)
  genes_viraldf$strand<-NULL
  
  #mmseq_viral_tab<-read.table(paste0(path2,i,"/run_mmseq_V/Matches_",i,"_summary_V.txt"),header = T)
  write.table(genes_viral, file=paste0(path2,i,"/run_mmseqs2_V/result_mmseqs2_summary_V.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote=F)
  
  #################################################################
  #### REMOVING THE MERGING PARTS BETWEEN THE TWO DataFrames ######
  #################################################################
  
  #genes_busco <- with(tab_busco, GRanges(seqnames=seqnames,ranges=IRanges(start=start,end = end)))
  #genes_viral <- with(mmseq_viral_tab1, GRanges(seqnames=seqnames,ranges=IRanges(start=start,end = end)))
  #genes_viral$numbers <- row(data.frame(genes_viral))[,1] #Allows to add a metadata columns with number from 1 to length (genes_viral)
  #genes_viral$type <- "V" #Allows to add a metadata columns with data type
  #genes_busco$type <- "H" #Allows to add a metadata columns with data type
  #genes_viral$strand2 <- mmseq_viral_tab1$strand
  #genes_viral$len <- mmseq_viral_tab1$len
  #
  #hits <- findOverlaps(genes_viral, genes_busco)
  #overlaps <- pintersect(genes_viral[queryHits(hits)], genes_busco[subjectHits(hits)]) #Allow to found overlapps between busco and viral locus 
  #overlaps$with <- width(overlaps)#Add the width of the overlapp between viral and busco locus
  #overlaps$with_viral <- width(genes_viral[queryHits(hits)]) #Add the width of the viral locus
  #overlaps$with_busco <- width(genes_busco[subjectHits(hits)]) #Add the width of the busco locus 
  #overlaps$percentOverlap <- width(overlaps) / width(genes_viral[queryHits(hits)]) #Add the percentage of the overlapp that correspond to the viral locus length
  
  #overlapsdf <- data.frame(overlaps)#Convert GR object to a df
 #overlapsdf <- overlapsdf[overlapsdf$percentOverlap %in% overlapsdf$percentOverlap[overlapsdf$percentOverlap>0.1],] #Allows to keep overlapping seq below 10%
  
  #genes_viraldf <- data.frame(genes_viral)#Convert GR object to a df
  #Without_overlapps <- genes_viraldf[!genes_viraldf$numbers %in% overlapsdf$numbers, ] #Allows to only keep sequence that do not overlapp
  
  #2Without_overlapps$len <- mmseq_viral_tab$qlen[match(Without_overlapps$seqnames, mmseq_viral_tab$query)]
  #2Without_overlapps2<-Without_overlapps %>%
  #2 mutate(Newqstart = case_when(strand2 == '+' ~ start, TRUE ~ len - end),
  #2       Newqend= case_when(strand2 == '-' ~ len - start, TRUE ~ end))
  
  #2Resume_scaffold_viral=count(Without_overlapps,Without_overlapps$seqnames)
  #2write.table(Resume_scaffold_viral, file=paste0(path2,i,"/run_mmseqs2_V/resume_scaffold_viral_V.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote=F)
  
 
  #write.table(overlapsdf , file=paste0(path2,i,"/Matches_",i,"_with_overlapping_sequences.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE,quote=F)
  #write.table(Without_overlapps, file=paste0(path2,i,"/Matches_",i,"_without_overlapping_sequences.m8"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE,quote=F)
  print(i)
}






