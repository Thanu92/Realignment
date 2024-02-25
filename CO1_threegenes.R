#The path of FullStats.tsv dataset is Desktop/FishData/New_dataset
setwd("~/Desktop/FishData/New_dataset")
NewFull_dataset<-readr::read_tsv("FullStats.tsv")
# library(readr)
# bb <- read_tsv("FullStats.tsv")
head(NewFull_dataset)
dim(NewFull_dataset)
names(NewFull_dataset)[111] <- "species_name"
names(NewFull_dataset)
# NewFull_ungapped_length <- NewFull_dataset[c("species","zic_ungapped_length","vcpip_ungapped_length","twelveS_ungapped_length","tbr_ungapped_length","svep_ungapped_length","sreb_ungapped_length","sixteenS_ungapped_length","sidkey_ungapped_length","shepex_ungapped_length","ripk_ungapped_length","rhod_ungapped_length","ragtwo_ungapped_length","ragone_ungapped_length","ptr_ungapped_length","plagl_ungapped_length","panx_ungapped_length","ndtwo_ungapped_length","ndfour_ungapped_length","myh_ungapped_length","kiaa_ungapped_length","hox_ungapped_length","glyt_ungapped_length","fourC_ungapped_length","ficd_ungapped_length","encone_ungapped_length","cytb_ungapped_length","coi_ungapped_length")]
# dim(NewFull_ungapped_length)
# head(NewFull_ungapped_length)
# names(NewFull_ungapped_length)[1] <- "species_name"
# names(NewFull_ungapped_length)
# #Read the fish dataset in R
fish <- read.table("final_alignment.phylip")
head(fish,1)
#remove the 1st row with numbers
fish <- fish[-1,]
head(fish,1)
class(fish)

#Filter out sequences without CO1 gene as we need only the multi-gene sequences including CO1 for downstream analysis 
#Rename columns
colnames(fish) <- c("species_name","seq")
#merge fish dataframe and NewFull_ungapped_length
NewFull_dataset_new <- merge(fish,NewFull_dataset,by="species_name")
dim(NewFull_dataset_new)

NewFull_ungapped_length <- NewFull_dataset_new[c("species_name","seq","zic_ungapped_length","vcpip_ungapped_length","twelveS_ungapped_length","tbr_ungapped_length","svep_ungapped_length","sreb_ungapped_length","sixteenS_ungapped_length","sidkey_ungapped_length","shepex_ungapped_length","ripk_ungapped_length","rhod_ungapped_length","ragtwo_ungapped_length","ragone_ungapped_length","ptr_ungapped_length","plagl_ungapped_length","panx_ungapped_length","ndtwo_ungapped_length","ndfour_ungapped_length","myh_ungapped_length","kiaa_ungapped_length","hox_ungapped_length","glyt_ungapped_length","fourC_ungapped_length","ficd_ungapped_length","encone_ungapped_length","cytb_ungapped_length","coi_ungapped_length")]
dim(NewFull_ungapped_length)
head(NewFull_ungapped_length)

#Remove species with less than 3 genes
NewFull_zeros_less_three <- NewFull_ungapped_length[rowSums(NewFull_ungapped_length == 0) <= 24, ]
library(readr)
#write the dataframe to a csv file
write_csv(NewFull_zeros_less_three,path = "~/Desktop/FishData/New_dataset/dataset_more_three_genes.csv")


#Extract CO1 nucleotides
fishCO1 <- substr(fish$seq,2292,2973)
#Check the class of fishCO1 dataset
class(fishCO1)
#Convert to dataframe
fishCO1 <- data.frame(fishCO1)
#Check the converted dataset
class(fishCO1)
#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(fish$species_name,fishCO1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(fish$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_co1_with_speciesname <- df_2
names(fishsamplefull_co1_with_speciesname) <- c("species_name","co1Seq")
#get the species with more than three genes including co1 gene 
Fish_CO1_threegenes <- merge(NewFull_zeros_less_three,fishsamplefull_co1_with_speciesname,by="species_name")
dim(Fish_CO1_threegenes)
names(Fish_CO1_threegenes)

#New dataframe with first two columns
fish_multigenenew_sp_seq_new <- Fish_CO1_threegenes[, 1:2]
names(fish_multigenenew_sp_seq_new)

library(devtools)
#install_github("helixcn/seqRFLP")
library(seqRFLP)
dataframe2fas(fish_multigenenew_sp_seq_new, file = "multigeneCO1_threeGene.fasta")
#write the data frame to a csv file
write_csv(Fish_CO1_threegenes,path = "~/Desktop/FishData/New_dataset/multigene_CO1_threegenes.csv")

#check whether co1 is present in the sequences

#Extract CO1 nucleotides
fishCO1 <- substr(Fish_CO1_threegenes$seq,2292,2973)
#Check the class of fishCO1 dataset
class(fishCO1)
#Convert to dataframe
fishCO1 <- data.frame(fishCO1)
#Check the converted dataset
class(fishCO1)
#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(Fish_CO1_threegenes$species_name,fishCO1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_CO1_threegenes$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_co1_with_speciesname <- df_2

#-----
# abc <- merge(Fish_CO1_threegenes,cc,by="species_name")
# names(abc)
# df1 <- subset(abc[,c("zic_ungapped_length.x" , "zic_ungapped_length.y")])
# df2 <- subset(abc[,c("seq","seq.x")])
# data_common1 <- generics::intersect(Fish_CO1_threegenes, cc1)
# names(cc)
# names(cc)[2] <- "seq"
# names(cc)[30] <- "co1Seq"
# cc1 <- cc
# cc1[cc1 == "Abalistes_stellaris"] <- "Abalistes_stellaris_new"
# cc1[cc1 == "Abalistes_stellaris"] <- "Abalistes_stellaris_new"
# 
# library(phylotools)
# fasta.df = read.fasta("multigene_with_co1.fasta")
# names(fasta.df) <- c("species_name","seq")
# identical(fish_multigenenew_sp_seq_new,fasta.df)

#----
#The data set I have to use for the realignment purpose is "fish_multigenenew_sp_seq_new"
names(fish_multigenenew_sp_seq_new)
Fish_4856 <- fish_multigenenew_sp_seq_new
dim(Fish_4856)
names(Fish_4856)

#----12S
library(tidyverse)

#Extract 12S nucleotides
fish12S <- substr(Fish_4856$seq,1,979)
#Check the class of fishCO1 dataset
class(fish12S)
#Convert to dataframe
fish12S <- data.frame(fish12S)
#Check the converted dataset
class(fish12S)
#Bind the column "species name" to the fishCO1 dataframe
fish12S_with_speciesname <- cbind(Fish_4856$species_name,fish12S)
#Replace "-" with nothing 
df_0 <- gsub('-','',fish12S_with_speciesname$fish12S)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_12S_with_speciesname <- df_2

#library(phylotools)
#install.packages("devtools")
library(devtools)
#install_github("helixcn/seqRFLP")
library(seqRFLP)
dataframe2fas(fishsamplefull_12S_with_speciesname, file = "Fish12S.fasta")

#Now we have to align the ungapped 12S sequences 
#install.packages("DECIPHER")

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   + install.packages("BiocManager")
# BiocManager::install("DECIPHER")

# load the DECIPHER library in R
library(DECIPHER)

# load the sequences from the file
# changed "DNA" to "RNA" as I am using 12S??
#seqs_RNA <- readRNAStringSet("Fish12S.fasta")
seqs_12S <- readDNAStringSet("Fish12S.fasta")
# look at some of the sequences (optional)
seqs_12S

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned12S <- AlignSeqs(seqs_12S,useStructures=TRUE)
#aligned12S_False <- AlignSeqs(seqs_12S,useStructures=FALSE)
# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned12S,
                file="Aligned_12S")
#writeXStringSet(aligned12S_False,
#                file="Aligned_12S_False")
#----16S
#Extract 16S nucleotides
fish16S <- substr(Fish_4856$seq,980,1756)
#Check the class of fishCO1 dataset
class(fish16S)
#Convert to dataframe
fish16S <- data.frame(fish16S)
#Check the converted dataset
class(fish16S)
#Bind the column "species name" to the fish_4856 dataframe
fish16S_with_speciesname <- cbind(Fish_4856$species_name,fish16S)
#Replace "-" with nothing 
df_0 <- gsub('-','',fish16S_with_speciesname$fish16S)
names(fish16S_with_speciesname)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_16S_with_speciesname <- df_2

#library(phylotools)
#install.packages("devtools")
library(devtools)
#install_github("helixcn/seqRFLP")
library(seqRFLP)
dataframe2fas(fishsamplefull_16S_with_speciesname, file = "Fish16S.fasta")
#align 16S
# load the sequences from the file
# changed "DNA" to "RNA" as I am using 12S??
#seqs_RNA <- readRNAStringSet("Fish12S.fasta")
seqs_16S <- readDNAStringSet("Fish16S.fasta")
# look at some of the sequences (optional)
seqs_16S

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned16S <- AlignSeqs(seqs_16S,useStructures=TRUE)

# view the alignment in a browser 
#BrowseSeqs(aligned16S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned16S,
                file="Aligned_16S")

#Now align Nd4 and CytB
#----ND4-----
#Extract ND4 nucleotides
fishND4 <- substr(Fish_4856$seq,9869,11401)
#Check the class of fishCO1 dataset
class(fishND4)
#Convert to dataframe
fishND4 <- data.frame(fishND4)
#Check the converted dataset
class(fishND4)
#Bind the column "species name" to the fishCO1 dataframe
fishND4_with_speciesname <- cbind(Fish_4856$species_name,fishND4)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishND4_with_speciesname$fishND4)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_ND4_with_speciesname <- df_2

#library(phylotools)
#install.packages("devtools")
library(devtools)
#install_github("helixcn/seqRFLP")
library(seqRFLP)
dataframe2fas(fishsamplefull_ND4_with_speciesname, file = "FishND4.fasta")

#Now we have to align the ungapped ND4 sequences 

# load the sequences from the file
seqs_ND4 <- readDNAStringSet("FishND4.fasta")
# look at some of the sequences (optional)
seqs_ND4

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
#seqs <- OrientNucleotides(seqs)

# perform the alignment
alignedND4 <- AlignSeqs(seqs_ND4)

# write the alignment to a new FASTA file
writeXStringSet(alignedND4,
                file="Aligned_ND4")

#----CytB----
#Extract CytB nucleotides
fishCytB <- substr(Fish_4856$seq,2974,4114)
#Check the class of fishCO1 dataset
class(fishCytB)
#Convert to dataframe
fishCytB <- data.frame(fishCytB)
#Check the converted dataset
class(fishCytB)
#Bind the column "species name" to the fishCO1 dataframe
fishCytB_with_speciesname <- cbind(Fish_4856$species_name,fishCytB)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishCytB_with_speciesname$fishCytB)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_CytB_with_speciesname <- df_2

dataframe2fas(fishsamplefull_CytB_with_speciesname, file = "FishCytB.fasta")

# load the sequences from the file
seqs_CytB <- readDNAStringSet("FishCytB.fasta")
# look at some of the sequences (optional)
seqs_CytB 

# perform the alignment
alignedCytB <- AlignSeqs(seqs_CytB)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedCytB,
                file="Aligned_CytB")

#Check the stats for ND4 And CytB
names(NewFull_dataset)
NewFull_dataset$twelveS_length_per
TwelveS_len_per <- NewFull_dataset[NewFull_dataset$twelveS_length_per<0.2 & NewFull_dataset$twelveS_length_per>0,c(111,15)]
write.csv(TwelveS_len_per,file="TwelveS_len_per.csv")
SixteenS_len_per <- NewFull_dataset[NewFull_dataset$sixteenS_length_per<0.3 & NewFull_dataset$sixteenS_length_per>0,c(111,35)]
write.csv(SixteenS_len_per,file="SixteenS_len_per.csv")

cytb_len_per <- NewFull_dataset[NewFull_dataset$cytb_length_per<0.2 & NewFull_dataset$cytb_length_per>0,c(111,131)]
write.csv(cytb_len_per,file="cytb_len_per.csv")
ndfour_len_per <- NewFull_dataset[NewFull_dataset$ndfour_length_per<0.2 & NewFull_dataset$ndfour_length_per>0,c(111,90)]
write.csv(ndfour_len_per,file="ndfour_len_per.csv")

#co1
#Extract co1 nucleotides
fishco1 <- substr(Fish_4856$seq,2292,2973)
#Check the class of fishCO1 dataset
class(fishco1)
#Convert to dataframe
fishco1  <- data.frame(fishco1)
#Check the converted dataset
class(fishco1)
#Bind the column "species name" to the fishCO1 dataframe
fishco1_with_speciesname <- cbind(Fish_4856$species_name,fishco1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishco1_with_speciesname$fishco1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_co1_with_speciesname <- df_2

#4c4 gene
#----4c4----
#Extract 4c4 nucleotides
fish4c4 <- substr(Fish_4856$seq,1757,2291)
#Check the class of fishCO1 dataset
class(fish4c4)
#Convert to dataframe
fish4c4 <- data.frame(fish4c4)
#Check the converted dataset
class(fish4c4)
#Bind the column "species name" to the fishCO1 dataframe
fish4c4_with_speciesname <- cbind(Fish_4856$species_name,fish4c4)
#Replace "-" with nothing 
df_0 <- gsub('-','',fish4c4_with_speciesname$fish4c4)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_4c4_with_speciesname <- df_2

dataframe2fas(fishsamplefull_4c4_with_speciesname, file = "Fish4c4.fasta")

# load the sequences from the file
seqs_4c4 <- readDNAStringSet("Fish4c4.fasta")
# look at some of the sequences (optional)
seqs_4c4 

# perform the alignment
aligned4c4 <- AlignSeqs(seqs_4c4)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned4c4,
                file="Aligned_4c4")

#enc1 gene
#----enc1----
#Extract enc1 nucleotides
fishenc1 <- substr(Fish_4856$seq,4115,4954)
#Check the class of fishCO1 dataset
class(fishenc1)
#Convert to dataframe
fishenc1 <- data.frame(fishenc1)
#Check the converted dataset
class(fishenc1)
#Bind the column "species name" to the fishCO1 dataframe
fishenc1_with_speciesname <- cbind(Fish_4856$species_name,fishenc1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishenc1_with_speciesname$fishenc1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_enc1_with_speciesname <- df_2

dataframe2fas(fishsamplefull_enc1_with_speciesname, file = "Fishenc1.fasta")

# load the sequences from the file
seqs_enc1 <- readDNAStringSet("Fishenc1.fasta")
# look at some of the sequences (optional)
seqs_enc1 

# perform the alignment
alignedenc1 <- AlignSeqs(seqs_enc1)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedenc1,
                file="Aligned_enc1")

#ficd gene
#----ficd----
#Extract ficd nucleotides
fishficd <- substr(Fish_4856$seq,4955,5680)
#Check the class of fishCO1 dataset
class(fishficd)
#Convert to dataframe
fishficd <- data.frame(fishficd)
#Check the converted dataset
class(fishficd)
#Bind the column "species name" to the fishCO1 dataframe
fishficd_with_speciesname <- cbind(Fish_4856$species_name,fishficd)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishficd_with_speciesname$fishficd)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_fishficd_with_speciesname <- df_2

dataframe2fas(fishsamplefull_fishficd_with_speciesname, file = "Fishficd.fasta")

# load the sequences from the file
seqs_ficd <- readDNAStringSet("Fishficd.fasta")
# look at some of the sequences (optional)
seqs_ficd 

# perform the alignment
alignedficd <- AlignSeqs(seqs_ficd)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedficd,
                file="Aligned_ficd")


#glyt gene
#----glyt----
#Extract glyt nucleotides
fishglyt <- substr(Fish_4856$seq,5681,6613)
#Check the class of fishglyt dataset
class(fishglyt)
#Convert to dataframe
fishglyt <- data.frame(fishglyt)
#Check the converted dataset
class(fishglyt)
#Bind the column "species name" to the fishCO1 dataframe
fishglyt_with_speciesname <- cbind(Fish_4856$species_name,fishglyt)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishglyt_with_speciesname$fishglyt)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_fishglyt_with_speciesname <- df_2

dataframe2fas(fishsamplefull_fishglyt_with_speciesname, file = "Fishglyt.fasta")

# load the sequences from the file
seqs_glyt <- readDNAStringSet("Fishglyt.fasta")
# look at some of the sequences (optional)
seqs_glyt 

# perform the alignment
alignedglyt <- AlignSeqs(seqs_glyt)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedglyt,
                file="Aligned_glyt")

#hoxc6a gene
#----hoxc6a----
#Extract hoxc6a nucleotides
fishhoxc6a <- substr(Fish_4856$seq,6614,7243)
#Check the class of fishglyt dataset
class(fishhoxc6a)
#Convert to dataframe
fishhoxc6a <- data.frame(fishhoxc6a)
#Check the converted dataset
class(fishhoxc6a)
#Bind the column "species name" to the fishCO1 dataframe
fishhoxc6a_with_speciesname <- cbind(Fish_4856$species_name,fishhoxc6a)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishhoxc6a_with_speciesname$fishhoxc6a)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_fishhoxc6a_with_speciesname <- df_2

dataframe2fas(fishsamplefull_fishhoxc6a_with_speciesname, file = "Fishhoxc6a.fasta")

# load the sequences from the file
seqs_hoxc6a <- readDNAStringSet("Fishhoxc6a.fasta")
# look at some of the sequences (optional)
seqs_hoxc6a

# perform the alignment
alignedhoxc6a <- AlignSeqs(seqs_hoxc6a)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedhoxc6a,
                file="Aligned_hoxc6a")

#----kiaa1239----
#Extract kiaa1239 nucleotides
fishkiaa1239 <- substr(Fish_4856$seq,7244,7987)
#Check the class of fishkiaa1239 dataset
class(fishkiaa1239)
#Convert to dataframe
fishkiaa1239 <- data.frame(fishkiaa1239)
#Check the converted dataset
class(fishkiaa1239)
#Bind the column "species name" to the fishCO1 dataframe
fishkiaa1239_with_speciesname <- cbind(Fish_4856$species_name,fishkiaa1239)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishkiaa1239_with_speciesname$fishkiaa1239)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(Fish_4856$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_fishkiaa1239_with_speciesname <- df_2

dataframe2fas(fishsamplefull_fishkiaa1239_with_speciesname, file = "Fishkiaa1239.fasta")

# load the sequences from the file
seqs_kiaa1239 <- readDNAStringSet("Fishkiaa1239.fasta")
# look at some of the sequences (optional)
seqs_kiaa1239

# perform the alignment
alignedkiaa1239 <- AlignSeqs(seqs_kiaa1239)

# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(alignedkiaa1239,
                file="Aligned_kiaa1239")
kiaa_length_per100
fourC_length_per 116
hox_length_per 105
glyt_length_per 110
ficd_length_per 121
encone_length_per 126
names(NewFull_dataset)
kiaa_length_per<- NewFull_dataset[NewFull_dataset$kiaa_length_per<0.5 & NewFull_dataset$kiaa_length_per>0,c(111,100)] #0 objects
fourC_length_per <- NewFull_dataset[NewFull_dataset$fourC_length_per<0.3 & NewFull_dataset$fourC_length_per>0,c(111,116)]
write.csv(fourC_length_per,file="fourC_length_per.csv")
glyt_length_per <- NewFull_dataset[NewFull_dataset$glyt_length_per<0.4 & NewFull_dataset$glyt_length_per>0,c(111,110)]#0 objects


hox_length_per <- NewFull_dataset[NewFull_dataset$hox_length_per<0.3 & NewFull_dataset$hox_length_per>0,c(111,105)]
write.csv(hox_length_per,file="hox_length_per.csv")

ficd_length_per <- NewFull_dataset[NewFull_dataset$ficd_length_per<0.4 & NewFull_dataset$ficd_length_per>0,c(111,121)]
write.csv(ficd_length_per,file="ficd_length_per.csv")

encone_length_per<- NewFull_dataset[NewFull_dataset$encone_length_per<0.4 & NewFull_dataset$encone_length_per>0,c(111,126)]
write.csv(encone_length_per,file="encone_length_per.csv")

