# ab <- read.csv(file.choose())
# ab <- read.csv(file.choose(),header=T)

#The path of FullStats.tsv dataset is Desktop/FishData/New_dataset
setwd("~/Desktop/FishData/New_dataset")
NewFull_dataset<-readr::read_tsv("FullStats.tsv")
# library(readr)
# bb <- read_tsv("FullStats.tsv")
head(NewFull_dataset)
dim(NewFull_dataset)
names(NewFull_dataset)
NewFull_ungapped_length <- NewFull_dataset[c("species","zic_ungapped_length","vcpip_ungapped_length","twelveS_ungapped_length","tbr_ungapped_length","svep_ungapped_length","sreb_ungapped_length","sixteenS_ungapped_length","sidkey_ungapped_length","shepex_ungapped_length","ripk_ungapped_length","rhod_ungapped_length","ragtwo_ungapped_length","ragone_ungapped_length","ptr_ungapped_length","plagl_ungapped_length","panx_ungapped_length","ndtwo_ungapped_length","ndfour_ungapped_length","myh_ungapped_length","kiaa_ungapped_length","hox_ungapped_length","glyt_ungapped_length","fourC_ungapped_length","ficd_ungapped_length","encone_ungapped_length","cytb_ungapped_length","coi_ungapped_length")]
dim(NewFull_ungapped_length)
head(NewFull_ungapped_length)
names(NewFull_ungapped_length)[1] <- "species_name"
names(NewFull_ungapped_length)
#Read the fish dataset in R
fish <- read.table("final_alignment.phylip")
head(fish,1)
#remove the 1st row with numbers
fish <- fish[-1,]
head(fish,1)
class(fish)

#Filter out sequences without CO1 gene as we need only the multi-gene sequences including CO1 for downstream analysis 
#Rename columns
colnames(fish) <- c("species_name","seq")

#check the dataset
dim(fish)
length(unique(fish$species_name))

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
#Convert to phylip format to be used in EPA
dat2phylip(fishCO1_with_speciesname,outfile= "fishsamplefull_co1_with_speciesname.phy")
#Rename column names
colnames(df_2) <- c("species_name","seq")
#Keep only rows with sequences
fish_multigene <- merge(fish,df_2, by= "species_name")
#New dataframe with first two columns
fish_multigene <- fish_multigene[, 1:2]
#Rename columns
colnames(fish_multigene) <- c("species_name","seq")


#-------
#This part is to check whether there are species without molecular data
#Didn't find any species without molecular data
fish_new <- fish
#Replace "-" with nothing 
df_0 <- gsub('-','',fish_new$seq)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
dim(df_0)
# #Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(fish_new$species_name,df_0)
dim(df_1)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
dim(df_2)

#-----
#combine fish dataset with NewFull_ungapped_length dataset to have species names and sequences for the downstream analysis
NewFull_ungapped_length <- cbind(fish,NewFull_ungapped_length)

#Remove species with less than 3 genes
NewFull_zeros_less_three <- NewFull_ungapped_length[rowSums(NewFull_ungapped_length == 0) <= 24, ]


library(readr)
#write the dataframe to a csv file
write_csv(NewFull_zeros_less_three,path = "~/Desktop/FishData/New_dataset/dataset_more_three_genes.csv")

#extract the species with co1 gene as I need co1 genes to place on the bb trees
# multigene_include_co1 <- NewFull_zeros_less_three[!(NewFull_zeros_less_three$coi_ungapped_length==0),]
#multigene_include_co1 <- merge(fish_multigene,NewFull_zeros_less_three, by= "species_name")
#Extract CO1 nucleotides
fishCO1 <- substr(NewFull_zeros_less_three$seq,2292,2973)
#Check the class of fishCO1 dataset
class(fishCO1)
#Convert to dataframe
fishCO1 <- data.frame(fishCO1)
#Check the converted dataset
class(fishCO1)
#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(NewFull_zeros_less_three$species_name,fishCO1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(NewFull_zeros_less_three$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
fishsamplefull_co1_with_speciesname <- df_2
#Convert to phylip format to be used in EPA
#dat2phylip(fishCO1_with_speciesname,outfile= "fishsamplefull_co1_with_speciesname.phy")
#Rename column names
colnames(df_2) <- c("species_name","seq")
#Keep only rows with sequences
#fish_multigene <- merge(fish,df_2, by= "species_name")
fish_multigenenew <- merge(NewFull_zeros_less_three,df_2, by= "species_name")
#New dataframe with first two columns
fish_multigenenew_sp_seq <- fish_multigenenew[, 1:2]
#Rename columns
colnames(fish_multigenenew_sp_seq) <- c("species_name","seq")

library(devtools)
#install_github("helixcn/seqRFLP")
library(seqRFLP)
dataframe2fas(fish_multigenenew_sp_seq, file = "multigene_with_co1.fasta")

#write the data frame to a csv file
write_csv(fish_multigenenew,path = "~/Desktop/FishData/New_dataset/multigene_with_co1_new.csv")

#----From what Tyler sent after checking my multigene_with_co1.csv file
FullStats_three<-readr::read_tsv("FullStats_three.tsv")
dim(FullStats_three)
names(FullStats_three)

#The data set I have to use for the realignment purpose is "multigene_include_co1"
names(multigene_include_co1)
Fish_4856 <- subset(multigene_include_co1[,1:2])
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
aligned12S_False <- AlignSeqs(seqs_12S,useStructures=FALSE)
# view the alignment in a browser 
#BrowseSeqs(aligned12S, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned12S,
                file="Aligned_12S")
writeXStringSet(aligned12S_False,
                file="Aligned_12S_False")
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
#-----
#After the alignment Tyler check the stats for 12S and 16S
#read the 12S and 16S stats tables in to R to compare Cr score with the Cr score of previous alignment done by FishTree of Life authors

Stats_12S<-readr::read_csv("12S.Table_1.csv")
dim(Stats_12S)
names(Stats_12S)
Stats_16S<-readr::read_csv("16S.Table_1.csv")
dim(Stats_16S)
names(Stats_16S)
#By checking the new alignment of 12S I found some species that have extra bases which makes all the other sequences gappy
#So I compare the Cr scores with the old alignment
#Stats_12S is the stats dataframe of new alignment
#FullStats_three is the stats dataframe of old alignment
Stats_12S[Stats_12S$Sequence=="Dicologlossa_hexophthalma"|Stats_12S$Sequence=="Zebrasoma_scopas"|Stats_12S$Sequence=="Microchirus_azevia"|Stats_12S$Sequence=="Triportheus_angulatus"|Stats_12S$Sequence=="Pagellus_bellottii"|Stats_12S$Sequence=="Nemipterus_japonicus",c(2,4)]
FullStats_three[FullStats_three$species=="Dicologlossa_hexophthalma"|FullStats_three$species=="Zebrasoma_scopas"|FullStats_three$species=="Microchirus_azevia"|FullStats_three$species=="Triportheus_angulatus"|FullStats_three$species=="Pagellus_bellottii"|FullStats_three$species=="Nemipterus_japonicus",c(111,138)]
names(FullStats_three)
#FullStats_three[FullStats_three$twelve_Cr<0.5 & FullStats_three$twelve_Cr>0,c(111,138)]

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
names(FullStats_three)
FullStats_three$ndfour_length_per
FullStats_three[FullStats_three$ndfour_length_per<0.3 & FullStats_three$ndfour_length_per>0,c(111,90)]
cytb_len_per <- FullStats_three[FullStats_three$cytb_length_per<0.3 & FullStats_three$cytb_length_per>0,c(111,131)]
write.csv(cytb_len_per,file="cytb_len_per.csv")
sixteenS_len_per <- FullStats_three[FullStats_three$sixteenS_length_per<0.4 & FullStats_three$sixteenS_length_per>0,c(111,35)]
write.csv(cytb_len_per,file="sixteenS_len_per.csv")


#extra----
Co1fasta <- read_fasta("COI_realign_with_replacements.fasta")


#--------
gg<-read.csv("multigene_with_co1.csv")#incorrect file
cc<-read.csv("multigene_with_co1_new.csv")#incorrect file
