#Nov 16,2023
#To get 99% from RAndom sampling (A part of this code extracted from my "Random_Large_4521.R" code)
#In FishData/New_dataset/RANDOM folder

#Read the fish dataset in R
fish_4520 <- read.table("Fish_Realigned_multigene.phy")

# Install from CRAN
#install.packages("tidyverse")
library(tidyverse)
#Rename columns
colnames(fish_4520) <- c("species_name","seq")

#New dataframe by removing first row as it consists of number of species and the multigene seq length
fish_multigene <- fish_4520[-1,]

head(fish_multigene)

#Extract the range of COI gene to be used in EPA
fishsamplefull_co1 <- substr(fish_multigene$seq,2710,3391)
#Check the class
class(fishsamplefull_co1)
#Convert to dataframe
fishsamplefull_co1 <- data.frame(fishsamplefull_co1)
#Check the class
class(fishsamplefull_co1)
#Combine columns of COI dataframe and fishsample100 dataframe to get the sepecies name
fishsamplefull_co1_with_speciesname <- cbind(fish_multigene$species_name,fishsamplefull_co1)
#Assign column names
colnames(fishsamplefull_co1_with_speciesname) <- c("species_name","seq")
library(phylotools)
#Convert to phylip format
dat2phylip(fishsamplefull_co1_with_speciesname,outfile= "fishsample4520_co1_with_speciesname.phy")
4520*(99/100)#4475

4520*(80/100)#3616

4520*(60/100)#2712

4520*(40/100)#1808

4520*(20/100)#904


fishsamples <- function (fishsample_n,y){
  #fishsample_y <- fishsample100[sample(nrow(fishsample100), y),]
  #generate 10 replicatesout without replacement
  fishsample_y <- as.data.frame(replicate(10, fishsample_n[sample(nrow(fishsample_n),y, replace=F),])) 
  return(fishsample_y)
}
set.seed(12)
#Use the function "fishsamples" to generate 80% dataset
fishsample80F <- fishsamples(fish_multigene,3616)
dim(fishsample80F)
#Separate dataframes in columns of mass fishsample80 dataframe
replicate_dataframesF <- c(fishsample80F$V1,fishsample80F$V2,fishsample80F$V3,fishsample80F$V4,fishsample80F$V5,fishsample80F$V6,fishsample80F$V7,fishsample80F$V8,fishsample80F$V9,fishsample80F$V10)
#function to generate dataframes from replicates
df <- function(rep_df){
  t8 <- as.data.frame(rep_df)
  return(t8)
}
# Use the function to generate replicate dataframes 
t9 <- df(replicate_dataframesF)
#Separte mass dataframe into a list
lst1 <- lapply(seq(1, ncol(t9), by=2), function(i) 
  t9[i: pmin((i+1), ncol(t9))])
#to get individual dataframes for each replicate
list2env(setNames(lst1, paste0("newdfF", seq_along(lst1))),
         envir=.GlobalEnv)
#Function to generate different level samples without replicates and replaces
fishsamples_no_rep <- function (fishsample_m,z){
  fishsample_z <- as.data.frame(fishsample_m[sample(nrow(fishsample_m), z),])
  #fishsample_y <- as.data.frame(replicate(10, fishsample_n[sample(nrow(fishsample_n),y, replace=F),])) 
  return(fishsample_z)
}
#Use the function to genetare phylip files from the dataframe
dff <- function(xx){
  
  return(dat2phylip(xx,outfile = "newphy.phy"))
  
}

#Check the dataset
length(unique(newdfF2$species_name))

#80% 1st dataframe convert to phylip
#dff(newdfF1)
dat2phylip(newdfF1,outfile= "fishsample80_1.phy")
library(seqRFLP)
dataframe2fas(newdfF1, file = "fishsample80_1.fasta")
fishsample60_1 <- fishsamples_no_rep(newdfF1,2712)
dat2phylip(fishsample60_1,outfile= "fishsample60_1.phy")
dataframe2fas(fishsample60_1, file = "fishsample60_1.fasta")
#dff(fishsample60_1)
fishsample40_1 <- fishsamples_no_rep(fishsample60_1,1808)
dataframe2fas(fishsample40_1, file = "fishsample40_1.fasta")
#dff(fishsample40_1)
dat2phylip(fishsample40_1,outfile= "fishsample40_1.phy")
fishsample20_1 <- fishsamples_no_rep(fishsample40_1,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_1,outfile= "fishsample20_1.phy")
dataframe2fas(fishsample20_1, file = "fishsample20_1.fasta")

#80% 2nd df
dat2phylip(newdfF2,outfile= "fishsample80_2.phy")
dataframe2fas(newdfF2, file = "fishsample80_2.fasta")
fishsample60_2 <- fishsamples_no_rep(newdfF2,2712)
dat2phylip(fishsample60_2,outfile= "fishsample60_2.phy")
dataframe2fas(fishsample60_2, file = "fishsample60_2.fasta")
#dff(fishsample60_1)
fishsample40_2 <- fishsamples_no_rep(fishsample60_2,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_2,outfile= "fishsample40_2.phy")
dataframe2fas(fishsample40_2, file = "fishsample40_2.fasta")
fishsample20_2 <- fishsamples_no_rep(fishsample40_2,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_2,outfile= "fishsample20_2.phy")
dataframe2fas(fishsample20_2, file = "fishsample20_2.fasta")

dat2phylip(newdfF3,outfile= "fishsample80_3.phy")
dataframe2fas(newdfF3, file = "fishsample80_3.fasta")
fishsample60_3 <- fishsamples_no_rep(newdfF3,2712)
dat2phylip(fishsample60_3,outfile= "fishsample60_3.phy")
dataframe2fas(fishsample60_3, file = "fishsample60_3.fasta")
#dff(fishsample60_1)
fishsample40_3 <- fishsamples_no_rep(fishsample60_3,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_3,outfile= "fishsample40_3.phy")
dataframe2fas(fishsample40_3, file = "fishsample40_3.fasta")
fishsample20_3 <- fishsamples_no_rep(fishsample40_3,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_3,outfile= "fishsample20_3.phy")
dataframe2fas(fishsample20_3, file = "fishsample20_3.fasta")

dat2phylip(newdfF4,outfile= "fishsample80_4.phy")
dataframe2fas(newdfF4, file = "fishsample80_4.fasta")
fishsample60_4 <- fishsamples_no_rep(newdfF4,2712)
dat2phylip(fishsample60_4,outfile= "fishsample60_4.phy")
dataframe2fas(fishsample60_4, file = "fishsample60_4.fasta")
#dff(fishsample60_1)
fishsample40_4 <- fishsamples_no_rep(fishsample60_4,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_4,outfile= "fishsample40_4.phy")
dataframe2fas(fishsample40_4, file = "fishsample40_4.fasta")
fishsample20_4 <- fishsamples_no_rep(fishsample40_4,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_4,outfile= "fishsample20_4.phy")
dataframe2fas(fishsample20_4, file = "fishsample20_4.fasta")

dat2phylip(newdfF5,outfile= "fishsample80_5.phy")
dataframe2fas(newdfF5, file = "fishsample80_5.fasta")
fishsample60_5 <- fishsamples_no_rep(newdfF5,2712)
dat2phylip(fishsample60_5,outfile= "fishsample60_5.phy")
dataframe2fas(fishsample60_5, file = "fishsample60_5.fasta")
#dff(fishsample60_1)
fishsample40_5 <- fishsamples_no_rep(fishsample60_5,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_5,outfile= "fishsample40_5.phy")
dataframe2fas(fishsample40_5, file = "fishsample40_5.fasta")
fishsample20_5 <- fishsamples_no_rep(fishsample40_5,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_5,outfile= "fishsample20_5.phy")
dataframe2fas(fishsample20_5, file = "fishsample20_5.fasta")

dat2phylip(newdfF6,outfile= "fishsample80_6.phy")
dataframe2fas(newdfF6, file = "fishsample80_6.fasta")
fishsample60_6 <- fishsamples_no_rep(newdfF6,2712)
dat2phylip(fishsample60_6,outfile= "fishsample60_6.phy")
dataframe2fas(fishsample60_6, file = "fishsample60_6.fasta")
#dff(fishsample60_1)
fishsample40_6 <- fishsamples_no_rep(fishsample60_6,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_6,outfile= "fishsample40_6.phy")
dataframe2fas(fishsample40_6, file = "fishsample40_6.fasta")
fishsample20_6 <- fishsamples_no_rep(fishsample40_6,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_6,outfile= "fishsample20_6.phy")
dataframe2fas(fishsample20_6, file = "fishsample20_6.fasta")
dat2phylip(newdfF7,outfile= "fishsample80_7.phy")
dataframe2fas(newdfF7, file = "fishsample80_7.fasta")
fishsample60_7 <- fishsamples_no_rep(newdfF7,2712)
dataframe2fas(fishsample60_7, file = "fishsample60_7.fasta")
dat2phylip(fishsample60_7,outfile= "fishsample60_7.phy")
#dff(fishsample60_1)
fishsample40_7 <- fishsamples_no_rep(fishsample60_7,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_7,outfile= "fishsample40_7.phy")
dataframe2fas(fishsample40_7, file = "fishsample40_7.fasta")
fishsample20_7 <- fishsamples_no_rep(fishsample40_7,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_7,outfile= "fishsample20_7.phy")
dataframe2fas(fishsample20_7, file = "fishsample20_7.fasta")
dat2phylip(newdfF8,outfile= "fishsample80_8.phy")
dataframe2fas(newdfF8, file = "fishsample80_8.fasta")
fishsample60_8 <- fishsamples_no_rep(newdfF8,2712)
dat2phylip(fishsample60_8,outfile= "fishsample60_8.phy")
dataframe2fas(fishsample60_8, file = "fishsample60_8.fasta")
#dff(fishsample60_1)
fishsample40_8 <- fishsamples_no_rep(fishsample60_8,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_8,outfile= "fishsample40_8.phy")
dataframe2fas(fishsample40_8, file = "fishsample40_8.fasta")
fishsample20_8 <- fishsamples_no_rep(fishsample40_8,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_8,outfile= "fishsample20_8.phy")
dataframe2fas(fishsample20_8, file = "fishsample20_8.fasta")
dat2phylip(newdfF9,outfile= "fishsample80_9.phy")
dataframe2fas(newdfF9, file = "fishsample80_9.fasta")
fishsample60_9 <- fishsamples_no_rep(newdfF9,2712)
dat2phylip(fishsample60_9,outfile= "fishsample60_9.phy")
dataframe2fas(fishsample60_9, file = "fishsample60_9.fasta")
#dff(fishsample60_1)
fishsample40_9 <- fishsamples_no_rep(fishsample60_9,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_9,outfile= "fishsample40_9.phy")
dataframe2fas(fishsample40_9, file = "fishsample40_9.fasta")
fishsample20_9 <- fishsamples_no_rep(fishsample40_9,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_9,outfile= "fishsample20_9.phy")
dataframe2fas(fishsample20_9, file = "fishsample20_9.fasta")

dat2phylip(newdfF10,outfile= "fishsample80_10.phy")
dataframe2fas(newdfF10, file = "fishsample80_10.fasta")
fishsample60_10 <- fishsamples_no_rep(newdfF10,2712)
dat2phylip(fishsample60_10,outfile= "fishsample60_10.phy")
dataframe2fas(fishsample60_10, file = "fishsample60_10.fasta")
#dff(fishsample60_1)
fishsample40_10 <- fishsamples_no_rep(fishsample60_10,1808)
#dff(fishsample40_1)
dat2phylip(fishsample40_10,outfile= "fishsample40_10.phy")
dataframe2fas(fishsample40_10, file = "fishsample40_10.fasta")
fishsample20_10 <- fishsamples_no_rep(fishsample40_10,904)
#dff(fishsample20_1)
dat2phylip(fishsample20_10,outfile= "fishsample20_10.phy")
dataframe2fas(fishsample20_10, file = "fishsample20_10.fasta")

#newdfF1 (multigene 80%) and fish_multigene (Reference dataset) antijoin to get 20% of sequences that exist in Fish_multigene reference dataset not in multigene 80%
fish_multi20_1 <- anti_join(fish_multigene,newdfF1,by="species_name")
names(newdfF2)
names(fishsamplefull_co1_with_speciesname)
colnames(newdfF2) <- c("species_name","seq")
colnames(newdfF3) <- c("species_name","seq")
colnames(newdfF4) <- c("species_name","seq")
colnames(newdfF5) <- c("species_name","seq")
colnames(newdfF6) <- c("species_name","seq")
colnames(newdfF7) <- c("species_name","seq")
colnames(newdfF8) <- c("species_name","seq")
colnames(newdfF9) <- c("species_name","seq")
colnames(newdfF10) <- c("species_name","seq")
fish_multi20_2 <- anti_join(fish_multigene,newdfF2,by="species_name")
fish_multi20_3 <- anti_join(fish_multigene,newdfF3,by="species_name")
fish_multi20_4 <- anti_join(fish_multigene,newdfF4,by="species_name")
fish_multi20_5 <- anti_join(fish_multigene,newdfF5,by="species_name")
fish_multi20_6 <- anti_join(fish_multigene,newdfF6,by="species_name")
fish_multi20_7 <- anti_join(fish_multigene,newdfF7,by="species_name")
fish_multi20_8 <- anti_join(fish_multigene,newdfF8,by="species_name")
fish_multi20_9 <- anti_join(fish_multigene,newdfF9,by="species_name")
fish_multi20_10 <- anti_join(fish_multigene,newdfF10,by="species_name")

#now to get 99% do a rbind to merge two dataframes 80% and 19% (this 19% collected from the removed 20% in the previous step)
set.seed(100)
fishsample19_1<- fishsamples_no_rep(fish_multi20_1,859)
fishsample19_2<- fishsamples_no_rep(fish_multi20_2,859)
fishsample19_3<- fishsamples_no_rep(fish_multi20_3,859)
fishsample19_4<- fishsamples_no_rep(fish_multi20_4,859)
fishsample19_5<- fishsamples_no_rep(fish_multi20_5,859)
fishsample19_6<- fishsamples_no_rep(fish_multi20_6,859)
fishsample19_7<- fishsamples_no_rep(fish_multi20_7,859)
fishsample19_8<- fishsamples_no_rep(fish_multi20_8,859)
fishsample19_9<- fishsamples_no_rep(fish_multi20_9,859)
fishsample19_10<- fishsamples_no_rep(fish_multi20_10,859)
#rowbind of two dataframes

fishsample99_1 <- rbind(newdfF1,fishsample19_1)
fishsample99_2 <- rbind(newdfF2,fishsample19_2)
fishsample99_3 <- rbind(newdfF3,fishsample19_3)
fishsample99_4 <- rbind(newdfF4,fishsample19_4)
fishsample99_5 <- rbind(newdfF5,fishsample19_5)
fishsample99_6 <- rbind(newdfF6,fishsample19_6)
fishsample99_7 <- rbind(newdfF7,fishsample19_7)
fishsample99_8 <- rbind(newdfF8,fishsample19_8)
fishsample99_9 <- rbind(newdfF9,fishsample19_9)
fishsample99_10 <- rbind(newdfF10,fishsample19_10)

#fishsample99 (multigene 99%) and fishsamplefull_co1_with_speciesname (Co1 100% marker) antijoin to get 1% of co1 that exist only in co1 100% not in multigene 99%
CO1_1_1 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_1,by="species_name")
CO1_1_2 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_2,by="species_name")
CO1_1_3 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_3,by="species_name")
CO1_1_4 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_4,by="species_name")
CO1_1_5 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_5,by="species_name")
CO1_1_6 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_6,by="species_name")
CO1_1_7 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_7,by="species_name")
CO1_1_8 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_8,by="species_name")
CO1_1_9 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_9,by="species_name")
CO1_1_10 <- anti_join(fishsamplefull_co1_with_speciesname,fishsample99_10,by="species_name")


#Genearte phylip files for aligned co1 dataset 
#Just change the dataset name (CO1_20_1...CO1_80_10) for the function dff (function dff leads to datarame to phylip format)
dff(CO1_1_10)

dff(fishsample99_10)

#Generate fasta files for aligned co1 dataset
dataframe2fas(CO1_1_1, file = "CO1_1_1.fasta")
dataframe2fas(CO1_1_2, file = "CO1_1_2.fasta")
dataframe2fas(CO1_1_3, file = "CO1_1_3.fasta")
dataframe2fas(CO1_1_4, file = "CO1_1_4.fasta")
dataframe2fas(CO1_1_5, file = "CO1_1_5.fasta")
dataframe2fas(CO1_1_6, file = "CO1_1_6.fasta")
dataframe2fas(CO1_1_7, file = "CO1_1_7.fasta")
dataframe2fas(CO1_1_8, file = "CO1_1_8.fasta")
dataframe2fas(CO1_1_9, file = "CO1_1_9.fasta")
dataframe2fas(CO1_1_10, file = "CO1_1_10.fasta")

dataframe2fas(fishsample99_1, file = "fishsample99_1.fasta")
dataframe2fas(fishsample99_2, file = "fishsample99_2.fasta")
dataframe2fas(fishsample99_3, file = "fishsample99_3.fasta")
dataframe2fas(fishsample99_4, file = "fishsample99_4.fasta")
dataframe2fas(fishsample99_5, file = "fishsample99_5.fasta")
dataframe2fas(fishsample99_6, file = "fishsample99_6.fasta")
dataframe2fas(fishsample99_7, file = "fishsample99_7.fasta")
dataframe2fas(fishsample99_8, file = "fishsample99_8.fasta")
dataframe2fas(fishsample99_9, file = "fishsample99_9.fasta")
dataframe2fas(fishsample99_10, file = "fishsample99_10.fasta")

#Convert fas to phylip format
#install.packages("seqmagick")
library(seqmagick)
library(seqinr)
library(phylotools)
#Here fish20_co1_80_1.fasta is the reference aligned sequence set which is an output of multigene 20% and co1 80% aligning by using muscle in Linux terminal
getwd()
