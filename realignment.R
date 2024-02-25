library(phylotools)
fasta.12S = read.fasta("part1.12S_realigned.fst")
fasta.16S = read.fasta("part2.16S_realigned.fst")
fasta.4c4 = read.fasta("part3.4c4_realigned_manual_edits.fst")
fasta.COI = read.fasta("part4.COI_manual_edits.fst")
fasta.CytB = read.fasta("part5.CytB_realigned.fst")
fasta.Enc1 = read.fasta("part6.Enc1_realigned.fst")
fasta.Ficd = read.fasta("part7.Ficd_manual_edits.fst")
fasta.Glyt = read.fasta("part8.Glyt_manual_edits.fst")
fasta.Hoxc6a = read.fasta("part9.Hoxc6a.fst")
fasta.Kiaa1239 = read.fasta("part10.Kiaa1239_manual_edits.fst")
fasta.Myh6 = read.fasta("part11.Myh6_manual_edits.fst")
fasta.Nd2 = read.fasta("part12.Nd2_manual_edits.fst")
fasta.ND4 = read.fasta("part13.ND4_realigned.fst")
fasta.Panx2 = read.fasta("part14.Panx2.fst")
fasta.Plagl2 = read.fasta("part15.Plagl2_manual_edits.fst")
fasta.Ptr = read.fasta("part16.Ptr.fst")
fasta.Rag1 = read.fasta("part17.Rag1_realigned_manual_edits.fst")
fasta.Rag2 = read.fasta("part18.Rag2_manual_edits.fst")
fasta.Rhodopsin = read.fasta("part19.Rhodopsin.fst")
fasta.Ripk4 = read.fasta("part20.Ripk4.fst")
fasta.Sh3px3 = read.fasta("part21.Sh3px3.fst")
fasta.Sidkey = read.fasta("part22.Sidkey_manual_edits.fst")
fasta.Sreb2 = read.fasta("part23.Sreb2.fst")
fasta.Svep1 = read.fasta("part24.Svep1.fst")
fasta.Tbr = read.fasta("part25.Tbr_manual_edits.fst")
fasta.Vcpip = read.fasta("part26.Vcpip.fst")
fasta.Zic1 = read.fasta("part27.Zic1_manual_edits.fst")

#get the length of the sequences with gaps 
#Important to generate the alignment partition file. 
#Now Fish Tol alignment partition file is useless as we have realigned all the genes.
sapply(fasta.12S,nchar)
sapply(fasta.COI,nchar)
sapply(fasta.Sidkey,nchar)
sapply(fasta.Zic1,nchar)

#Check whether the species in the data frame are unique
length(unique(fasta.12S$species_name))
setdiff(fasta.12S$species_name,fasta.COI$seq.name)
# dim(fish_multigenenew_sp_seq_new)
# names(fish_multigenenew_sp_seq_new)
# names(fasta.12S) <- c("species_name","seq" )
# equ <- dplyr::all.equal(fasta.12S$species_name,fish_multigenenew_sp_seq_new$species_name)
# setdiff(fasta.12S$species_name,fish_multigenenew_sp_seq_new$species_name)
# anti_join(fasta.12S$species_name,fish_multigenenew_sp_seq_new$species_name)
# extra <- merge(fasta.12S,fish_multigenenew_sp_seq_new,by="species_name")
my_string <- "-"  
species <- strrep(my_string, 5) 

#----
#First merge all 27 dataframes into a one datframe using seq.name
#multigene_4856 <- merge(fasta.12S,fasta.16S,fasta.4c4,fasta.COI,fasta.CytB,fasta.Enc1,fasta.Ficd,fasta.Glyt,fasta.Hoxc6a,fasta.Kiaa1239,fasta.Myh6,fasta.Nd2,fasta.ND4,fasta.Panx2,fasta.Plagl2,fasta.Ptr,fasta.Rag1,fasta.Rag2,fasta.Rhodopsin,fasta.Ripk4,fasta.Sh3px3,fasta.Sidkey,fasta.Sreb2,fasta.Svep1,fasta.Tbr,fasta.Vcpip,fasta.Zic1,by="seq.name")
# Join multiple data.frames
library(tidyverse)
list_df = list(fasta.12S,fasta.16S,fasta.4c4,fasta.COI,fasta.CytB,fasta.Enc1,fasta.Ficd,fasta.Glyt,fasta.Hoxc6a,fasta.Kiaa1239,fasta.Myh6,fasta.Nd2,fasta.ND4,fasta.Panx2,fasta.Plagl2,fasta.Ptr,fasta.Rag1,fasta.Rag2,fasta.Rhodopsin,fasta.Ripk4,fasta.Sh3px3,fasta.Sidkey,fasta.Sreb2,fasta.Svep1,fasta.Tbr,fasta.Vcpip,fasta.Zic1)
#multigene_4856 <- list_df %>% reduce(merge, by="seq.name")
multigene_4856 <- list_df %>% reduce(inner_join, by="seq.name")
#assign column names 
colnames(multigene_4856) <- c("species_name","fasta.12S","fasta.16S","fasta.4c4","fasta.COI","fasta.CytB","fasta.Enc1","fasta.Ficd","fasta.Glyt","fasta.Hoxc6a","fasta.Kiaa1239","fasta.Myh6","fasta.Nd2","fasta.ND4","fasta.Panx2","fasta.Plagl2","fasta.Ptr","fasta.Rag1","fasta.Rag2","fasta.Rhodopsin","fasta.Ripk4","fasta.Sh3px3","fasta.Sidkey","fasta.Sreb2","fasta.Svep1","fasta.Tbr","fasta.Vcpip","fasta.Zic1")
# Concatanate the 27 genes in toa one column to make the sequence a multigene one

#multigene_4856$Seq <- multigene_4856 %>% paste(fasta.12S,fasta.16S,fasta.4c4,fasta.COI,fasta.CytB,fasta.Enc1,fasta.Ficd,fasta.Glyt,fasta.Hoxc6a,fasta.Kiaa1239,fasta.Myh6,fasta.Nd2,fasta.ND4,fasta.Panx2,fasta.Plagl2,fasta.Ptr,fasta.Rag1,fasta.Rag2,fasta.Rhodopsin,fasta.Ripk4,fasta.Sh3px3,fasta.Sidkey,fasta.Sreb2,fasta.Svep1,fasta.Tbr,fasta.Vcpip,fasta.Zic1)

library(tidyr)
multigene_4856 <- multigene_4856 %>% 
  unite("Merged",fasta.12S:fasta.Zic1)
multigene_4856$Merged

#remove _ in sequences of the multigene as they appear when concatanating 27 genes
# sapply(multigene_4856$Merged,nchar)
multigene_4856_new <- multigene_4856

multigene_4856_new$Merged <- gsub('_','',multigene_4856_new$Merged)
multigene_4856_new$Merged
sapply(multigene_4856_new$Merged,nchar)

library(readr)
#write the dataframe to a csv file
write_csv(multigene_4856_new,path = "~/Desktop/FishData/New_dataset/Realignments/multigene_4856_new.csv")

#Now remove the species with issues
#Read the csv file with species to be removed
Species_wt_issues <- read_csv("Species_wt_issues.csv")

class(Species_wt_issues)
Species_wt_issues <- as.data.frame(Species_wt_issues)
length(unique(Species_wt_issues$Species_name))#334
colnames(multigene_4856_new) <- c("Species_name","Seq")
Realigned_multigene <- anti_join(multigene_4856_new,Species_wt_issues,by="Species_name")#4521
names(multigene_4856_new)
class(Realigned_multigene)

library(stringr)
library(dplyr)

data1 <- Realigned_multigene[str_detect(Realigned_multigene$Seq, "z"), ] 
#As I found a bad base (z) for Dascyllus_aruanus, I remove that species from the dataset
Realigned_multigene <- Realigned_multigene[-3213,]

#Keep a csv file of Realigned_multigene
#write_csv(Realigned_multigene,path = "~/Desktop/FishData/New_dataset/Realignments/Fish_Realigned_multigene.csv")
#The `path` argument of `write_csv()` is deprecated as of readr 1.4.0.
#So, I used the `file` argument instead of path. Got the same results
write_csv(Realigned_multigene,file = "~/Desktop/FishData/New_dataset/Realignments/Fish_Realigned_multigene.csv")
#Keep a fasta file of Realigned_multigene
library(seqRFLP)
dataframe2fas(Realigned_multigene, file = "Fish_Realigned_multigene.fasta")
#Keep a phylip file of Realigned_multigene
dat2phylip(Realigned_multigene,outfile= "Fish_Realigned_multigene.phy")
getwd()
head(Realigned_multigene)
