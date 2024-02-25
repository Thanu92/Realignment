library(phylotools)
fasta.12S = read.fasta(file = "part1.12S_realigned -WITH_ROOT_SEQUENCES.fst")
sapply(fasta.12S,nchar)
# Using subset
sub12S <- subset(fasta.12S, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))
class(fasta.12S)
fasta.16S = read.fasta(file = "part2.16S_realigned_WITH_ROOT_SEQUENCES.fst")
sub16S <- subset(fasta.16S, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))
sapply(fasta.16S,nchar)
fasta.4c4 = read.fasta("part3.4c4_realigned_manual_edits.fst")
sapply(fasta.4c4,nchar)
fasta.COI = read.fasta(file = "part4.COI_manual_edits_WITH_ROOT_SEQUENCES.fst")
subCOI <- subset(fasta.COI, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))

fasta.CytB = read.fasta(file = "part5.CytB_realigned_WITH_ROOT_SEQUENCES.fst")
subCytB <- subset(fasta.CytB, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))
fasta.Enc1 = read.fasta("part6.Enc1_realigned.fst")
sapply(fasta.Enc1,nchar)
fasta.Ficd = read.fasta("part7.Ficd_manual_edits.fst")
sapply(fasta.Ficd,nchar)
fasta.Glyt = read.fasta("part8.Glyt_manual_edits.fst")
sapply(fasta.Glyt,nchar)
fasta.Hoxc6a = read.fasta("part9.Hoxc6a.fst")
sapply(fasta.Hoxc6a,nchar)
fasta.Kiaa1239 = read.fasta("part10.Kiaa1239_manual_edits.fst")
sapply(fasta.Kiaa1239,nchar)
fasta.Myh6 = read.fasta("part11.Myh6_manual_edits.fst")
sapply(fasta.Myh6,nchar)
fasta.Nd2 = read.fasta(file = "part12.Nd2_manual_edits_WITH_ROOT_SEQUENCES.fst")
sapply(fasta.Nd2,nchar)
subNd2 <- subset(fasta.Nd2, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))

fasta.ND4 = read.fasta("part13.ND4_realigned.fst")
sapply(fasta.ND4,nchar)
fasta.Panx2 = read.fasta("part14.Panx2.fst")
sapply(fasta.Panx2,nchar)
fasta.Plagl2 = read.fasta("part15.Plagl2_manual_edits.fst")
sapply(fasta.Plagl2,nchar)
fasta.Ptr = read.fasta("part16.Ptr.fst")
sapply(fasta.Ptr,nchar)
fasta.Rag1 = read.fasta(file = "part17.Rag1_realigned_manual_edits_WITH_ROOT_SEQUENCES.fst")
sapply(fasta.Rag1,nchar)
subRag1 <- subset(fasta.Rag1, seq.name %in% c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus"))

fasta.Rag2 = read.fasta("part18.Rag2_manual_edits.fst")
sapply(fasta.Rag2,nchar)
fasta.Rhodopsin = read.fasta("part19.Rhodopsin.fst")
sapply(fasta.Rhodopsin,nchar)
fasta.Ripk4 = read.fasta("part20.Ripk4.fst")
sapply(fasta.Ripk4,nchar)
fasta.Sh3px3 = read.fasta("part21.Sh3px3.fst")
sapply(fasta.Sh3px3,nchar)
fasta.Sidkey = read.fasta("part22.Sidkey_manual_edits.fst")
sapply(fasta.Sidkey,nchar)
fasta.Sreb2 = read.fasta("part23.Sreb2.fst")
sapply(fasta.Sreb2,nchar)
fasta.Svep1 = read.fasta("part24.Svep1.fst")
sapply(fasta.Svep1,nchar)
fasta.Tbr = read.fasta("part25.Tbr_manual_edits.fst")
sapply(fasta.Tbr,nchar)
fasta.Vcpip = read.fasta("part26.Vcpip.fst")
sapply(fasta.Vcpip,nchar)
fasta.Zic1 = read.fasta("part27.Zic1_manual_edits.fst")


#get the length of the sequences with gaps 
#Important to generate the alignment partition file. 
#Now Fish Tol alignment partition file is useless as we have realigned all the genes.
sapply(fasta.CytB,nchar)
sapply(fasta.COI,nchar)
sapply(fasta.Sidkey,nchar)
sapply(fasta.Zic1,nchar)

#create dataframes with dashes for outgroup sequences to show genes without bases to have the 27 aligned genes when combine 27 dataframes
my_string <- "-"  

txt_4c4 <- strrep(my_string, 535) 
txt_Enc1 <- strrep(my_string, 840) 
txt_Sidkey <- strrep(my_string, 1216) 
txt_Zic1 <- strrep(my_string, 1026) 
txt_vcpip <- strrep(my_string, 765) 
txt_Tbr <- strrep(my_string, 798) 
txt_Svep1 <- strrep(my_string, 807) 
txt_Sreb2 <- strrep(my_string, 987) 
txt_Sh3px3 <- strrep(my_string, 723) 
txt_Ripk4 <- strrep(my_string, 645) 
txt_Phodopsin <- strrep(my_string, 924) 
txt_Rag2 <- strrep(my_string, 1224) 
txt_Ptr <- strrep(my_string, 711) 
txt_Plagl2 <- strrep(my_string, 804) 
txt_Panx2 <- strrep(my_string, 740)
txt_ND4 <- strrep(my_string, 1667) 
txt_Myh6 <- strrep(my_string, 834) 
txt_Kiaa1239 <- strrep(my_string, 744) 
txt_Hoxc6a <- strrep(my_string, 630) 
txt_Glyt <- strrep(my_string, 933) 
txt_Ficd <- strrep(my_string, 726) 
 



seq.name <- c("Xenopus_laevis","Gallus_gallus","Homo_sapiens","Isurus_paucus","Pseudotriakis_microdon","Alopias_pelagicus","Triakis_semifasciata","Carcharhinus_plumbeus")
seq.txt_12S <- c(replicate(8,txt_4c4))
seq.txt_4c4 <- c(replicate(8,txt_4c4)) 
seq.txt_Enc1 <- c(replicate(8,txt_Enc1)) 
seq.txt_Sidkey <- c(replicate(8,txt_Sidkey)) 
seq.txt_Zic1 <- c(replicate(8,txt_Zic1))
seq.txt_vcpip <- c(replicate(8,txt_vcpip)) 
seq.txt_Tbr <- c(replicate(8,txt_Tbr))
seq.txt_Svep1 <- c(replicate(8,txt_Svep1)) 
seq.txt_Sreb2 <- c(replicate(8,txt_Sreb2))
seq.txt_Sh3px3 <- c(replicate(8,txt_Sh3px3)) 
seq.txt_Ripk4 <- c(replicate(8,txt_Ripk4))
seq.txt_Phodopsin <- c(replicate(8,txt_Phodopsin)) 
seq.txt_Rag2 <- c(replicate(8,txt_Rag2))
seq.txt_Ptr <- c(replicate(8,txt_Ptr))
seq.txt_Plagl2 <- c(replicate(8,txt_Plagl2)) 
seq.txt_Panx2 <- c(replicate(8,txt_Panx2))
seq.txt_ND4 <- c(replicate(8,txt_ND4))
seq.txt_Myh6 <- c(replicate(8,txt_Myh6)) 
seq.txt_Kiaa1239 <- c(replicate(8,txt_Kiaa1239)) 
seq.txt_Hoxc6a <- c(replicate(8,txt_Hoxc6a))
seq.txt_Glyt <- c(replicate(8,txt_Glyt))
seq.txt_Ficd <- c(replicate(8,txt_Ficd))



sub4c4 <- data.frame(seq.name,seq.txt_4c4)
subEnc1 <- data.frame(seq.name,seq.txt_Enc1)
subSidkey <- data.frame(seq.name,seq.txt_Sidkey)
subZic1 <- data.frame(seq.name,seq.txt_Zic1)
subvcpip <- data.frame(seq.name,seq.txt_vcpip)
subTbr <- data.frame(seq.name,seq.txt_Tbr)
subSvep1 <- data.frame(seq.name,seq.txt_Svep1)
subSreb2 <- data.frame(seq.name,seq.txt_Sreb2)
subSh3px3 <- data.frame(seq.name,seq.txt_Sh3px3)
subRipk4 <- data.frame(seq.name,seq.txt_Ripk4)
subPhodopsin <- data.frame(seq.name,seq.txt_Phodopsin) 
subRag2 <- data.frame(seq.name,seq.txt_Rag2)
subPtr <- data.frame(seq.name,seq.txt_Ptr)
subPlagl2 <- data.frame(seq.name,seq.txt_Plagl2)
subPanx2 <- data.frame(seq.name,seq.txt_Panx2)
subND4 <- data.frame(seq.name,seq.txt_ND4)
subMyh6 <- data.frame(seq.name,seq.txt_Myh6)
subKiaa1239 <- data.frame(seq.name,seq.txt_Kiaa1239) 
subHoxc6a <- data.frame(seq.name,seq.txt_Hoxc6a)
subGlyt <- data.frame(seq.name,seq.txt_Glyt)
subFicd <- data.frame(seq.name,seq.txt_Ficd)

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
multigene_4851 <- list_df %>% reduce(inner_join, by="seq.name")
#assign column names 
colnames(multigene_4851) <- c("species_name","fasta.12S","fasta.16S","fasta.4c4","fasta.COI","fasta.CytB","fasta.Enc1","fasta.Ficd","fasta.Glyt","fasta.Hoxc6a","fasta.Kiaa1239","fasta.Myh6","fasta.Nd2","fasta.ND4","fasta.Panx2","fasta.Plagl2","fasta.Ptr","fasta.Rag1","fasta.Rag2","fasta.Rhodopsin","fasta.Ripk4","fasta.Sh3px3","fasta.Sidkey","fasta.Sreb2","fasta.Svep1","fasta.Tbr","fasta.Vcpip","fasta.Zic1")
# Concatanate the 27 genes of 8 outgroups in to a one column to make the sequence a multigene one
list_df_outgroup = list(sub12S,sub16S,sub4c4,subCOI,subCytB,subEnc1,subFicd,subGlyt,subHoxc6a,subKiaa1239,subMyh6,subNd2,subND4,subPanx2,subPlagl2,subPtr,subRag1,subRag2,subPhodopsin,subRipk4,subSh3px3,subSidkey,subSreb2,subSvep1,subTbr,subvcpip,subZic1)
#multigene_4856 <- list_df %>% reduce(merge, by="seq.name")
multigene_outgroup <- list_df_outgroup %>% reduce(inner_join, by="seq.name")
#assign column names 
colnames(multigene_outgroup) <- c("species_name","fasta.12S","fasta.16S","fasta.4c4","fasta.COI","fasta.CytB","fasta.Enc1","fasta.Ficd","fasta.Glyt","fasta.Hoxc6a","fasta.Kiaa1239","fasta.Myh6","fasta.Nd2","fasta.ND4","fasta.Panx2","fasta.Plagl2","fasta.Ptr","fasta.Rag1","fasta.Rag2","fasta.Rhodopsin","fasta.Ripk4","fasta.Sh3px3","fasta.Sidkey","fasta.Sreb2","fasta.Svep1","fasta.Tbr","fasta.Vcpip","fasta.Zic1")

multigene_4859 <- rbind(multigene_4851,multigene_outgroup)
#multigene_4856$Seq <- multigene_4856 %>% paste(fasta.12S,fasta.16S,fasta.4c4,fasta.COI,fasta.CytB,fasta.Enc1,fasta.Ficd,fasta.Glyt,fasta.Hoxc6a,fasta.Kiaa1239,fasta.Myh6,fasta.Nd2,fasta.ND4,fasta.Panx2,fasta.Plagl2,fasta.Ptr,fasta.Rag1,fasta.Rag2,fasta.Rhodopsin,fasta.Ripk4,fasta.Sh3px3,fasta.Sidkey,fasta.Sreb2,fasta.Svep1,fasta.Tbr,fasta.Vcpip,fasta.Zic1)

library(tidyr)
multigene_4859 <- multigene_4859 %>% 
  unite("Merged",fasta.12S:fasta.Zic1)
#multigene_4859$Merged

#remove _ in sequences of the multigene as they appear when concatanating 27 genes
# sapply(multigene_4856$Merged,nchar)
multigene_4859_new <- multigene_4859

multigene_4859_new$Merged <- gsub('_','',multigene_4859_new$Merged)
#multigene_4859_new$Merged
sapply(multigene_4859_new$Merged,nchar)

library(readr)
#write the dataframe to a csv file
write_csv(multigene_4859_new,path = "~/Desktop/FishData/New_dataset/Realignments/multigene_4859_new.csv")

#Now remove the species with issues
#Read the csv file with species to be removed
Species_wt_issues <- read_csv("Species_wt_issues.csv")

class(Species_wt_issues)
Species_wt_issues <- as.data.frame(Species_wt_issues)
length(unique(Species_wt_issues$Species_name))#334
colnames(multigene_4859_new) <- c("Species_name","Seq")
Realigned_multigene <- anti_join(multigene_4859_new,Species_wt_issues,by="Species_name")#4526
names(multigene_4859_new)
class(Realigned_multigene)

library(stringr)
library(dplyr)

data1 <- Realigned_multigene[str_detect(Realigned_multigene$Seq, "z"), ] 
#As I found a bad base (z) for Dascyllus_aruanus, I remove that species from the dataset
Realigned_multigene <- Realigned_multigene[-3211,]
#Jinzhong asked to use just 3-5 outgroups, so I am using 4 species 2 from tetrapods and 2 from sharks for the down stream analysis.
#So, now my dataset contains 4517 ray fin fishes and 5 out groups sequences. All together 4521
Realigned_multigene <- Realigned_multigene[-4522:-4524,]
Realigned_multigene <- Realigned_multigene[-4522,]

#Keep a csv file of Realigned_multigene
write_csv(Realigned_multigene,path = "~/Desktop/FishData/New_dataset/Realignments/Fish_Realigned_multigene.csv")
#Keep a fasta file of Realigned_multigene
library(seqRFLP)
dataframe2fas(Realigned_multigene, file = "Fish_Realigned_multigene.fasta")
#Keep a phylip file of Realigned_multigene
dat2phylip(Realigned_multigene,outfile= "Fish_Realigned_multigene.phy")
getwd()

sapply(Realigned_multigene,nchar)#24855
