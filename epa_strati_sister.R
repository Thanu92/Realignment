#EPAng stratified sister position (EPAng using raxml info file as the model)
#Last updated Jan 19, 2024
#Jan 15, 2024
library(ape)
#Sister species table
getwd()
#/home/thanu/Desktop/FishData/New_dataset/Realignments/mafft-linux64"
library(foreach)
#install.packages("phytools")
library(phytools)
library(maps)
library(phangorn)
library(phylotools)
#Read fasta files of placed co1 on the bb tree# in stratified folder
R_co1_20_1=read.fasta(file = "CO1_20_1S.fasta")
R_co1_20_2=read.fasta(file = "CO1_20_2S.fasta")
R_co1_20_3=read.fasta(file = "CO1_20_3S.fasta")
R_co1_20_4=read.fasta(file = "CO1_20_4S.fasta")
R_co1_20_5=read.fasta(file = "CO1_20_5S.fasta")
R_co1_20_6=read.fasta(file = "CO1_20_6S.fasta")
R_co1_20_7=read.fasta(file = "CO1_20_7S.fasta")
R_co1_20_8=read.fasta(file = "CO1_20_8S.fasta")
R_co1_20_9=read.fasta(file = "CO1_20_9S.fasta")
R_co1_20_10=read.fasta(file = "CO1_20_10S.fasta")
R_co1_40_1=read.fasta(file = "CO1_40_1S.fasta")
R_co1_40_2=read.fasta(file = "CO1_40_2S.fasta")
R_co1_40_3=read.fasta(file = "CO1_40_3S.fasta")
R_co1_40_4=read.fasta(file = "CO1_40_4S.fasta")
R_co1_40_5=read.fasta(file = "CO1_40_4S.fasta")
R_co1_40_6=read.fasta(file = "CO1_40_6S.fasta")
R_co1_40_7=read.fasta(file = "CO1_40_7S.fasta")
R_co1_40_8=read.fasta(file = "CO1_40_8S.fasta")
R_co1_40_9=read.fasta(file = "CO1_40_9S.fasta")
R_co1_40_10=read.fasta(file = "CO1_40_10S.fasta")
R_co1_60_1=read.fasta(file = "CO1_60_1S.fasta")
R_co1_60_2=read.fasta(file = "CO1_60_2S.fasta")
R_co1_60_3=read.fasta(file = "CO1_60_2S.fasta")
R_co1_60_4=read.fasta(file = "CO1_60_1S.fasta")
R_co1_60_5=read.fasta(file = "CO1_60_7S.fasta")
R_co1_60_6=read.fasta(file = "CO1_60_8S.fasta")
R_co1_60_7=read.fasta(file = "CO1_60_7S.fasta")
R_co1_60_8=read.fasta(file = "CO1_60_8S.fasta")
R_co1_60_9=read.fasta(file = "CO1_60_9S.fasta")
R_co1_60_10=read.fasta(file = "CO1_60_9S.fasta")
R_co1_80_1=read.fasta(file = "CO1_80_1S.fasta")
R_co1_80_2=read.fasta(file = "CO1_80_2S.fasta")
R_co1_80_3=read.fasta(file = "CO1_80_2S.fasta")
R_co1_80_4=read.fasta(file = "CO1_80_4S.fasta")
R_co1_80_5=read.fasta(file = "CO1_80_7S.fasta")
R_co1_80_6=read.fasta(file = "CO1_80_7S.fasta")
R_co1_80_7=read.fasta(file = "CO1_80_7S.fasta")
R_co1_80_8=read.fasta(file = "CO1_80_8S.fasta")
R_co1_80_9=read.fasta(file = "CO1_80_8S.fasta")
R_co1_80_10=read.fasta(file = "CO1_80_8S.fasta")
R_co1_1_1=read.fasta(file = "CO1_99_2S.fasta")
R_co1_1_2=read.fasta(file = "CO1_99_3S.fasta")
R_co1_1_3=read.fasta(file = "CO1_99_4S.fasta")
R_co1_1_4=read.fasta(file = "CO1_99_5S.fasta")
R_co1_1_5=read.fasta(file = "CO1_99_6S.fasta")
R_co1_1_6=read.fasta(file = "CO1_99_7S.fasta")
R_co1_1_7=read.fasta(file = "CO1_99_8S.fasta")
R_co1_1_8=read.fasta(file = "CO1_99_9S.fasta")
R_co1_1_9=read.fasta(file = "CO1_99_10S.fasta")
R_co1_1_10=read.fasta(file = "CO1_99_1S.fasta")

getwd()

#read the refernce tree
ReferenceTree <- read.tree("RAxML_bestTree.FR100_new")
#Read the 50 trees generated using stratified samples#in epa_strati folder
R20_1_tree<-read.tree("EPA20_1S_new.nwk")
R20_2_tree<-read.tree("EPA20_2S_new.nwk")
R20_3_tree<-read.tree("EPA20_3S_new.nwk")
R20_4_tree<-read.tree("EPA20_4S_new.nwk")
R20_5_tree<-read.tree("EPA20_5S_new.nwk")
R20_6_tree<-read.tree("EPA20_6S_new.nwk")
R20_7_tree<-read.tree("EPA20_7S_new.nwk")
R20_8_tree<-read.tree("EPA20_8S_new.nwk")
R20_9_tree<-read.tree("EPA20_9S_new.nwk")
R20_10_tree<-read.tree("EPA20_10S_new.nwk")
R40_1_tree<-read.tree("EPA40_1S_new.nwk")
R40_2_tree<-read.tree("EPA40_2S_new.nwk")
R40_3_tree<-read.tree("EPA40_3S_new.nwk")
R40_4_tree<-read.tree("EPA40_4S_new.nwk")
R40_5_tree<-read.tree("EPA40_5S_new.nwk")
R40_6_tree<-read.tree("EPA40_6S_new.nwk")
R40_7_tree<-read.tree("EPA40_7S_new.nwk")
R40_8_tree<-read.tree("EPA40_8S_new.nwk")
R40_9_tree<-read.tree("EPA40_9S_new.nwk")
R40_10_tree<-read.tree("EPA40_10S_new.nwk")
R60_1_tree<-read.tree("EPA60_1S_new.nwk")
R60_2_tree<-read.tree("EPA60_2S_new.nwk")
R60_3_tree<-read.tree("EPA60_3S_new.nwk")
R60_4_tree<-read.tree("EPA60_4S_new.nwk")
R60_5_tree<-read.tree("EPA60_5S_new.nwk")
R60_6_tree<-read.tree("EPA60_6S_new.nwk")
R60_7_tree<-read.tree("EPA60_7S_new.nwk")
R60_8_tree<-read.tree("EPA60_8S_new.nwk")
R60_9_tree<-read.tree("EPA60_9S_new.nwk")
R60_10_tree<-read.tree("EPA60_10S_new.nwk")
R80_1_tree<-read.tree("EPA80_1S_new.nwk")
R80_2_tree<-read.tree("EPA80_2S_new.nwk")
R80_3_tree<-read.tree("EPA80_3S_new.nwk")
R80_4_tree<-read.tree("EPA80_4S_new.nwk")
R80_5_tree<-read.tree("EPA80_5S_new.nwk")
R80_6_tree<-read.tree("EPA80_6S_new.nwk")
R80_7_tree<-read.tree("EPA80_7S_new.nwk")
R80_8_tree<-read.tree("EPA80_8S_new.nwk")
R80_9_tree<-read.tree("EPA80_9S_new.nwk")
R80_10_tree<-read.tree("EPA80_10S_new.nwk")
R80_1_tree<-read.tree("EPA80_1S_new.nwk")
class(R80_1_tree)
R99_1_tree<-read.tree("EPA99_1S_new.nwk")
R99_2_tree<-read.tree("EPA99_2S_new.nwk")
R99_3_tree<-read.tree("EPA99_3S_new.nwk")
R99_4_tree<-read.tree("EPA99_4S_new.nwk")
R99_5_tree<-read.tree("EPA99_5S_new.nwk")
R99_6_tree<-read.tree("EPA99_6S_new.nwk")
R99_7_tree<-read.tree("EPA99_7S_new.nwk")
R99_8_tree<-read.tree("EPA99_8S_new.nwk")
R99_9_tree<-read.tree("EPA99_9S_new.nwk")
R99_10_tree<-read.tree("EPA99_10S_new.nwk")
#------
names(R_co1_20_1)
R_co1_20_1_species <- R_co1_20_1$seq.name
class(R_co1_20_1_species )
R_co1_20_2_species <- R_co1_20_2$seq.name
R_co1_20_3_species <- R_co1_20_3$seq.name
R_co1_20_4_species <- R_co1_20_4$seq.name
R_co1_20_5_species <- R_co1_20_5$seq.name
R_co1_20_6_species <- R_co1_20_6$seq.name
R_co1_20_7_species <- R_co1_20_7$seq.name
R_co1_20_8_species <- R_co1_20_8$seq.name
R_co1_20_9_species <- R_co1_20_9$seq.name
R_co1_20_10_species <- R_co1_20_10$seq.name
R_co1_40_1_species <- R_co1_40_1$seq.name
R_co1_40_2_species <- R_co1_40_2$seq.name
R_co1_40_3_species <- R_co1_40_3$seq.name
R_co1_40_4_species <- R_co1_40_4$seq.name
R_co1_40_5_species <- R_co1_40_5$seq.name
R_co1_40_6_species <- R_co1_40_6$seq.name
R_co1_40_7_species <- R_co1_40_7$seq.name
R_co1_40_8_species <- R_co1_40_8$seq.name
R_co1_40_9_species <- R_co1_40_9$seq.name
R_co1_40_10_species <- R_co1_40_10$seq.name
R_co1_60_1_species <- R_co1_60_1$seq.name
R_co1_60_2_species <- R_co1_60_2$seq.name
R_co1_60_3_species <- R_co1_60_3$seq.name
R_co1_60_4_species <- R_co1_60_4$seq.name
R_co1_60_5_species <- R_co1_60_5$seq.name
R_co1_60_6_species <- R_co1_60_6$seq.name
R_co1_60_7_species <- R_co1_60_7$seq.name
R_co1_60_8_species <- R_co1_60_8$seq.name
R_co1_60_9_species <- R_co1_60_9$seq.name
R_co1_60_10_species <- R_co1_60_10$seq.name
R_co1_80_1_species <- R_co1_80_1$seq.name
R_co1_80_2_species <- R_co1_80_2$seq.name
R_co1_80_3_species <- R_co1_80_3$seq.name
R_co1_80_4_species <- R_co1_80_4$seq.name
R_co1_80_5_species <- R_co1_80_5$seq.name
R_co1_80_6_species <- R_co1_80_6$seq.name
R_co1_80_7_species <- R_co1_80_7$seq.name
R_co1_80_8_species <- R_co1_80_8$seq.name
R_co1_80_9_species <- R_co1_80_9$seq.name
R_co1_80_10_species <- R_co1_80_10$seq.name
R_co1_1_1_species <- R_co1_1_1$seq.name
R_co1_1_2_species <- R_co1_1_2$seq.name
R_co1_1_3_species <- R_co1_1_3$seq.name
R_co1_1_4_species <- R_co1_1_4$seq.name
R_co1_1_5_species <- R_co1_1_5$seq.name
R_co1_1_6_species <- R_co1_1_6$seq.name
R_co1_1_7_species <- R_co1_1_7$seq.name
R_co1_1_8_species <- R_co1_1_8$seq.name
R_co1_1_9_species <- R_co1_1_9$seq.name
R_co1_1_10_species <- R_co1_1_10$seq.name



#--------------------------
# R80_1_tree <- as.phylo(R80_1_tree)
# class(R80_1_tree)
##bb 20% tree when placed 80% co1
library(ips)
#sister(tree, "t4", type = "terminal", label = T)
RefTree_R_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(ReferenceTree,R_co1_20_1_species[i],type="terminal",label=T)
R_80_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(R80_1_tree,R_co1_20_1_species[i],type="terminal",label=T) #Error in sister(R80_1_tree, R_co1_20_1_species[i], type = "terminal",  : 
#task 7 failed - "argument is of length zero"

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_1_List [(RefTree_R_co1_20_1_List %in% R_80_co1_20_1_List )] #2018
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_1_List [!(RefTree_R_co1_20_1_List %in% R_80_co1_20_1_List)] #1265

RefTree_R_co1_20_2_List <- foreach(i=1:length(R_co1_20_2_species)) %do% sister(ReferenceTree,R_co1_20_2_species[i],type="terminal",label=T)
R_80_co1_20_2_List <- foreach(i=1:length(R_co1_20_2_species)) %do% sister(R80_2_tree,R_co1_20_2_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_2_List [(RefTree_R_co1_20_2_List %in% R_80_co1_20_2_List )] #1939
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_2_List [!(RefTree_R_co1_20_2_List %in% R_80_co1_20_2_List)] #1344

RefTree_R_co1_20_3_List <- foreach(i=1:length(R_co1_20_3_species)) %do% sister(ReferenceTree,R_co1_20_3_species[i],type="terminal",label=T)
R_80_co1_20_3_List <- foreach(i=1:length(R_co1_20_3_species)) %do% sister(R80_3_tree,R_co1_20_3_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_3_List [(RefTree_R_co1_20_3_List %in% R_80_co1_20_3_List )] #1971
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_3_List [!(RefTree_R_co1_20_3_List %in% R_80_co1_20_3_List)] #1312

RefTree_R_co1_20_4_List <- foreach(i=1:length(R_co1_20_4_species)) %do% sister(ReferenceTree,R_co1_20_4_species[i],type="terminal",label=T)
R_80_co1_20_4_List <- foreach(i=1:length(R_co1_20_4_species)) %do% sister(R80_4_tree,R_co1_20_4_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_4_List [(RefTree_R_co1_20_4_List %in% R_80_co1_20_4_List )] #2016
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_4_List [!(RefTree_R_co1_20_4_List %in% R_80_co1_20_4_List)] #1267

RefTree_R_co1_20_5_List <- foreach(i=1:length(R_co1_20_5_species)) %do% sister(ReferenceTree,R_co1_20_5_species[i],type="terminal",label=T)
R_80_co1_20_5_List <- foreach(i=1:length(R_co1_20_5_species)) %do% sister(R80_5_tree,R_co1_20_5_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_5_List [(RefTree_R_co1_20_5_List %in% R_80_co1_20_5_List )] #1973
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_5_List [!(RefTree_R_co1_20_5_List %in% R_80_co1_20_5_List)] #1310

RefTree_R_co1_20_6_List <- foreach(i=1:length(R_co1_20_6_species)) %do% sister(ReferenceTree,R_co1_20_6_species[i],type="terminal",label=T)
R_80_co1_20_6_List <- foreach(i=1:length(R_co1_20_6_species)) %do% sister(R80_6_tree,R_co1_20_6_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_6_List [(RefTree_R_co1_20_6_List %in% R_80_co1_20_6_List )] #1910
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_6_List [!(RefTree_R_co1_20_6_List %in% R_80_co1_20_6_List)] #1373

RefTree_R_co1_20_7_List <- foreach(i=1:length(R_co1_20_7_species)) %do% sister(ReferenceTree,R_co1_20_7_species[i],type="terminal",label=T)
R_80_co1_20_7_List <- foreach(i=1:length(R_co1_20_7_species)) %do% sister(R80_7_tree,R_co1_20_7_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_7_List [(RefTree_R_co1_20_7_List %in% R_80_co1_20_7_List )] #1975
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_7_List [!(RefTree_R_co1_20_7_List %in% R_80_co1_20_7_List)] #1308

RefTree_R_co1_20_8_List <- foreach(i=1:length(R_co1_20_8_species)) %do% sister(ReferenceTree,R_co1_20_8_species[i],type="terminal",label=T)
R_80_co1_20_8_List <- foreach(i=1:length(R_co1_20_8_species)) %do% sister(R80_8_tree,R_co1_20_8_species[i],type="terminal",label=T)



#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_8_List [(RefTree_R_co1_20_8_List %in% R_80_co1_20_8_List )] #1999
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_8_List [!(RefTree_R_co1_20_8_List %in% R_80_co1_20_8_List)] #1284

RefTree_R_co1_20_9_List <- foreach(i=1:length(R_co1_20_9_species)) %do% sister(ReferenceTree,R_co1_20_9_species[i],type="terminal",label=T)
R_80_co1_20_9_List <- foreach(i=1:length(R_co1_20_9_species)) %do% sister(R80_9_tree,R_co1_20_9_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_9_List [(RefTree_R_co1_20_9_List %in% R_80_co1_20_9_List )] #1964
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_9_List [!(RefTree_R_co1_20_9_List %in% R_80_co1_20_9_List)] #1319

RefTree_R_co1_20_10_List <- foreach(i=1:length(R_co1_20_10_species)) %do% sister(ReferenceTree,R_co1_20_10_species[i],type="terminal",label=T)
R_80_co1_20_10_List <- foreach(i=1:length(R_co1_20_10_species)) %do% sister(R80_10_tree,R_co1_20_10_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_10_List [(RefTree_R_co1_20_10_List %in% R_80_co1_20_10_List )] #1993
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_10_List [!(RefTree_R_co1_20_10_List %in% R_80_co1_20_10_List)] #1290

#bb 60% tree when placed 40% co1

RefTree_R_co1_40_1_List <- foreach(i=1:length(R_co1_40_1_species)) %do% sister(ReferenceTree,R_co1_40_1_species[i],type="terminal",label=T)
R_60_co1_40_1_List <- foreach(i=1:length(R_co1_40_1_species)) %do% sister(R60_1_tree,R_co1_40_1_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_1_List [(RefTree_R_co1_40_1_List %in% R_60_co1_40_1_List )] #1061
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_1_List [!(RefTree_R_co1_40_1_List %in% R_60_co1_40_1_List )] #1408

RefTree_R_co1_40_2_List <- foreach(i=1:length(R_co1_40_2_species)) %do% sister(ReferenceTree,R_co1_40_2_species[i],type="terminal",label=T)
R_60_co1_40_2_List <- foreach(i=1:length(R_co1_40_2_species)) %do% sister(R60_2_tree,R_co1_40_2_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_2_List [(RefTree_R_co1_40_2_List %in% R_60_co1_40_2_List )] #1067
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_2_List [!(RefTree_R_co1_40_2_List %in% R_60_co1_40_2_List )] #1402

RefTree_R_co1_40_3_List <- foreach(i=1:length(R_co1_40_3_species)) %do% sister(ReferenceTree,R_co1_40_3_species[i],type="terminal",label=T)
R_60_co1_40_3_List <- foreach(i=1:length(R_co1_40_3_species)) %do% sister(R60_3_tree,R_co1_40_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_3_List [(RefTree_R_co1_40_3_List %in% R_60_co1_40_3_List )] #1077
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_3_List [!(RefTree_R_co1_40_3_List %in% R_60_co1_40_3_List )] #1392

RefTree_R_co1_40_4_List <- foreach(i=1:length(R_co1_40_4_species)) %do% sister(ReferenceTree,R_co1_40_4_species[i],type="terminal",label=T)
R_60_co1_40_4_List <- foreach(i=1:length(R_co1_40_4_species)) %do% sister(R60_4_tree,R_co1_40_4_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_4_List [(RefTree_R_co1_40_4_List %in% R_60_co1_40_4_List )] #1065
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_4_List [!(RefTree_R_co1_40_4_List %in% R_60_co1_40_4_List )] #1404

RefTree_R_co1_40_5_List <- foreach(i=1:length(R_co1_40_5_species)) %do% sister(ReferenceTree,R_co1_40_5_species[i],type="terminal",label=T)
R_60_co1_40_5_List <- foreach(i=1:length(R_co1_40_5_species)) %do% sister(R60_5_tree,R_co1_40_5_species[i],type="terminal",label=T)
#
#
# #get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_5_List [(RefTree_R_co1_40_5_List %in% R_60_co1_40_5_List )] #1158
# #get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_5_List [!(RefTree_R_co1_40_5_List %in% R_60_co1_40_5_List )] #1311

RefTree_R_co1_40_6_List <- foreach(i=1:length(R_co1_40_6_species)) %do% sister(ReferenceTree,R_co1_40_6_species[i],type="terminal",label=T)
R_60_co1_40_6_List <- foreach(i=1:length(R_co1_40_6_species)) %do% sister(R60_6_tree,R_co1_40_6_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_6_List [(RefTree_R_co1_40_6_List %in% R_60_co1_40_6_List )] #1050
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_6_List [!(RefTree_R_co1_40_6_List %in% R_60_co1_40_6_List )] #1419

RefTree_R_co1_40_7_List <- foreach(i=1:length(R_co1_40_7_species)) %do% sister(ReferenceTree,R_co1_40_7_species[i],type="terminal",label=T)
R_60_co1_40_7_List <- foreach(i=1:length(R_co1_40_7_species)) %do% sister(R60_7_tree,R_co1_40_7_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_7_List [(RefTree_R_co1_40_7_List %in% R_60_co1_40_7_List )] #981
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_7_List [!(RefTree_R_co1_40_7_List %in% R_60_co1_40_7_List )] #1488

RefTree_R_co1_40_8_List <- foreach(i=1:length(R_co1_40_8_species)) %do% sister(ReferenceTree,R_co1_40_8_species[i],type="terminal",label=T)
R_60_co1_40_8_List <- foreach(i=1:length(R_co1_40_8_species)) %do% sister(R60_8_tree,R_co1_40_8_species[i],type="terminal",label=T)
#sister(R60_8_tree,"Synodontis_serratus",type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_8_List [(RefTree_R_co1_40_8_List %in% R_60_co1_40_8_List )] #1067
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_8_List [!(RefTree_R_co1_40_8_List %in% R_60_co1_40_8_List )] #1402

RefTree_R_co1_40_9_List <- foreach(i=1:length(R_co1_40_9_species)) %do% sister(ReferenceTree,R_co1_40_9_species[i],type="terminal",label=T)
R_60_co1_40_9_List <- foreach(i=1:length(R_co1_40_9_species)) %do% sister(R60_9_tree,R_co1_40_9_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_9_List [(RefTree_R_co1_40_9_List %in% R_60_co1_40_9_List )] #986
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_9_List [!(RefTree_R_co1_40_9_List %in% R_60_co1_40_9_List )] #1483

RefTree_R_co1_40_10_List <- foreach(i=1:length(R_co1_40_10_species)) %do% sister(ReferenceTree,R_co1_40_10_species[i],type="terminal",label=T)
R_60_co1_40_10_List <- foreach(i=1:length(R_co1_40_10_species)) %do% sister(R60_10_tree,R_co1_40_10_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_10_List [(RefTree_R_co1_40_10_List %in% R_60_co1_40_10_List )] #1038
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_10_List [!(RefTree_R_co1_40_10_List %in% R_60_co1_40_10_List )] #1431

#RefTree_biasco1List[sapply(names(RefTree_biasco1List), function(x) !identical(RefTree_biasco1List[[x]], biasco1List[[x]]))] 

#bb 40% tree when placed 60% co1

RefTree_R_co1_60_1_List <- foreach(i=1:length(R_co1_60_1_species)) %do% sister(ReferenceTree,R_co1_60_1_species[i],type="terminal",label=T)
R_40_co1_60_1_List <- foreach(i=1:length(R_co1_60_1_species)) %do% sister(R40_1_tree,R_co1_60_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_1_List [(RefTree_R_co1_60_1_List %in% R_40_co1_60_1_List )] #393
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_1_List [!(RefTree_R_co1_60_1_List %in% R_40_co1_60_1_List )] #1261

RefTree_R_co1_60_2_List <- foreach(i=1:length(R_co1_60_2_species)) %do% sister(ReferenceTree,R_co1_60_2_species[i],type="terminal",label=T)
R_40_co1_60_2_List <- foreach(i=1:length(R_co1_60_2_species)) %do% sister(R40_2_tree,R_co1_60_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_2_List [(RefTree_R_co1_60_2_List %in% R_40_co1_60_2_List )] #383
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_2_List [!(RefTree_R_co1_60_2_List %in% R_40_co1_60_2_List )] #1271


RefTree_R_co1_60_3_List <- foreach(i=1:length(R_co1_60_3_species)) %do% sister(ReferenceTree,R_co1_60_3_species[i],type="terminal",label=T)
R_40_co1_60_3_List <- foreach(i=1:length(R_co1_60_3_species)) %do% sister(R40_3_tree,R_co1_60_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_3_List [(RefTree_R_co1_60_3_List %in% R_40_co1_60_3_List )] #583
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_3_List [!(RefTree_R_co1_60_3_List %in% R_40_co1_60_3_List )] #1071

RefTree_R_co1_60_4_List <- foreach(i=1:length(R_co1_60_4_species)) %do% sister(ReferenceTree,R_co1_60_4_species[i],type="terminal",label=T)
R_40_co1_60_4_List <- foreach(i=1:length(R_co1_60_4_species)) %do% sister(R40_4_tree,R_co1_60_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_4_List [(RefTree_R_co1_60_4_List %in% R_40_co1_60_4_List )] #490
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_4_List [!(RefTree_R_co1_60_4_List %in% R_40_co1_60_4_List )] #1164

RefTree_R_co1_60_5_List <- foreach(i=1:length(R_co1_60_5_species)) %do% sister(ReferenceTree,R_co1_60_5_species[i],type="terminal",label=T)
R_40_co1_60_5_List <- foreach(i=1:length(R_co1_60_5_species)) %do% sister(R40_5_tree,R_co1_60_5_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_5_List [(RefTree_R_co1_60_5_List %in% R_40_co1_60_5_List )] #553
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_5_List [!(RefTree_R_co1_60_5_List %in% R_40_co1_60_5_List )] #1101

RefTree_R_co1_60_6_List <- foreach(i=1:length(R_co1_60_6_species)) %do% sister(ReferenceTree,R_co1_60_6_species[i],type="terminal",label=T)
R_40_co1_60_6_List <- foreach(i=1:length(R_co1_60_6_species)) %do% sister(R40_6_tree,R_co1_60_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_6_List [(RefTree_R_co1_60_6_List %in% R_40_co1_60_6_List )] #561
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_6_List [!(RefTree_R_co1_60_6_List %in% R_40_co1_60_6_List )] #1093

RefTree_R_co1_60_7_List <- foreach(i=1:length(R_co1_60_7_species)) %do% sister(ReferenceTree,R_co1_60_7_species[i],type="terminal",label=T)
R_40_co1_60_7_List <- foreach(i=1:length(R_co1_60_7_species)) %do% sister(R40_7_tree,R_co1_60_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_7_List [(RefTree_R_co1_60_7_List %in% R_40_co1_60_7_List )] #383
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_7_List [!(RefTree_R_co1_60_7_List %in% R_40_co1_60_7_List )] #1271

RefTree_R_co1_60_8_List <- foreach(i=1:length(R_co1_60_8_species)) %do% sister(ReferenceTree,R_co1_60_8_species[i],type="terminal",label=T)
R_40_co1_60_8_List <- foreach(i=1:length(R_co1_60_8_species)) %do% sister(R40_8_tree,R_co1_60_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_8_List [(RefTree_R_co1_60_8_List %in% R_40_co1_60_8_List )] #387
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_8_List [!(RefTree_R_co1_60_8_List %in% R_40_co1_60_8_List )] #1267

RefTree_R_co1_60_9_List <- foreach(i=1:length(R_co1_60_9_species)) %do% sister(ReferenceTree,R_co1_60_9_species[i],type="terminal",label=T)
R_40_co1_60_9_List <- foreach(i=1:length(R_co1_60_9_species)) %do% sister(R40_9_tree,R_co1_60_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_9_List [(RefTree_R_co1_60_9_List %in% R_40_co1_60_9_List )] #384
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_9_List [!(RefTree_R_co1_60_9_List %in% R_40_co1_60_9_List )] #1270

RefTree_R_co1_60_10_List <- foreach(i=1:length(R_co1_60_10_species)) %do% sister(ReferenceTree,R_co1_60_10_species[i],type="terminal",label=T)
R_40_co1_60_10_List <- foreach(i=1:length(R_co1_60_10_species)) %do% sister(R40_10_tree,R_co1_60_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_10_List [(RefTree_R_co1_60_10_List %in% R_40_co1_60_10_List )] #565
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_10_List [!(RefTree_R_co1_60_10_List %in% R_40_co1_60_10_List )] #1089


#bb 20% tree when placed 80% co1

RefTree_R_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(ReferenceTree,R_co1_80_1_species[i],type="terminal",label=T)
R_20_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(R20_1_tree,R_co1_80_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_1_List [(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #102
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_1_List [!(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #711

RefTree_R_co1_80_2_List <- foreach(i=1:length(R_co1_80_2_species)) %do% sister(ReferenceTree,R_co1_80_2_species[i],type="terminal",label=T)
R_20_co1_80_2_List <- foreach(i=1:length(R_co1_80_2_species)) %do% sister(R20_2_tree,R_co1_80_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_2_List [(RefTree_R_co1_80_2_List %in% R_20_co1_80_2_List )] #93
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_2_List [!(RefTree_R_co1_80_2_List %in% R_20_co1_80_2_List )] #720

RefTree_R_co1_80_3_List <- foreach(i=1:length(R_co1_80_3_species)) %do% sister(ReferenceTree,R_co1_80_3_species[i],type="terminal",label=T)
R_20_co1_80_3_List <- foreach(i=1:length(R_co1_80_3_species)) %do% sister(R20_3_tree,R_co1_80_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_3_List [(RefTree_R_co1_80_3_List %in% R_20_co1_80_3_List )] #146
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_3_List [!(RefTree_R_co1_80_3_List %in% R_20_co1_80_3_List )] #667

RefTree_R_co1_80_4_List <- foreach(i=1:length(R_co1_80_4_species)) %do% sister(ReferenceTree,R_co1_80_4_species[i],type="terminal",label=T)
R_20_co1_80_4_List <- foreach(i=1:length(R_co1_80_4_species)) %do% sister(R20_4_tree,R_co1_80_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_4_List [(RefTree_R_co1_80_4_List %in% R_20_co1_80_4_List )] #90
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_4_List [!(RefTree_R_co1_80_4_List %in% R_20_co1_80_4_List )] #723

RefTree_R_co1_80_5_List <- foreach(i=1:length(R_co1_80_5_species)) %do% sister(ReferenceTree,R_co1_80_5_species[i],type="terminal",label=T)
R_20_co1_80_5_List <- foreach(i=1:length(R_co1_80_5_species)) %do% sister(R20_5_tree,R_co1_80_5_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_5_List [(RefTree_R_co1_80_5_List %in% R_20_co1_80_5_List )] #159
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_5_List [!(RefTree_R_co1_80_5_List %in% R_20_co1_80_5_List )] #654

RefTree_R_co1_80_6_List <- foreach(i=1:length(R_co1_80_6_species)) %do% sister(ReferenceTree,R_co1_80_6_species[i],type="terminal",label=T)
R_20_co1_80_6_List <- foreach(i=1:length(R_co1_80_6_species)) %do% sister(R20_6_tree,R_co1_80_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_6_List [(RefTree_R_co1_80_6_List %in% R_20_co1_80_6_List )] #155
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_6_List [!(RefTree_R_co1_80_6_List %in% R_20_co1_80_6_List )] #658

RefTree_R_co1_80_7_List <- foreach(i=1:length(R_co1_80_7_species)) %do% sister(ReferenceTree,R_co1_80_7_species[i],type="terminal",label=T)
R_20_co1_80_7_List <- foreach(i=1:length(R_co1_80_7_species)) %do% sister(R20_7_tree,R_co1_80_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_7_List [(RefTree_R_co1_80_7_List %in% R_20_co1_80_7_List )] #83
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_7_List [!(RefTree_R_co1_80_7_List %in% R_20_co1_80_7_List )] #730

RefTree_R_co1_80_8_List <- foreach(i=1:length(R_co1_80_8_species)) %do% sister(ReferenceTree,R_co1_80_8_species[i],type="terminal",label=T)
R_20_co1_80_8_List <- foreach(i=1:length(R_co1_80_8_species)) %do% sister(R20_8_tree,R_co1_80_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_8_List [(RefTree_R_co1_80_8_List %in% R_20_co1_80_8_List )] #92
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_8_List [!(RefTree_R_co1_80_8_List %in% R_20_co1_80_8_List )] #721

RefTree_R_co1_80_9_List <- foreach(i=1:length(R_co1_80_9_species)) %do% sister(ReferenceTree,R_co1_80_9_species[i],type="terminal",label=T)
R_20_co1_80_9_List <- foreach(i=1:length(R_co1_80_9_species)) %do% sister(R20_9_tree,R_co1_80_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_9_List [(RefTree_R_co1_80_9_List %in% R_20_co1_80_9_List )] #140
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_9_List [!(RefTree_R_co1_80_9_List %in% R_20_co1_80_9_List )] #673

RefTree_R_co1_80_10_List <- foreach(i=1:length(R_co1_80_10_species)) %do% sister(ReferenceTree,R_co1_80_10_species[i],type="terminal",label=T)
R_20_co1_80_10_List <- foreach(i=1:length(R_co1_80_10_species)) %do% sister(R20_10_tree,R_co1_80_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_10_List [(RefTree_R_co1_80_10_List %in% R_20_co1_80_10_List )] #132
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_10_List [!(RefTree_R_co1_80_10_List %in% R_20_co1_80_10_List )] #681
#---------------------


#bb 99% tree when placed 1% co1
RefTree_R_co1_1_1_List <- foreach(i=1:length(R_co1_1_1_species)) %do% sister(ReferenceTree,R_co1_1_1_species[i],type="terminal",label=T)
R_99_co1_1_1_List <- foreach(i=1:length(R_co1_1_1_species)) %do% sister(R99_1_tree,R_co1_1_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_1_List [(RefTree_R_co1_1_1_List %in% R_99_co1_1_1_List )] #20
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_1_List [!(RefTree_R_co1_1_1_List %in% R_99_co1_1_1_List)] #7

RefTree_R_co1_1_2_List <- foreach(i=1:length(R_co1_1_2_species)) %do% sister(ReferenceTree,R_co1_1_2_species[i],type="terminal",label=T)
R_99_co1_1_2_List <- foreach(i=1:length(R_co1_1_2_species)) %do% sister(R99_2_tree,R_co1_1_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_2_List [(RefTree_R_co1_1_2_List %in% R_99_co1_1_2_List )] #20
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_2_List [!(RefTree_R_co1_1_2_List %in% R_99_co1_1_2_List)] #7

RefTree_R_co1_1_3_List <- foreach(i=1:length(R_co1_1_3_species)) %do% sister(ReferenceTree,R_co1_1_3_species[i],type="terminal",label=T)
R_99_co1_1_3_List <- foreach(i=1:length(R_co1_1_3_species)) %do% sister(R99_3_tree,R_co1_1_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_3_List [(RefTree_R_co1_1_3_List %in% R_99_co1_1_3_List )] #20
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_3_List [!(RefTree_R_co1_1_3_List %in% R_99_co1_1_3_List)] #7

RefTree_R_co1_1_4_List <- foreach(i=1:length(R_co1_1_4_species)) %do% sister(ReferenceTree,R_co1_1_4_species[i],type="terminal",label=T)
R_99_co1_1_4_List <- foreach(i=1:length(R_co1_1_4_species)) %do% sister(R99_4_tree,R_co1_1_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_4_List [(RefTree_R_co1_1_4_List %in% R_99_co1_1_4_List )] #23
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_4_List [!(RefTree_R_co1_1_4_List %in% R_99_co1_1_4_List)] #4

RefTree_R_co1_1_5_List <- foreach(i=1:length(R_co1_1_5_species)) %do% sister(ReferenceTree,R_co1_1_5_species[i],type="terminal",label=T)
R_99_co1_1_5_List <- foreach(i=1:length(R_co1_1_5_species)) %do% sister(R99_5_tree,R_co1_1_5_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_5_List [(RefTree_R_co1_1_5_List %in% R_99_co1_1_5_List )] #16
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_5_List [!(RefTree_R_co1_1_5_List %in% R_99_co1_1_5_List)] #11

RefTree_R_co1_1_6_List <- foreach(i=1:length(R_co1_1_6_species)) %do% sister(ReferenceTree,R_co1_1_6_species[i],type="terminal",label=T)
R_99_co1_1_6_List <- foreach(i=1:length(R_co1_1_6_species)) %do% sister(R99_6_tree,R_co1_1_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_6_List [(RefTree_R_co1_1_6_List %in% R_99_co1_1_6_List )] #18
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_6_List [!(RefTree_R_co1_1_6_List %in% R_99_co1_1_6_List)] #9

RefTree_R_co1_1_7_List <- foreach(i=1:length(R_co1_1_7_species)) %do% sister(ReferenceTree,R_co1_1_7_species[i],type="terminal",label=T)
R_99_co1_1_7_List <- foreach(i=1:length(R_co1_1_7_species)) %do% sister(R99_7_tree,R_co1_1_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_7_List [(RefTree_R_co1_1_7_List %in% R_99_co1_1_7_List )] #19
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_7_List [!(RefTree_R_co1_1_7_List %in% R_99_co1_1_7_List)] #8

RefTree_R_co1_1_8_List <- foreach(i=1:length(R_co1_1_8_species)) %do% sister(ReferenceTree,R_co1_1_8_species[i],type="terminal",label=T)
R_99_co1_1_8_List <- foreach(i=1:length(R_co1_1_8_species)) %do% sister(R99_8_tree,R_co1_1_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_8_List [(RefTree_R_co1_1_8_List %in% R_99_co1_1_8_List )] #19
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_8_List [!(RefTree_R_co1_1_8_List %in% R_99_co1_1_8_List)] #8

RefTree_R_co1_1_9_List <- foreach(i=1:length(R_co1_1_9_species)) %do% sister(ReferenceTree,R_co1_1_9_species[i],type="terminal",label=T)
R_99_co1_1_9_List <- foreach(i=1:length(R_co1_1_9_species)) %do% sister(R99_9_tree,R_co1_1_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_9_List [(RefTree_R_co1_1_9_List %in% R_99_co1_1_9_List )] #22
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_9_List [!(RefTree_R_co1_1_9_List %in% R_99_co1_1_9_List)] #5


RefTree_R_co1_1_10_List <- foreach(i=1:length(R_co1_1_10_species)) %do% sister(ReferenceTree,R_co1_1_10_species[i],type="terminal",label=T)
R_99_co1_1_10_List <- foreach(i=1:length(R_co1_1_10_species)) %do% sister(R99_10_tree,R_co1_1_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_10_List [(RefTree_R_co1_1_10_List %in% R_99_co1_1_10_List )] #21
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_10_List [!(RefTree_R_co1_1_10_List %in% R_99_co1_1_10_List)] #6
