#Aug 1, 2024
#epang_stratified sampling
#working in stratified sampling folder
#epang_stratified
#Sister species table
getwd()
library(foreach)
#install.packages("phytools")
library(phytools)
library(maps)
library(phangorn)
library(phylotools)
library(ips)
library(dplyr)
library(tidyr)
library(stringr)

### Functions  

# Placing my code in one function extract_Genus (takes one input - RefTree_R_co1_20_1_List or R_80_co1_20_1_List)
extract_Genus <- function(list){
  # use a regex expression to grab the first part of every species name before the underscore character and place in a new list
  # The new list will retain the number of list elements in each sublist that existed previously
  foreach(i=1:length(list)) %do% {
    list[[i]] <- str_extract(list[[i]], "[^_]+")
  }
  
  # Iterate through each list and convert to dataframe format using lapply (very fast function)
  dfList <- lapply(list, function(x) { x <- as.data.frame(x) })
  
  # Need the data.table package for this step, can use rbindlist with idcol = TRUE (this will give an id column corresponding to which list element was used, ex: id=1 for list element 1)
  # Rbindlist is very fast and will combine all the dataframes together quickly
  df <- rbindlist(dfList, idcol=TRUE)
  
  # Edit the column names
  # id column
  colnames(df)[1] <- "id"
  # genus column
  colnames(df)[2] <- "genus"
  
  # remove unneeded vars
  rm(dfList, list)
  
  # return dataframe
  return(df)
}


#read 50 placement fasta files
placement_files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/STRATIFIED", pattern="CO1_[0-9]{1,2}_[0-9]{1,2}S.fasta", full.names=TRUE, recursive=FALSE)

placementList <- foreach(i=1:length(placement_files)) %do% read.fasta(placement_files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file

names(placementList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/STRATIFIED/", "", placement_files)

# This will let us keep track of which list element corresponds to which tree file

names(placementList)

# [1] "CO1_20_10S.fasta" "CO1_20_1S.fasta"  "CO1_20_2S.fasta"  "CO1_20_3S.fasta"  "CO1_20_4S.fasta"  "CO1_20_5S.fasta"  "CO1_20_6S.fasta"  "CO1_20_7S.fasta"  "CO1_20_8S.fasta"  "CO1_20_9S.fasta" 
# [11] "CO1_40_10S.fasta" "CO1_40_1S.fasta"  "CO1_40_2S.fasta"  "CO1_40_3S.fasta"  "CO1_40_4S.fasta"  "CO1_40_5S.fasta"  "CO1_40_6S.fasta"  "CO1_40_7S.fasta"  "CO1_40_8S.fasta"  "CO1_40_9S.fasta" 
# [21] "CO1_60_10S.fasta" "CO1_60_1S.fasta"  "CO1_60_2S.fasta"  "CO1_60_3S.fasta"  "CO1_60_4S.fasta"  "CO1_60_5S.fasta"  "CO1_60_6S.fasta"  "CO1_60_7S.fasta"  "CO1_60_8S.fasta"  "CO1_60_9S.fasta" 
# [31] "CO1_80_10S.fasta" "CO1_80_1S.fasta"  "CO1_80_2S.fasta"  "CO1_80_3S.fasta"  "CO1_80_4S.fasta"  "CO1_80_5S.fasta"  "CO1_80_6S.fasta"  "CO1_80_7S.fasta"  "CO1_80_8S.fasta"  "CO1_80_9S.fasta" 
# [41] "CO1_99_10S.fasta" "CO1_99_1S.fasta"  "CO1_99_2S.fasta"  "CO1_99_3S.fasta"  "CO1_99_4S.fasta"  "CO1_99_5S.fasta"  "CO1_99_6S.fasta"  "CO1_99_7S.fasta"  "CO1_99_8S.fasta"  "CO1_99_9S.fasta" # we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:

#placementList$`CO1_20_1S.fasta`

#Read 50 tree files
tree_files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/STRATIFIED", pattern="EPA[0-9]{1,2}_[0-9]{1,2}S_new.nwk", full.names=TRUE, recursive=FALSE)

treeList <- foreach(i=1:length(tree_files)) %do% read.tree(tree_files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file

names(treeList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/STRATIFIED/", "", tree_files)

# This will let us keep track of which list element corresponds to which tree file

names(treeList)
# [1] "EPA20_10S_new.nwk" "EPA20_1S_new.nwk"  "EPA20_2S_new.nwk"  "EPA20_3S_new.nwk"  "EPA20_4S_new.nwk"  "EPA20_5S_new.nwk"  "EPA20_6S_new.nwk"  "EPA20_7S_new.nwk"  "EPA20_8S_new.nwk"  "EPA20_9S_new.nwk" 
# [11] "EPA40_10S_new.nwk" "EPA40_1S_new.nwk"  "EPA40_2S_new.nwk"  "EPA40_3S_new.nwk"  "EPA40_4S_new.nwk"  "EPA40_5S_new.nwk"  "EPA40_6S_new.nwk"  "EPA40_7S_new.nwk"  "EPA40_8S_new.nwk"  "EPA40_9S_new.nwk" 
# [21] "EPA60_10S_new.nwk" "EPA60_1S_new.nwk"  "EPA60_2S_new.nwk"  "EPA60_3S_new.nwk"  "EPA60_4S_new.nwk"  "EPA60_5S_new.nwk"  "EPA60_6S_new.nwk"  "EPA60_7S_new.nwk"  "EPA60_8S_new.nwk"  "EPA60_9S_new.nwk" 
# [31] "EPA80_10S_new.nwk" "EPA80_1S_new.nwk"  "EPA80_2S_new.nwk"  "EPA80_3S_new.nwk"  "EPA80_4S_new.nwk"  "EPA80_5S_new.nwk"  "EPA80_6S_new.nwk"  "EPA80_7S_new.nwk"  "EPA80_8S_new.nwk"  "EPA80_9S_new.nwk" 
# [41] "EPA99_10S_new.nwk" "EPA99_1S_new.nwk"  "EPA99_2S_new.nwk"  "EPA99_3S_new.nwk"  "EPA99_4S_new.nwk"  "EPA99_5S_new.nwk"  "EPA99_6S_new.nwk"  "EPA99_7S_new.nwk"  "EPA99_8S_new.nwk"  "EPA99_9S_new.nwk" 


# we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:

treeList$`EPA80_1S_new.nwk`

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)

species_List <- foreach(i=1:length(placementList)) %do% placementList[[i]]$seq.name
species_List

# #Read fasta files of placed co1 on the bb tree
# R_co1_20_1=read.fasta(file = "CO1_20_1.fasta")
# R_co1_20_2=read.fasta(file = "CO1_20_2.fasta")

#read the refernce tree
ReferenceTree <- read.tree("RAxML_bestTree.FR100_new")
getwd()

# R20_1_tree<-read.tree("EPA20_1R_new.nwk")
# R20_2_tree<-read.tree("EPA20_2R_new.nwk")
# R80_2_tree<-read.tree("EPA80_2R_new.nwk")
# R80_1_tree<-read.tree("EPA80_1R_new.nwk")

#------
#I need the seq.names for downstream analysis. Hence, assign for variables 

# R_co1_20_1_species <- R_co1_20_1$seq.name
# class(R_co1_20_1_species)
# R_co1_20_2_species <- R_co1_20_2$seq.name

getwd()
##bb 80% tree when placed 20% co1
library(ips)
library(dplyr)
library(tidyr)


#new_code-----------------
#To remove single species in genus groups
Ref_Tree <- as.list(ReferenceTree)
Ref_Tree_Tips <- Ref_Tree$tip.label
class(Ref_Tree_Tips)
Ref_Tree_Tips <- as.data.frame(Ref_Tree_Tips)
df_ReferenceTree<-data.frame(str_extract(Ref_Tree_Tips$Ref_Tree_Tips , "[^_]+"))
colnames(df_ReferenceTree) <- c("genus")
names(df_ReferenceTree)
#Get ununique species by removing single species
genus_more_1<- df_ReferenceTree %>% 
  group_by(genus) %>% 
  filter(n()>=2)
#Get distinct names
Ref_Tree_distinct <- distinct(genus_more_1)#766 obs.



#get the genus name list of co1 placements ()To replace ids with co1 placement genus names

names(placementList$`CO1_20_1S.fasta`)
#head(R_co1_20_1$seq.name)
#As I need to get the genus name (1st part of two parts of species name), follow the below
#check placementList order
names(placementList)
# nw <- as.data.frame(str_extract(species_List[[11]],"[^_]+"))
#for all the placement list with genus name
placement_genus_list <- foreach(i=1:length(species_List)) %do% as.data.frame(str_extract(species_List[[i]],"[^_]+"))

library(purrr) # using purrr for map function
# map setNames to each dataframe in the list
placement_genus_list_new <- map(placement_genus_list , setNames, nm = "genus")
# # then list2env should work perfectly
# list2env(placement_genus_list_new,envir = .GlobalEnv)

#to get individual dataframes for each replicate
list2env(setNames(placement_genus_list_new, paste0("newdfF", seq_along(placement_genus_list_new))),
         envir=.GlobalEnv)

#Adding ids to the dataframes. id 1 means list element 1 
id <- c(1:2469)
co1_60_list <- c(newdfF11,newdfF12,newdfF13,newdfF14,newdfF15,newdfF16,newdfF17,newdfF18,newdfF19,newdfF20)
co1_60_id_list <- foreach(i=1:length(co1_60_list)) %do% data.frame(id,co1_60_list[i])



id <- c(1:3283)
co1_80_list <- c(newdfF1,newdfF2,newdfF3,newdfF4,newdfF5,newdfF6,newdfF7,newdfF8,newdfF9,newdfF10)
co1_80_id_list <- foreach(i=1:length(co1_80_list)) %do% data.frame(id,co1_80_list[i])


id <- c(1:1654)
co1_40_list <- c(newdfF21,newdfF22,newdfF23,newdfF24,newdfF25,newdfF26,newdfF27,newdfF28,newdfF29,newdfF30)
co1_40_id_list <- foreach(i=1:length(co1_40_list)) %do% data.frame(id,co1_40_list[i])
id <- c(1:813)
co1_20_list <- c(newdfF31,newdfF32,newdfF33,newdfF34,newdfF35,newdfF36,newdfF37,newdfF38,newdfF39,newdfF40)
co1_20_id_list <- foreach(i=1:length(co1_20_list)) %do% data.frame(id,co1_20_list[i])

id <- c(1:27)
co1_1_list <- c(newdfF41,newdfF42,newdfF43,newdfF44,newdfF45,newdfF46,newdfF47,newdfF48,newdfF49,newdfF50)
co1_1_id_list <- foreach(i=1:length(co1_1_list)) %do% data.frame(id,co1_1_list[i])
# placement_genus_list<- as.data.frame(str_extract(species_List[[i]],"[^_]+"))
# R_co1_20_1_G<-as.data.frame(str_extract(R_co1_20_1$seq.name , "[^_]+"))
#R_co1_20_1_G<-as.data.frame(str_extract(R_co1_20_1 , "[^_]+"))
# class(R_co1_20_1_G)
# colnames(R_co1_20_1_G) <- "genus"
# class(R_co1_20_1_G)

# 
# #Sister species of placement sequences in reference tree
# RefTree_R_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(ReferenceTree,R_co1_20_1_species[i],type="terminal",label=T)
# #Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
# R_80_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(R80_1_tree,R_co1_20_1_species[i],type="terminal",label=T)

# [1] "CO1_20_10S.fasta" "CO1_20_1S.fasta"  "CO1_20_2S.fasta"  "CO1_20_3S.fasta"  "CO1_20_4S.fasta"  "CO1_20_5S.fasta"  "CO1_20_6S.fasta"  "CO1_20_7S.fasta"  "CO1_20_8S.fasta"  "CO1_20_9S.fasta" 
# [11] "CO1_40_10S.fasta" "CO1_40_1S.fasta"  "CO1_40_2S.fasta"  "CO1_40_3S.fasta"  "CO1_40_4S.fasta"  "CO1_40_5S.fasta"  "CO1_40_6S.fasta"  "CO1_40_7S.fasta"  "CO1_40_8S.fasta"  "CO1_40_9S.fasta" 
# [21] "CO1_60_10S.fasta" "CO1_60_1S.fasta"  "CO1_60_2S.fasta"  "CO1_60_3S.fasta"  "CO1_60_4S.fasta"  "CO1_60_5S.fasta"  "CO1_60_6S.fasta"  "CO1_60_7S.fasta"  "CO1_60_8S.fasta"  "CO1_60_9S.fasta" 
# [31] "CO1_80_10S.fasta" "CO1_80_1S.fasta"  "CO1_80_2S.fasta"  "CO1_80_3S.fasta"  "CO1_80_4S.fasta"  "CO1_80_5S.fasta"  "CO1_80_6S.fasta"  "CO1_80_7S.fasta"  "CO1_80_8S.fasta"  "CO1_80_9S.fasta" 
# [41] "CO1_99_10S.fasta" "CO1_99_1S.fasta"  "CO1_99_2S.fasta"  "CO1_99_3S.fasta"  "CO1_99_4S.fasta"  "CO1_99_5S.fasta"  "CO1_99_6S.fasta"  "CO1_99_7S.fasta"  "CO1_99_8S.fasta"  "CO1_99_9S.fasta" # we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:

#Sister species of placement sequences in reference tree
RefTree_R_co1_20_1_List <- foreach(i=1:length(species_List[[1]])) %do% sister(ReferenceTree,species_List[[1]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_80_co1_20_1_List <- foreach(i=1:length(species_List[[1]])) %do% sister(treeList$EPA80_10S_new.nwk,species_List[[1]][i],type="terminal",label=T)


#####
#EPAng stratified sister position (EPAng using raxml info file as the model)
#Last updated Jan 19, 2024
#Jan 15, 2024
library(ape)
#Sister species table
getwd()
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
R_co1_1_1=read.fasta(file = "CO1_99_1S.fasta")
R_co1_1_2=read.fasta(file = "CO1_99_2S.fasta")
R_co1_1_3=read.fasta(file = "CO1_99_3S.fasta")
R_co1_1_4=read.fasta(file = "CO1_99_4S.fasta")
R_co1_1_5=read.fasta(file = "CO1_99_5S.fasta")
R_co1_1_6=read.fasta(file = "CO1_99_6S.fasta")
R_co1_1_7=read.fasta(file = "CO1_99_7S.fasta")
R_co1_1_8=read.fasta(file = "CO1_99_8S.fasta")
R_co1_1_9=read.fasta(file = "CO1_99_9S.fasta")
R_co1_1_10=read.fasta(file = "CO1_99_10S.fasta")

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
# R60_1_tree<-read.tree("epa80_1S.nwk")
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
##bb 80% tree when placed 20% co1
library(ips)
#sister(tree, "t4", type = "terminal", label = T)
RefTree_R_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(ReferenceTree,R_co1_20_1_species[i],type="terminal",label=T)
R_80_co1_20_1_List <- foreach(i=1:length(R_co1_20_1_species)) %do% sister(R80_1_tree,R_co1_20_1_species[i],type="terminal",label=T) #Error in sister(R80_1_tree, R_co1_20_1_species[i], type = "terminal",  : 
#task 7 failed - "argument is of length zero"

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_1_List [(RefTree_R_co1_20_1_List %in% R_80_co1_20_1_List )] #489
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_1_List [!(RefTree_R_co1_20_1_List %in% R_80_co1_20_1_List)] #339

RefTree_R_co1_20_2_List <- foreach(i=1:length(R_co1_20_2_species)) %do% sister(ReferenceTree,R_co1_20_2_species[i],type="terminal",label=T)
R_80_co1_20_2_List <- foreach(i=1:length(R_co1_20_2_species)) %do% sister(R80_2_tree,R_co1_20_2_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_2_List [(RefTree_R_co1_20_2_List %in% R_80_co1_20_2_List )] #474
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_2_List [!(RefTree_R_co1_20_2_List %in% R_80_co1_20_2_List)] #354

RefTree_R_co1_20_3_List <- foreach(i=1:length(R_co1_20_3_species)) %do% sister(ReferenceTree,R_co1_20_3_species[i],type="terminal",label=T)
R_80_co1_20_3_List <- foreach(i=1:length(R_co1_20_3_species)) %do% sister(R80_3_tree,R_co1_20_3_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_3_List [(RefTree_R_co1_20_3_List %in% R_80_co1_20_3_List )] #477
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_3_List [!(RefTree_R_co1_20_3_List %in% R_80_co1_20_3_List)] #351

RefTree_R_co1_20_4_List <- foreach(i=1:length(R_co1_20_4_species)) %do% sister(ReferenceTree,R_co1_20_4_species[i],type="terminal",label=T)
R_80_co1_20_4_List <- foreach(i=1:length(R_co1_20_4_species)) %do% sister(R80_4_tree,R_co1_20_4_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_4_List [(RefTree_R_co1_20_4_List %in% R_80_co1_20_4_List )] #488
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_4_List [!(RefTree_R_co1_20_4_List %in% R_80_co1_20_4_List)] #340

RefTree_R_co1_20_5_List <- foreach(i=1:length(R_co1_20_5_species)) %do% sister(ReferenceTree,R_co1_20_5_species[i],type="terminal",label=T)
R_80_co1_20_5_List <- foreach(i=1:length(R_co1_20_5_species)) %do% sister(R80_5_tree,R_co1_20_5_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_5_List [(RefTree_R_co1_20_5_List %in% R_80_co1_20_5_List )] #458
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_5_List [!(RefTree_R_co1_20_5_List %in% R_80_co1_20_5_List)] #370

RefTree_R_co1_20_6_List <- foreach(i=1:length(R_co1_20_6_species)) %do% sister(ReferenceTree,R_co1_20_6_species[i],type="terminal",label=T)
R_80_co1_20_6_List <- foreach(i=1:length(R_co1_20_6_species)) %do% sister(R80_6_tree,R_co1_20_6_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_6_List [(RefTree_R_co1_20_6_List %in% R_80_co1_20_6_List )] #426
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_6_List [!(RefTree_R_co1_20_6_List %in% R_80_co1_20_6_List)] #402

RefTree_R_co1_20_7_List <- foreach(i=1:length(R_co1_20_7_species)) %do% sister(ReferenceTree,R_co1_20_7_species[i],type="terminal",label=T)
R_80_co1_20_7_List <- foreach(i=1:length(R_co1_20_7_species)) %do% sister(R80_7_tree,R_co1_20_7_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_7_List [(RefTree_R_co1_20_7_List %in% R_80_co1_20_7_List )] #467
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_7_List [!(RefTree_R_co1_20_7_List %in% R_80_co1_20_7_List)] #361

RefTree_R_co1_20_8_List <- foreach(i=1:length(R_co1_20_8_species)) %do% sister(ReferenceTree,R_co1_20_8_species[i],type="terminal",label=T)
R_80_co1_20_8_List <- foreach(i=1:length(R_co1_20_8_species)) %do% sister(R80_8_tree,R_co1_20_8_species[i],type="terminal",label=T)



#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_8_List [(RefTree_R_co1_20_8_List %in% R_80_co1_20_8_List )] #485
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_8_List [!(RefTree_R_co1_20_8_List %in% R_80_co1_20_8_List)] #343

RefTree_R_co1_20_9_List <- foreach(i=1:length(R_co1_20_9_species)) %do% sister(ReferenceTree,R_co1_20_9_species[i],type="terminal",label=T)
R_80_co1_20_9_List <- foreach(i=1:length(R_co1_20_9_species)) %do% sister(R80_9_tree,R_co1_20_9_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_9_List [(RefTree_R_co1_20_9_List %in% R_80_co1_20_9_List )] #461
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_9_List [!(RefTree_R_co1_20_9_List %in% R_80_co1_20_9_List)] #367

RefTree_R_co1_20_10_List <- foreach(i=1:length(R_co1_20_10_species)) %do% sister(ReferenceTree,R_co1_20_10_species[i],type="terminal",label=T)
R_80_co1_20_10_List <- foreach(i=1:length(R_co1_20_10_species)) %do% sister(R80_10_tree,R_co1_20_10_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_20_10_List [(RefTree_R_co1_20_10_List %in% R_80_co1_20_10_List )] #468
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_20_10_List [!(RefTree_R_co1_20_10_List %in% R_80_co1_20_10_List)] #360

#bb 60% tree when placed 40% co1

RefTree_R_co1_40_1_List <- foreach(i=1:length(R_co1_40_1_species)) %do% sister(ReferenceTree,R_co1_40_1_species[i],type="terminal",label=T)
R_60_co1_40_1_List <- foreach(i=1:length(R_co1_40_1_species)) %do% sister(R40_1_tree,R_co1_40_1_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_1_List [(RefTree_R_co1_40_1_List %in% R_60_co1_40_1_List )] #640
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_1_List [!(RefTree_R_co1_40_1_List %in% R_60_co1_40_1_List )] #1829

RefTree_R_co1_40_2_List <- foreach(i=1:length(R_co1_40_2_species)) %do% sister(ReferenceTree,R_co1_40_2_species[i],type="terminal",label=T)
R_60_co1_40_2_List <- foreach(i=1:length(R_co1_40_2_species)) %do% sister(R40_2_tree,R_co1_40_2_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_2_List [(RefTree_R_co1_40_2_List %in% R_60_co1_40_2_List )] #607
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_2_List [!(RefTree_R_co1_40_2_List %in% R_60_co1_40_2_List )] #1862

RefTree_R_co1_40_3_List <- foreach(i=1:length(R_co1_40_3_species)) %do% sister(ReferenceTree,R_co1_40_3_species[i],type="terminal",label=T)
R_60_co1_40_3_List <- foreach(i=1:length(R_co1_40_3_species)) %do% sister(R40_3_tree,R_co1_40_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_3_List [(RefTree_R_co1_40_3_List %in% R_60_co1_40_3_List )] #644
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_3_List [!(RefTree_R_co1_40_3_List %in% R_60_co1_40_3_List )] #1825

RefTree_R_co1_40_4_List <- foreach(i=1:length(R_co1_40_4_species)) %do% sister(ReferenceTree,R_co1_40_4_species[i],type="terminal",label=T)
R_60_co1_40_4_List <- foreach(i=1:length(R_co1_40_4_species)) %do% sister(R40_4_tree,R_co1_40_4_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_4_List [(RefTree_R_co1_40_4_List %in% R_60_co1_40_4_List )] #518
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_4_List [!(RefTree_R_co1_40_4_List %in% R_60_co1_40_4_List )] #1951

# RefTree_R_co1_40_5_List <- foreach(i=1:length(R_co1_40_5_species)) %do% sister(ReferenceTree,R_co1_40_5_species[i],type="terminal",label=T)
# R_60_co1_40_5_List <- foreach(i=1:length(R_co1_40_5_species)) %do% sister(R40_5_tree,R_co1_40_5_species[i],type="terminal",label=T)
# # 
# # 
# # #get the number of species have similar sisters by comparing with reference tree
# RefTree_R_co1_40_5_List [(RefTree_R_co1_40_5_List %in% R_60_co1_40_5_List )] #578
# # #get the number of species doesn't have similar sisters by comparing with reference tree
# RefTree_R_co1_40_5_List [!(RefTree_R_co1_40_5_List %in% R_60_co1_40_5_List )] #1230

RefTree_R_co1_40_6_List <- foreach(i=1:length(R_co1_40_6_species)) %do% sister(ReferenceTree,R_co1_40_6_species[i],type="terminal",label=T)
R_60_co1_40_6_List <- foreach(i=1:length(R_co1_40_6_species)) %do% sister(R40_6_tree,R_co1_40_6_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_6_List [(RefTree_R_co1_40_6_List %in% R_60_co1_40_6_List )] #598
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_6_List [!(RefTree_R_co1_40_6_List %in% R_60_co1_40_6_List )] #1871

RefTree_R_co1_40_7_List <- foreach(i=1:length(R_co1_40_7_species)) %do% sister(ReferenceTree,R_co1_40_7_species[i],type="terminal",label=T)
R_60_co1_40_7_List <- foreach(i=1:length(R_co1_40_7_species)) %do% sister(R40_7_tree,R_co1_40_7_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_7_List [(RefTree_R_co1_40_7_List %in% R_60_co1_40_7_List )] #591
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_7_List [!(RefTree_R_co1_40_7_List %in% R_60_co1_40_7_List )] #1878

RefTree_R_co1_40_8_List <- foreach(i=1:length(R_co1_40_8_species)) %do% sister(ReferenceTree,R_co1_40_8_species[i],type="terminal",label=T)
R_60_co1_40_8_List <- foreach(i=1:length(R_co1_40_8_species)) %do% sister(R40_8_tree,R_co1_40_8_species[i],type="terminal",label=T)
#sister(R60_8_tree,"Synodontis_serratus",type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_8_List [(RefTree_R_co1_40_8_List %in% R_60_co1_40_8_List )] #626
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_8_List [!(RefTree_R_co1_40_8_List %in% R_60_co1_40_8_List )] #1843

RefTree_R_co1_40_9_List <- foreach(i=1:length(R_co1_40_9_species)) %do% sister(ReferenceTree,R_co1_40_9_species[i],type="terminal",label=T)
R_60_co1_40_9_List <- foreach(i=1:length(R_co1_40_9_species)) %do% sister(R40_9_tree,R_co1_40_9_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_9_List [(RefTree_R_co1_40_9_List %in% R_60_co1_40_9_List )] #629
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_9_List [!(RefTree_R_co1_40_9_List %in% R_60_co1_40_9_List )] #1840

RefTree_R_co1_40_10_List <- foreach(i=1:length(R_co1_40_10_species)) %do% sister(ReferenceTree,R_co1_40_10_species[i],type="terminal",label=T)
R_60_co1_40_10_List <- foreach(i=1:length(R_co1_40_10_species)) %do% sister(R40_10_tree,R_co1_40_10_species[i],type="terminal",label=T)


#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_40_10_List [(RefTree_R_co1_40_10_List %in% R_60_co1_40_10_List )] #630
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_40_10_List [!(RefTree_R_co1_40_10_List %in% R_60_co1_40_10_List )] #1839

#RefTree_biasco1List[sapply(names(RefTree_biasco1List), function(x) !identical(RefTree_biasco1List[[x]], biasco1List[[x]]))] 

#bb 60% tree when placed 40% co1

RefTree_R_co1_60_1_List <- foreach(i=1:length(R_co1_60_1_species)) %do% sister(ReferenceTree,R_co1_60_1_species[i],type="terminal",label=T)
R_40_co1_60_1_List <- foreach(i=1:length(R_co1_60_1_species)) %do% sister(R60_1_tree,R_co1_60_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_1_List [(RefTree_R_co1_60_1_List %in% R_40_co1_60_1_List )] #536
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_1_List [!(RefTree_R_co1_60_1_List %in% R_40_co1_60_1_List )] #1118

RefTree_R_co1_60_2_List <- foreach(i=1:length(R_co1_60_2_species)) %do% sister(ReferenceTree,R_co1_60_2_species[i],type="terminal",label=T)
R_40_co1_60_2_List <- foreach(i=1:length(R_co1_60_2_species)) %do% sister(R60_2_tree,R_co1_60_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_2_List [(RefTree_R_co1_60_2_List %in% R_40_co1_60_2_List )] #566
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_2_List [!(RefTree_R_co1_60_2_List %in% R_40_co1_60_2_List )] #1088


RefTree_R_co1_60_3_List <- foreach(i=1:length(R_co1_60_3_species)) %do% sister(ReferenceTree,R_co1_60_3_species[i],type="terminal",label=T)
R_40_co1_60_3_List <- foreach(i=1:length(R_co1_60_3_species)) %do% sister(R60_3_tree,R_co1_60_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_3_List [(RefTree_R_co1_60_3_List %in% R_40_co1_60_3_List )] #753
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_3_List [!(RefTree_R_co1_60_3_List %in% R_40_co1_60_3_List )] #901

RefTree_R_co1_60_4_List <- foreach(i=1:length(R_co1_60_4_species)) %do% sister(ReferenceTree,R_co1_60_4_species[i],type="terminal",label=T)
R_40_co1_60_4_List <- foreach(i=1:length(R_co1_60_4_species)) %do% sister(R60_4_tree,R_co1_60_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_4_List [(RefTree_R_co1_60_4_List %in% R_40_co1_60_4_List )] #744
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_4_List [!(RefTree_R_co1_60_4_List %in% R_40_co1_60_4_List )] #910

RefTree_R_co1_60_5_List <- foreach(i=1:length(R_co1_60_5_species)) %do% sister(ReferenceTree,R_co1_60_5_species[i],type="terminal",label=T)
R_40_co1_60_5_List <- foreach(i=1:length(R_co1_60_5_species)) %do% sister(R60_5_tree,R_co1_60_5_species[i],type="terminal",label=T)
# 
# #get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_5_List [(RefTree_R_co1_60_5_List %in% R_40_co1_60_5_List )] #763
# #get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_5_List [!(RefTree_R_co1_60_5_List %in% R_40_co1_60_5_List )] #891

RefTree_R_co1_60_6_List <- foreach(i=1:length(R_co1_60_6_species)) %do% sister(ReferenceTree,R_co1_60_6_species[i],type="terminal",label=T)
R_40_co1_60_6_List <- foreach(i=1:length(R_co1_60_6_species)) %do% sister(R60_6_tree,R_co1_60_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_6_List [(RefTree_R_co1_60_6_List %in% R_40_co1_60_6_List )] #750
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_6_List [!(RefTree_R_co1_60_6_List %in% R_40_co1_60_6_List )] #904

RefTree_R_co1_60_7_List <- foreach(i=1:length(R_co1_60_7_species)) %do% sister(ReferenceTree,R_co1_60_7_species[i],type="terminal",label=T)
R_40_co1_60_7_List <- foreach(i=1:length(R_co1_60_7_species)) %do% sister(R60_7_tree,R_co1_60_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_7_List [(RefTree_R_co1_60_7_List %in% R_40_co1_60_7_List )] #505
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_7_List [!(RefTree_R_co1_60_7_List %in% R_40_co1_60_7_List )] #1149

RefTree_R_co1_60_8_List <- foreach(i=1:length(R_co1_60_8_species)) %do% sister(ReferenceTree,R_co1_60_8_species[i],type="terminal",label=T)
R_40_co1_60_8_List <- foreach(i=1:length(R_co1_60_8_species)) %do% sister(R60_8_tree,R_co1_60_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_8_List [(RefTree_R_co1_60_8_List %in% R_40_co1_60_8_List )] #561
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_8_List [!(RefTree_R_co1_60_8_List %in% R_40_co1_60_8_List )] #1093

RefTree_R_co1_60_9_List <- foreach(i=1:length(R_co1_60_9_species)) %do% sister(ReferenceTree,R_co1_60_9_species[i],type="terminal",label=T)
R_40_co1_60_9_List <- foreach(i=1:length(R_co1_60_9_species)) %do% sister(R60_9_tree,R_co1_60_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_9_List [(RefTree_R_co1_60_9_List %in% R_40_co1_60_9_List )] #513
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_9_List [!(RefTree_R_co1_60_9_List %in% R_40_co1_60_9_List )] #1141

RefTree_R_co1_60_10_List <- foreach(i=1:length(R_co1_60_10_species)) %do% sister(ReferenceTree,R_co1_60_10_species[i],type="terminal",label=T)
R_40_co1_60_10_List <- foreach(i=1:length(R_co1_60_10_species)) %do% sister(R60_10_tree,R_co1_60_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_60_10_List [(RefTree_R_co1_60_10_List %in% R_40_co1_60_10_List )] #748
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_60_10_List [!(RefTree_R_co1_60_10_List %in% R_40_co1_60_10_List )] #906


#bb 80% tree when placed 20% co1

RefTree_R_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(ReferenceTree,R_co1_80_1_species[i],type="terminal",label=T)
R_20_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(R80_1_tree,R_co1_80_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_1_List [(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #342
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_1_List [!(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #471

RefTree_R_co1_80_2_List <- foreach(i=1:length(R_co1_80_2_species)) %do% sister(ReferenceTree,R_co1_80_2_species[i],type="terminal",label=T)
R_20_co1_80_2_List <- foreach(i=1:length(R_co1_80_2_species)) %do% sister(R80_2_tree,R_co1_80_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_2_List [(RefTree_R_co1_80_2_List %in% R_20_co1_80_2_List )] #301
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_2_List [!(RefTree_R_co1_80_2_List %in% R_20_co1_80_2_List )] #512

RefTree_R_co1_80_3_List <- foreach(i=1:length(R_co1_80_3_species)) %do% sister(ReferenceTree,R_co1_80_3_species[i],type="terminal",label=T)
R_20_co1_80_3_List <- foreach(i=1:length(R_co1_80_3_species)) %do% sister(R80_3_tree,R_co1_80_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_3_List [(RefTree_R_co1_80_3_List %in% R_20_co1_80_3_List )] #425
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_3_List [!(RefTree_R_co1_80_3_List %in% R_20_co1_80_3_List )] #388

RefTree_R_co1_80_4_List <- foreach(i=1:length(R_co1_80_4_species)) %do% sister(ReferenceTree,R_co1_80_4_species[i],type="terminal",label=T)
R_20_co1_80_4_List <- foreach(i=1:length(R_co1_80_4_species)) %do% sister(R80_4_tree,R_co1_80_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_4_List [(RefTree_R_co1_80_4_List %in% R_20_co1_80_4_List )] #334
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_4_List [!(RefTree_R_co1_80_4_List %in% R_20_co1_80_4_List )] #479

RefTree_R_co1_80_5_List <- foreach(i=1:length(R_co1_80_5_species)) %do% sister(ReferenceTree,R_co1_80_5_species[i],type="terminal",label=T)
R_20_co1_80_5_List <- foreach(i=1:length(R_co1_80_5_species)) %do% sister(R80_5_tree,R_co1_80_5_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_5_List [(RefTree_R_co1_80_5_List %in% R_20_co1_80_5_List )] #460
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_5_List [!(RefTree_R_co1_80_5_List %in% R_20_co1_80_5_List )] #353

RefTree_R_co1_80_6_List <- foreach(i=1:length(R_co1_80_6_species)) %do% sister(ReferenceTree,R_co1_80_6_species[i],type="terminal",label=T)
R_20_co1_80_6_List <- foreach(i=1:length(R_co1_80_6_species)) %do% sister(R80_6_tree,R_co1_80_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_6_List [(RefTree_R_co1_80_6_List %in% R_20_co1_80_6_List )] #458
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_6_List [!(RefTree_R_co1_80_6_List %in% R_20_co1_80_6_List )] #355

RefTree_R_co1_80_7_List <- foreach(i=1:length(R_co1_80_7_species)) %do% sister(ReferenceTree,R_co1_80_7_species[i],type="terminal",label=T)
R_20_co1_80_7_List <- foreach(i=1:length(R_co1_80_7_species)) %do% sister(R80_7_tree,R_co1_80_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_7_List [(RefTree_R_co1_80_7_List %in% R_20_co1_80_7_List )] #343
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_7_List [!(RefTree_R_co1_80_7_List %in% R_20_co1_80_7_List )] #470

RefTree_R_co1_80_8_List <- foreach(i=1:length(R_co1_80_8_species)) %do% sister(ReferenceTree,R_co1_80_8_species[i],type="terminal",label=T)
R_20_co1_80_8_List <- foreach(i=1:length(R_co1_80_8_species)) %do% sister(R80_8_tree,R_co1_80_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_8_List [(RefTree_R_co1_80_8_List %in% R_20_co1_80_8_List )] #326
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_8_List [!(RefTree_R_co1_80_8_List %in% R_20_co1_80_8_List )] #487

RefTree_R_co1_80_9_List <- foreach(i=1:length(R_co1_80_9_species)) %do% sister(ReferenceTree,R_co1_80_9_species[i],type="terminal",label=T)
R_20_co1_80_9_List <- foreach(i=1:length(R_co1_80_9_species)) %do% sister(R80_9_tree,R_co1_80_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_9_List [(RefTree_R_co1_80_9_List %in% R_20_co1_80_9_List )] #460
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_9_List [!(RefTree_R_co1_80_9_List %in% R_20_co1_80_9_List )] #353

RefTree_R_co1_80_10_List <- foreach(i=1:length(R_co1_80_10_species)) %do% sister(ReferenceTree,R_co1_80_10_species[i],type="terminal",label=T)
R_20_co1_80_10_List <- foreach(i=1:length(R_co1_80_10_species)) %do% sister(R80_10_tree,R_co1_80_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_80_10_List [(RefTree_R_co1_80_10_List %in% R_20_co1_80_10_List )] #463
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_80_10_List [!(RefTree_R_co1_80_10_List %in% R_20_co1_80_10_List )] #350
#---------------------


#bb 99% tree when placed 1% co1
RefTree_R_co1_1_1_List <- foreach(i=1:length(R_co1_1_1_species)) %do% sister(ReferenceTree,R_co1_1_1_species[i],type="terminal",label=T)
R_99_co1_1_1_List <- foreach(i=1:length(R_co1_1_1_species)) %do% sister(R99_1_tree,R_co1_1_1_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_1_List [(RefTree_R_co1_1_1_List %in% R_99_co1_1_1_List )] #15
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_1_List [!(RefTree_R_co1_1_1_List %in% R_99_co1_1_1_List)] #12

RefTree_R_co1_1_2_List <- foreach(i=1:length(R_co1_1_2_species)) %do% sister(ReferenceTree,R_co1_1_2_species[i],type="terminal",label=T)
R_99_co1_1_2_List <- foreach(i=1:length(R_co1_1_2_species)) %do% sister(R99_2_tree,R_co1_1_2_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_2_List [(RefTree_R_co1_1_2_List %in% R_99_co1_1_2_List )] #10
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_2_List [!(RefTree_R_co1_1_2_List %in% R_99_co1_1_2_List)] #17

RefTree_R_co1_1_3_List <- foreach(i=1:length(R_co1_1_3_species)) %do% sister(ReferenceTree,R_co1_1_3_species[i],type="terminal",label=T)
R_99_co1_1_3_List <- foreach(i=1:length(R_co1_1_3_species)) %do% sister(R99_3_tree,R_co1_1_3_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_3_List [(RefTree_R_co1_1_3_List %in% R_99_co1_1_3_List )] #12
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_3_List [!(RefTree_R_co1_1_3_List %in% R_99_co1_1_3_List)] #15

RefTree_R_co1_1_4_List <- foreach(i=1:length(R_co1_1_4_species)) %do% sister(ReferenceTree,R_co1_1_4_species[i],type="terminal",label=T)
R_99_co1_1_4_List <- foreach(i=1:length(R_co1_1_4_species)) %do% sister(R99_4_tree,R_co1_1_4_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_4_List [(RefTree_R_co1_1_4_List %in% R_99_co1_1_4_List )] #14
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_4_List [!(RefTree_R_co1_1_4_List %in% R_99_co1_1_4_List)] #13

RefTree_R_co1_1_5_List <- foreach(i=1:length(R_co1_1_5_species)) %do% sister(ReferenceTree,R_co1_1_5_species[i],type="terminal",label=T)
R_99_co1_1_5_List <- foreach(i=1:length(R_co1_1_5_species)) %do% sister(R99_5_tree,R_co1_1_5_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_5_List [(RefTree_R_co1_1_5_List %in% R_99_co1_1_5_List )] #13
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_5_List [!(RefTree_R_co1_1_5_List %in% R_99_co1_1_5_List)] #14

RefTree_R_co1_1_6_List <- foreach(i=1:length(R_co1_1_6_species)) %do% sister(ReferenceTree,R_co1_1_6_species[i],type="terminal",label=T)
R_99_co1_1_6_List <- foreach(i=1:length(R_co1_1_6_species)) %do% sister(R99_6_tree,R_co1_1_6_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_6_List [(RefTree_R_co1_1_6_List %in% R_99_co1_1_6_List )] #9
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_6_List [!(RefTree_R_co1_1_6_List %in% R_99_co1_1_6_List)] #21

RefTree_R_co1_1_7_List <- foreach(i=1:length(R_co1_1_7_species)) %do% sister(ReferenceTree,R_co1_1_7_species[i],type="terminal",label=T)
R_99_co1_1_7_List <- foreach(i=1:length(R_co1_1_7_species)) %do% sister(R99_7_tree,R_co1_1_7_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_7_List [(RefTree_R_co1_1_7_List %in% R_99_co1_1_7_List )] #10
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_7_List [!(RefTree_R_co1_1_7_List %in% R_99_co1_1_7_List)] #17

RefTree_R_co1_1_8_List <- foreach(i=1:length(R_co1_1_8_species)) %do% sister(ReferenceTree,R_co1_1_8_species[i],type="terminal",label=T)
R_99_co1_1_8_List <- foreach(i=1:length(R_co1_1_8_species)) %do% sister(R99_8_tree,R_co1_1_8_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_8_List [(RefTree_R_co1_1_8_List %in% R_99_co1_1_8_List )] #11
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_8_List [!(RefTree_R_co1_1_8_List %in% R_99_co1_1_8_List)] #16

RefTree_R_co1_1_9_List <- foreach(i=1:length(R_co1_1_9_species)) %do% sister(ReferenceTree,R_co1_1_9_species[i],type="terminal",label=T)
R_99_co1_1_9_List <- foreach(i=1:length(R_co1_1_9_species)) %do% sister(R99_9_tree,R_co1_1_9_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_9_List [(RefTree_R_co1_1_9_List %in% R_99_co1_1_9_List )] #12
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_9_List [!(RefTree_R_co1_1_9_List %in% R_99_co1_1_9_List)] #15


RefTree_R_co1_1_10_List <- foreach(i=1:length(R_co1_1_10_species)) %do% sister(ReferenceTree,R_co1_1_10_species[i],type="terminal",label=T)
R_99_co1_1_10_List <- foreach(i=1:length(R_co1_1_10_species)) %do% sister(R99_10_tree,R_co1_1_10_species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_R_co1_1_10_List [(RefTree_R_co1_1_10_List %in% R_99_co1_1_10_List )] #12
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_R_co1_1_10_List [!(RefTree_R_co1_1_10_List %in% R_99_co1_1_10_List)] #15







#####

#Load the stringr package, data.tabl, tidyverse & dplyr packages
library(stringr)
library(data.table)
library(tidyverse)
library(dplyr)

# # extract_Genus function first
# df_RefTree_R_co1_20_1 <- extract_Genus(RefTree_R_co1_20_1_List)
# df_R_80_co1_20_1 <- extract_Genus(R_80_co1_20_1_List)

# generics::intersect(df_RefTree_R_co1_20_1,df_R_80_co1_20_1 )
# generics::setdiff(df_RefTree_R_co1_20_1,df_R_80_co1_20_1 )
#Function for bb80% (co1 20%) #(For the bb60% have to change only one word in co1_20_id_list to co1_40_id_list)
genus_count <- function(RefTree,placementTree,N){
  df_RefTree <- extract_Genus(RefTree)
  df_placement <- extract_Genus(placementTree)
  #get the rows of placement genus that's not matching with reference tree
  # mismatches <- setdiff(df_RefTree,df_placement)
  nn <- union(df_RefTree ,df_placement)
  inter <- intersect(df_RefTree ,df_placement)
  mismatches <- anti_join(nn,inter)
  #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
  uniq <- count(unique(mismatches,by="id"))
  #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
  #First I get the genus name to particular the ids 
  mismatch_wt_names <- inner_join(mismatches,co1_20_id_list[[N]],by="id")
  #Get the number of placement genus without matches to the placement genus in reference tree
  unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
  unique_mismatches_count <- count(unique_mismatches)
  #genus.y col in unique_mismatches dataframe has the genus list of mismatches
  #Using that col and set diff function get the monotypic genus
  monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
  monotypic_genus_count <- length(monotypic_genus)
  #Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
  match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
  
  mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
  #Let's see whats common in both (mixture set)
  #intersect of id
  mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10
  mixtures_count <- length(mixtures)
  #Only in match_monotypic
  All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
  All_match_monotypic_count <- length(All_match_monotypic)
  #Only in mismatch_monotypic
  
  All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144
  All_mismatch_monotypic_count <- length(All_mismatch_monotypic)
  results <- data.frame(genus_mismatches=uniq,unique_mismatches_count,monotypic_genus=monotypic_genus_count,mixed_monotypic=mixtures_count,All_match_monotypic=All_match_monotypic_count,All_mismatch_monotypic=All_mismatch_monotypic_count)
  return(results)
}

#Sister species of placement sequences in reference tree
RefTree_R_co1_20_1_List <- foreach(i=1:length(species_List[[39]])) %do% sister(ReferenceTree,species_List[[39]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_80_co1_20_1_List <- foreach(i=1:length(species_List[[39]])) %do% sister(treeList$EPA80_8S_new.nwk,species_List[[39]][i],type="terminal",label=T)
#For the N (order is,1= 80_10,2=80_1,3=80_2,4=80_3,5=80_4,6=80_5,7=80_6,8=80_7,9=80_8,10=80_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_20_1_List,R_80_co1_20_1_List,9)

#--------------
#For bb60% and co1 40%
genus_count <- function(RefTree,placementTree,N){
  df_RefTree <- extract_Genus(RefTree)
  df_placement <- extract_Genus(placementTree)
  #get the rows of placement genus that's not matching with reference tree
  # mismatches <- setdiff(df_RefTree,df_placement)
  nn <- union(df_RefTree ,df_placement)
  inter <- intersect(df_RefTree ,df_placement)
  mismatches <- anti_join(nn,inter)
  #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
  uniq <- count(unique(mismatches,by="id"))
  #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
  #First I get the genus name to particular the ids 
  mismatch_wt_names <- inner_join(mismatches,co1_40_id_list[[N]],by="id")
  #Get the number of placement genus without matches to the placement genus in reference tree
  unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
  unique_mismatches_count <- count(unique_mismatches)
  #genus.y col in unique_mismatches dataframe has the genus list of mismatches
  #Using that col and set diff function get the monotypic genus
  monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
  monotypic_genus_count <- length(monotypic_genus)
  #Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
  match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
  
  mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
  #Let's see whats common in both (mixture set)
  #intersect of id
  mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10
  mixtures_count <- length(mixtures)
  #Only in match_monotypic
  All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
  All_match_monotypic_count <- length(All_match_monotypic)
  #Only in mismatch_monotypic
  
  All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144
  All_mismatch_monotypic_count <- length(All_mismatch_monotypic)
  results <- data.frame(genus_mismatches=uniq,unique_mismatches_count,monotypic_genus=monotypic_genus_count,mixed_monotypic=mixtures_count,All_match_monotypic=All_match_monotypic_count,All_mismatch_monotypic=All_mismatch_monotypic_count)
  return(results)
}
#Sister species of placement sequences in reference tree#for bb60%
RefTree_R_co1_40_1_List <- foreach(i=1:length(species_List[[29]])) %do% sister(ReferenceTree,species_List[[29]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_60_co1_40_1_List <- foreach(i=1:length(species_List[[29]])) %do% sister(treeList$EPA60_8S_new.nwk,species_List[[29]][i],type="terminal",label=T)

#Call the function for bb60% and co1 40%
#For the N (order is,1= 60_1,2=60_10,3=60_2,4=60_3,5=60_4,6=60_5,7=60_6,8=60_7,9=60_8,10=60_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_40_1_List,R_60_co1_40_1_List,9)

#for bb40% and co1 60%
genus_count <- function(RefTree,placementTree,N){
  df_RefTree <- extract_Genus(RefTree)
  df_placement <- extract_Genus(placementTree)
  #get the rows of placement genus that's not matching with reference tree
  # mismatches <- setdiff(df_RefTree,df_placement)
  nn <- union(df_RefTree ,df_placement)
  inter <- intersect(df_RefTree ,df_placement)
  mismatches <- anti_join(nn,inter)
  #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
  uniq <- count(unique(mismatches,by="id"))
  #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
  #First I get the genus name to particular the ids 
  mismatch_wt_names <- inner_join(mismatches,co1_60_id_list[[N]],by="id")
  #Get the number of placement genus without matches to the placement genus in reference tree
  unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
  unique_mismatches_count <- count(unique_mismatches)
  #genus.y col in unique_mismatches dataframe has the genus list of mismatches
  #Using that col and set diff function get the monotypic genus
  monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
  monotypic_genus_count <- length(monotypic_genus)
  #Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
  match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
  
  mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
  #Let's see whats common in both (mixture set)
  #intersect of id
  mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10
  mixtures_count <- length(mixtures)
  #Only in match_monotypic
  All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
  All_match_monotypic_count <- length(All_match_monotypic)
  #Only in mismatch_monotypic
  
  All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144
  All_mismatch_monotypic_count <- length(All_mismatch_monotypic)
  results <- data.frame(genus_mismatches=uniq,unique_mismatches_count,monotypic_genus=monotypic_genus_count,mixed_monotypic=mixtures_count,All_match_monotypic=All_match_monotypic_count,All_mismatch_monotypic=All_mismatch_monotypic_count)
  return(results)
}
#Sister species of placement sequences in reference tree#for bb40%
RefTree_R_co1_60_1_List <- foreach(i=1:length(species_List[[19]])) %do% sister(ReferenceTree,species_List[[19]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_40_co1_60_1_List <- foreach(i=1:length(species_List[[19]])) %do% sister(treeList$EPA40_8S_new.nwk,species_List[[19]][i],type="terminal",label=T)

#Call the function for bb40% and co1 60%
#For the N (order is,1= 40_10,2=40_1,3=40_2,4=40_3,5=40_4,6=40_5,7=40_6,8=40_7,9=40_8,10=40_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_60_1_List,R_40_co1_60_1_List,9)

#for bb20% and co1 80%
genus_count <- function(RefTree,placementTree,N){
  df_RefTree <- extract_Genus(RefTree)
  df_placement <- extract_Genus(placementTree)
  #get the rows of placement genus that's not matching with reference tree
  #mismatches <- setdiff(df_RefTree,df_placement)
  nn <- union(df_RefTree ,df_placement)
  inter <- intersect(df_RefTree ,df_placement)
  mismatches <- anti_join(nn,inter)
  #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
  uniq <- count(unique(mismatches,by="id"))
  #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
  #First I get the genus name to particular the ids 
  mismatch_wt_names <- inner_join(mismatches,co1_80_id_list[[N]],by="id")
  #Get the number of placement genus without matches to the placement genus in reference tree
  unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
  unique_mismatches_count <- count(unique_mismatches)
  #genus.y col in unique_mismatches dataframe has the genus list of mismatches
  #Using that col and set diff function get the monotypic genus
  monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
  monotypic_genus_count <- length(monotypic_genus)
  #Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
  match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
  
  mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
  #Let's see whats common in both (mixture set)
  #intersect of id
  mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10
  mixtures_count <- length(mixtures)
  #Only in match_monotypic
  All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
  All_match_monotypic_count <- length(All_match_monotypic)
  #Only in mismatch_monotypic
  
  All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144
  All_mismatch_monotypic_count <- length(All_mismatch_monotypic)
  results <- data.frame(genus_mismatches=uniq,unique_mismatches_count,monotypic_genus=monotypic_genus_count,mixed_monotypic=mixtures_count,All_match_monotypic=All_match_monotypic_count,All_mismatch_monotypic=All_mismatch_monotypic_count)
  return(results)
}
#Sister species of placement sequences in reference tree#for bb20%
RefTree_R_co1_80_1_List <- foreach(i=1:length(species_List[[9]])) %do% sister(ReferenceTree,species_List[[9]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_20_co1_80_1_List <- foreach(i=1:length(species_List[[9]])) %do% sister(treeList$EPA20_8S_new.nwk,species_List[[9]][i],type="terminal",label=T)

#Call the function for bb40% and co1 60%
#For the N (order is,1= 20_1,2=20_10,3=20_2,4=20_3,5=20_4,6=20_5,7=20_6,8=20_7,9=20_8,10=20_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_80_1_List,R_20_co1_80_1_List,9)

##for bb99% and co1 1%
genus_count <- function(RefTree,placementTree,N){
  df_RefTree <- extract_Genus(RefTree)
  df_placement <- extract_Genus(placementTree)
  #get the rows of placement genus that's not matching with reference tree
  nn <- union(df_RefTree ,df_placement)
  inter <- intersect(df_RefTree ,df_placement)
  mismatches <- anti_join(nn,inter)
  #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
  uniq <- count(unique(mismatches,by="id"))
  #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
  #First I get the genus name to particular the ids 
  mismatch_wt_names <- inner_join(mismatches,co1_1_id_list[[N]],by="id")
  #Get the number of placement genus without matches to the placement genus in reference tree
  unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
  unique_mismatches_count <- count(unique_mismatches)
  #genus.y col in unique_mismatches dataframe has the genus list of mismatches
  #Using that col and set diff function get the monotypic genus
  monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
  monotypic_genus_count <- length(monotypic_genus)
  #Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
  match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
  
  mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
  #Let's see whats common in both (mixture set)
  #intersect of id
  mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10
  mixtures_count <- length(mixtures)
  #Only in match_monotypic
  All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
  All_match_monotypic_count <- length(All_match_monotypic)
  #Only in mismatch_monotypic
  
  All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144
  All_mismatch_monotypic_count <- length(All_mismatch_monotypic)
  results <- data.frame(genus_mismatches=uniq,unique_mismatches_count,monotypic_genus=monotypic_genus_count,mixed_monotypic=mixtures_count,All_match_monotypic=All_match_monotypic_count,All_mismatch_monotypic=All_mismatch_monotypic_count)
  return(results)
}
#Sister species of placement sequences in reference tree#for bb99%
RefTree_R_co1_1_1_List <- foreach(i=1:length(species_List[[50]])) %do% sister(ReferenceTree,species_List[[50]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_99_co1_1_1_List <- foreach(i=1:length(species_List[[50]])) %do% sister(treeList$EPA99_9S_new.nwk,species_List[[50]][i],type="terminal",label=T)

#Call the function for bb99% and co1 1%
#For the N (order is,1= 99_1,2=99_10,3=99_2,4=99_3,5=99_4,6=99_5,7=99_6,8=99_7,9=99_8,10=99_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_1_1_List,R_99_co1_1_1_List,10)


RefTree_R_co1_1_1_List <- extract_Genus(RefTree_R_co1_1_1_List)
R_99_co1_1_1_List <- extract_Genus(R_99_co1_1_1_List)

nn <- union(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)
inter <- intersect(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)

mismatches <- anti_join(nn,inter)
#mismatches_1 <- setdiff(nn,inter)
mismatches <- setdiff(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)
count(unique_mismatches)
#Get the number of placement genus without matches to the placement genus in reference tree
count(unique(mismatches,by="id"))#212

#In this case since we have 904 co1 placements  in 80% tree the there are 692 matches
904-212

#for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
#First I get the genus name to particular the ids
mismatch_wt_names <- inner_join(mismatches,co1_1_id_list[[2]],by="id")#580 rows

#Get the number of placement genus without matches to the placement genus in reference tree
unique_mismatches <- unique(mismatch_wt_names,by="genus.y")#196
count(unique_mismatches)
#R_co1_20_1_G_vec col in unique_mismatches dataframe has the genus list of mismatches
#Using that col and set diff function get the monotypic genus

monotypic_genus <- setdiff(unique_mismatches$genus.y,Ref_Tree_distinct)#196
length(monotypic_genus)
#Now let's see the number of monotypic placement genus with clade genus that all matches, zero matches and mixtures (close enough monotypic)
match_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x==mismatch_wt_names$genus.y, ]#61
count(match_monotypic)
mismatch_monotypic <- mismatch_wt_names[mismatch_wt_names$genus.x!=mismatch_wt_names$genus.y, ]#519
#Let's see whats common in both (mixture set)
#intersect of id
mixtures <- intersect(match_monotypic$genus.y,mismatch_monotypic$genus.y)#10

#Only in match_monotypic
All_match_monotypic <- setdiff(match_monotypic$genus.y,mismatch_monotypic$genus.y)#42
#Only in mismatch_monotypic

All_mismatch_monotypic <- setdiff(mismatch_monotypic$genus.y,match_monotypic$genus.y)#144

