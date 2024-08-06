
#Last updated on Aug1, 2024
#All the files in stratified folder
#epang_stratified_genus
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
RefTree_R_co1_20_1_List <- foreach(i=1:length(species_List[[37]])) %do% sister(ReferenceTree,species_List[[37]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_80_co1_20_1_List <- foreach(i=1:length(species_List[[37]])) %do% sister(treeList$EPA80_6S_new.nwk,species_List[[37]][i],type="terminal",label=T)
#For the N (order is,1= 80_10,2=80_1,3=80_2,4=80_3,5=80_4,6=80_5,7=80_6,8=80_7,9=80_8,10=80_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_20_1_List,R_80_co1_20_1_List,7)
monotypic_genus_result
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
monotypic_genus_result 

#---------------
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
RefTree_R_co1_60_1_List <- foreach(i=1:length(species_List[[17]])) %do% sister(ReferenceTree,species_List[[17]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_40_co1_60_1_List <- foreach(i=1:length(species_List[[17]])) %do% sister(treeList$EPA40_6S_new.nwk,species_List[[17]][i],type="terminal",label=T)

#Call the function for bb40% and co1 60%
#For the N (order is,1= 40_10,2=40_1,3=40_2,4=40_3,5=40_4,6=40_5,7=40_6,8=40_7,9=40_8,10=40_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_60_1_List,R_40_co1_60_1_List,7)
monotypic_genus_result
#-----------
# df_RefTree <- extract_Genus(RefTree_R_co1_60_1_List)
# df_placement <- extract_Genus(R_40_co1_60_1_List)
# #get the rows of placement genus that's not matching with reference tree
# # mismatches <- setdiff(df_RefTree,df_placement)
# nn <- union(df_RefTree ,df_placement)
# inter <- intersect(df_RefTree ,df_placement)
# mismatches <- anti_join(nn,inter)
# #We need the number of genus not rows. So, get the number of placement genus that doesn't match to the placement genus in reference tree
# uniq <- count(unique(mismatches,by="id"))
# #for the downstream analysis, I need mismatches and check the genus with single species in mismatch genus list
# #First I get the genus name to particular the ids
# mismatch_wt_names <- inner_join(mismatches,co1_60_id_list[[7]],by="id")
# #Get the number of placement genus without matches to the placement genus in reference tree
# unique_mismatches <- unique(mismatch_wt_names,by="genus.y")
# unique_mismatches_count <- count(unique_mismatches)
#---------------
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
monotypic_genus_result
#------------------
# #get the number of species have similar sisters by comparing with reference tree
# RefTree_R_co1_80_1_List [(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #473
# #get the number of species doesn't have similar sisters by comparing with reference tree
# RefTree_R_co1_80_1_List [!(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #2810
#=---------

# RefTree_R_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(ReferenceTree,R_co1_80_1_species[i],type="terminal",label=T)
# R_20_co1_80_1_List <- foreach(i=1:length(R_co1_80_1_species)) %do% sister(R20_1_tree,R_co1_80_1_species[i],type="terminal",label=T)
# 
# #get the number of species have similar sisters by comparing with reference tree
# RefTree_R_co1_80_1_List [(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #102 #12.5%
# #get the number of species doesn't have similar sisters by comparing with reference tree
# RefTree_R_co1_80_1_List [!(RefTree_R_co1_80_1_List %in% R_20_co1_80_1_List )] #711

#-------

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
RefTree_R_co1_1_1_List <- foreach(i=1:length(species_List[[46]])) %do% sister(ReferenceTree,species_List[[46]][i],type="terminal",label=T)
#Sister species for placement sequeneces(20%) after placing on bb tree(80% bb tree)
R_99_co1_1_1_List <- foreach(i=1:length(species_List[[46]])) %do% sister(treeList$EPA99_5S_new.nwk,species_List[[46]][i],type="terminal",label=T)

#Call the function for bb99% and co1 1%
#For the N (order is,1= 99_1,2=99_10,3=99_2,4=99_3,5=99_4,6=99_5,7=99_6,8=99_7,9=99_8,10=99_9)
monotypic_genus_result <- genus_count(RefTree_R_co1_1_1_List,R_99_co1_1_1_List,6)
monotypic_genus_result

RefTree_R_co1_1_1_List <- extract_Genus(RefTree_R_co1_1_1_List)
R_99_co1_1_1_List <- extract_Genus(R_99_co1_1_1_List)

nn <- union(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)
inter <- intersect(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)
length(inter)
mismatches <- anti_join(nn,inter)
#mismatches_1 <- setdiff(nn,inter)
mismatches <- setdiff(RefTree_R_co1_1_1_List ,R_99_co1_1_1_List)
length(mismatches)
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

