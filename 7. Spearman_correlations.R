library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)

# load normalized data 
load("~/cluster_jobs/UMAP_Objects.RData")


#        filtering and subseting 
#*************************************************************************************************************************************************


#RSP_seurat_list_log<- seurat_list_log[11]


HIP_seurat_list_log<- seurat_list_log[1]



MOs_FRP_seurat_list_log<- seurat_list_log[14]


MOp_seurat_list_log<- seurat_list_log[2]

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************




#RSP_seurat_list_log[[1]]@meta.data$number <- (RSP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

#RSP_seurat_list_log[[1]]@meta.data$CellNumber<- length(RSP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

#RSP_subnumber<- as.data.frame(table(RSP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

#colnames(RSP_subnumber)<- c("number","subclassCount")

#RSP_seurat_list_log[[1]]@meta.data<- merge(RSP_seurat_list_log[[1]]@meta.data,RSP_subnumber,by="number") 






#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************



#HIP_seurat_list_log[[1]]@meta.data$number <- (HIP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

HIP_seurat_list_log[[1]]@meta.data$CellNumber<- length(HIP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

HIP_subnumber<- as.data.frame(table(HIP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(HIP_subnumber)<- c("BrainRegion","subclassCount")

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************


#MOs_FRP_seurat_list_log[[1]]@meta.data$number <- (MOs_FRP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOs_FRP_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOs_FRP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOs_FRP_subnumber<- as.data.frame(table(MOs_FRP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOs_FRP_subnumber)<- c("BrainRegion","subclassCount")

#MOs_FRP_seurat_list_log[[1]]@meta.data<- merge(MOs_FRP_seurat_list_log[[1]]@meta.data,MOs_FRP_subnumber,by="number") 


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

#MOp_seurat_list_log[[1]]@meta.data$number <- (MOp_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOp_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOp_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOp_subnumber<- as.data.frame(table(MOp_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOp_subnumber)<- c("BrainRegion","subclassCount")

#MOp_seurat_list_log[[1]]@meta.data<- merge(MOp_seurat_list_log[[1]]@meta.data,MOp_subnumber,by="number") 


#*************************************************************************************************************************************************


#Idents(object = RSP_seurat_list_log[[1]])

Idents(object = MOs_FRP_seurat_list_log[[1]])

Idents(object = MOp_seurat_list_log[[1]])

Idents(object = HIP_seurat_list_log[[1]])




#Idents(object = RSP_seurat_list_log[[1]]) <- "subclass"


#*************************************************************************************************************************************************

Idents(object = MOs_FRP_seurat_list_log[[1]]) <- "subclass"

Idents(object = MOp_seurat_list_log[[1]]) <- "subclass"

Idents(object = HIP_seurat_list_log[[1]]) <- "subclass"
#*************************************************************************************************************************************************


#RSP_subnumber

MOs_FRP_subnumber
HIP_subnumber
MOp_subnumber



#as.list(subset(RSP_subnumber, subclassCount >100)["BrainRegion" ] )  


#*******************  Filter out clusters that have less than 100 cells  ******************************************************************************************************************************


as.list(subset(MOs_FRP_subnumber, subclassCount >100)["BrainRegion" ] )  


as.list(subset(HIP_subnumber, subclassCount >100)["BrainRegion" ] )  


as.list(subset(MOp_subnumber, subclassCount >100)["BrainRegion" ] )  



#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

#RSP_subnumber


#cutoff 100 to cell nunmber df

cutoffMOs_FRP_subnumber<-MOs_FRP_subnumber[MOs_FRP_subnumber$subclassCount >100, ]

cutoffHIP_subnumber <-HIP_subnumber[HIP_subnumber$subclassCount >100, ]

cutoffMOp_subnumber <- MOp_subnumber[MOp_subnumber$subclassCount >100, ]


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

#RSP_subsetOver100 <-  subset(x = RSP_seurat_list_log[[1]], idents =  c(  "Astro",  "L2/3 IT CTX", "L2/3 IT PPP" ,"L4 RSP-ACA" , 
#                                                                         "L4/5 IT CTX", "L5 IT CTX",   "L5 PT CTX",   "L5/6 NP CTX", "L6 CT CTX",  "L6 IT CTX",
#                                                                         "Lamp5","Oligo","Pvalb","Sncg","Sst",  "Vip" ))

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

MOs_FRPsubsetOver100 <-  subset(x = MOs_FRP_seurat_list_log[[1]], idents =  c(  "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX"  , 
                                                                                "L5/6 NP CTX" ,"L6 CT CTX" ,  "L6 IT CTX" , 
                                                                                "L6b CTX"   ,  "Lamp5"      , "Oligo"  ,   
                                                                                "Pvalb"   ,    "Sncg"    ,  "Sst"    ,
                                                                                "Vip"   ))


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

MOpsubsetOver100 <-  subset(x = MOp_seurat_list_log[[1]], idents =  c(  "Astro" ,      "Car3"     ,   "Endo"     , 
                                                                        "L2/3 IT CTX" ,"L2/3 IT PPP" ,"L4/5 IT CTX", "L5 IT CTX" , 
                                                                        "L5 PT CTX"  , "L5/6 NP CTX","L6 CT CTX"  , "L6 IT CTX"  ,
                                                                        "L6b CTX"  ,   "Lamp5"     ,  "Micro-PVM" ,  "Oligo"    ,  
                                                                        "Pvalb"    ,   "Sncg"    ,    "Sst"      ,   "Sst Chodl" , 
                                                                        "Vip"   ))

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

HIP_FRPsubsetOver100 <-  subset(x = HIP_seurat_list_log[[1]], idents =  c(  "Astro","CA1-ProS" , "CA2-IG-FC" ,
                                                                            "CA3"   ,    "DG"    ,"Lamp5"   ,  
                                                                            "Oligo"    , "Sncg"     , "SUB-ProS" ,
                                                                            "Vip"    ))

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

listsubsetOver100<-c()

listsubsetOver100 <- append(listsubsetOver100, c( HIP=  HIP_FRPsubsetOver100 ,MOp= MOpsubsetOver100,MOs_FRP=MOs_FRPsubsetOver100))#, RSP= RSP_subsetOver100))

saveRDS(listsubsetOver100,file = "Cutoff100.RData")


#*** shosh go to PandR its OK ***************


#************** Division of large clusters (hip and mop) for the purpose of analysis ************************************************************
listsubsetOver100<- read_rds("Cutoff100.RData")


table(Idents(listsubsetOver100[[1]]))
table(Idents(listsubsetOver100[[2]]))
table(Idents(listsubsetOver100[[3]]))

write.csv(table(Idents(listsubsetOver100[[1]])), "HIPCutoff100Table .csv")
write.csv(table(Idents(listsubsetOver100[[2]])),"MopCutoff100Table .csv")
write.csv(table(Idents(listsubsetOver100[[3]])), "MOs_FRPCutoff100Table .csv")


list_subclass<-unique(listsubsetOver100[[1]]@meta.data[["subclass"]]) 

LIST_noDG<-list_subclass[-4]



DG<-subset(x = listsubsetOver100[[1]], idents = "DG")

MINOS_DG<- subset(x = listsubsetOver100[[1]], idents = LIST_noDG)


listHIP<- list(DG=DG,MINOS_DG=MINOS_DG)

listHIP

test<-head(listHIP)

#*******************************  Mop  *********************************** 

Idents(object = listsubsetOver100[[2]]) <- "subclass"
table(Idents(listsubsetOver100[[2]]))


list_subclass<-unique(listsubsetOver100[[2]]@meta.data[["subclass"]])
#L6 IT CTX       L5 IT CTX     L2/3 IT CTX     L4/5 IT CTX 

LIST_noIT_CTX<-list_subclass[- c(8,9, 10,11)]



IT_CTX<-subset(x = listsubsetOver100[[2]], idents = c("L6 IT CTX","L5 IT CTX","L2/3 IT CTX","L4/5 IT CTX" ))

#unique(IT_CTX@meta.data[["subclass"]])

MINOS_IT_CTX<- subset(x = listsubsetOver100[[2]], idents = LIST_noIT_CTX)

#unique(MINOS_IT_CTX@meta.data[["subclass"]])


listMop<- list(IT_CTX=IT_CTX,MINOS_IT_CTX=MINOS_IT_CTX)


MOs_FRP<- listsubsetOver100[3]



###########  the standard deviation is zero ##############################################
# warning message in cor(X) : the standard deviation is zero” 
#The reason for this is that one of our  vector object contains only one value 
#for example 0 (or different nus identical value) expression of a gene across all cell [horizontal]
# or cell that all gens express the same value [vertical]

#*************************************************************************************************************************************************
#********************************** Create a reduced matrix containing my desired genes with Spearman correlation values ***************************************************************************************************************

coralationMatrix <- function(list) {
  
  for (i in 1:length(list)) {
    
    
    
    fivegeneMatrix= FetchData(object = list[[i]], vars = c("Camkmt","Calm1","Calm2","Calm3","Actb"))
    fivegeneMatrix<- as.matrix(fivegeneMatrix)
    print(dim(fivegeneMatrix))
    #print(fivegeneMatrix)
    
    allgenesMatrix <- t(as.matrix(list[[i]]@assays$RNA@data))
    print(dim(allgenesMatrix))
    #print(allgenesMatrix)
    
    print(nrow(fivegeneMatrix))
    print(nrow(allgenesMatrix))
    
    cor<-cor(fivegeneMatrix,allgenesMatrix,method ="spearman")
    
    saveRDS(cor, file= paste( names(list[i]) ,"cor.rds"))
    
    
  }
}

coralationMatrix(test)


coralationMatrix(listHIP)
coralationMatrix(listMop)
coralationMatrix(MOs_FRP)


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
HIP_DG <- readRDS("~/cluster_jobs/HIP_DG cor.rds")
HIP_NO_DG <- readRDS("~/cluster_jobs/HIP_MINOS_DG cor.rds")
Mop_IT_CTX <- readRDS("~/cluster_jobs/MOp_IT_CTX cor.rds")
Mop_NO_IT_CTX <- readRDS("~/cluster_jobs/MOp_MINOS_IT_CTX cor.rds")
MOs_FRP <- readRDS("~/cluster_jobs/MOs_FRP cor.rds")



#*************************************************************************************************************************************************

CellNumber <- data.frame (BrainRegion  = c("HIP_DG", "HIP_NO_DG", "Mop_IT_CTX","Mop_NO_IT_CTX","MOs_FRP"),
                          subclassCount = c( 17360, 54518,62811 ,59531,33792)
)


#cutoffMOs_FRP_subnumber<-MOs_FRP_subnumber[MOs_FRP_subnumber$subclassCount >100, ]

#cutoffHIP_subnumber <-HIP_subnumber[HIP_subnumber$subclassCount >100, ]

#cutoffMOp_subnumber <- MOp_subnumber[MOp_subnumber$subclassCount >100, ]

#L6 IT CTX 7093      L5 IT CTX 10531    L2/3 IT CTX 11844    L4/5 IT CTX 33343

##sum(cutoffMOp_subnumber$subclassCount)-(7093+10531+11844+33343)



#*************************************************************************************************************************************************

corlist<-list(HIP_DG=HIP_DG,  HIP_NO_DG=HIP_NO_DG, Mop_IT_CTX=Mop_IT_CTX,
              
              Mop_NO_IT_CTX=Mop_NO_IT_CTX, MOs_FRP=MOs_FRP )

#******************************* p value calculation ******************************************************************************************************************
#t = r√(n-2) / √(1-r2)   t-score
# p= pt(% above and under the t value ,amount sampels freedom degree =2 ,bolth tailes)            

for (i in 1:nrow(CellNumber))
{
  t <- corlist[[ CellNumber$BrainRegion[i] ]] * ((sqrt(CellNumber$subclassCount[i]) -2) / (1-(corlist[[CellNumber$BrainRegion[i]]]^2)))
  p <- pt(t,(CellNumber$subclassCount[i] -2),lower.tail = FALSE)
  #print(CellNumber$BrainRegion[i])
  saveRDS(p, file= paste(CellNumber$BrainRegion[i] ,"p_value.rds"))
}

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

#MOs_FRP
MOs_FRP_cor <- readRDS("~/cluster_jobs/MOs_FRP cor.rds")
MOs_FRP_p_value <- readRDS("~/cluster_jobs/MOs_FRP p_value.rds")

#Mop

Mop_IT_CTX_cor <- readRDS("~/cluster_jobs/MOp_IT_CTX cor.rds")
Mop_IT_CTX_p_value <- readRDS("~/cluster_jobs/Mop_IT_CTX p_value.rds")

Mop_NO_IT_CTX_cor <-  readRDS("~/cluster_jobs/MOp_MINOS_IT_CTX cor.rds")
Mop_NO_IT_CTX_p_value <- readRDS("~/cluster_jobs/Mop_NO_IT_CTX p_value.rds")

#HIP
HIP_DG_cor <-  readRDS("~/cluster_jobs/HIP_DG cor.rds")
HIP_DG_p_value <- readRDS("~/cluster_jobs/HIP_DG p_value.rds")

HIP_NO_DG_cor <-  readRDS("~/cluster_jobs/HIP_MINOS_DG cor.rds")
HIP_NO_DG_p_value <- readRDS("~/cluster_jobs/HIP_NO_DG p_value.rds")



#**************** melt df to create one table with p values and Spearman values *********************************************************************************************************************************


#MOs_FRP

spearman_MOs_FRP_melt<- reshape2::melt(MOs_FRP_cor,  value.name  ="spearman", varnames=c("genes1","genes2")) 
p_value_MOs_FRP_melt<- reshape2::melt(MOs_FRP_p_value,  value.name  ="p_value", varnames=c("genes1","genes2"))


#Mop
spearman_Mop_IT_CTX_melt<- reshape2::melt(Mop_IT_CTX_cor,  value.name  ="spearman", varnames=c("genes1","genes2")) 
p_value_Mop_IT_CTX_melt<- reshape2::melt(Mop_IT_CTX_p_value,  value.name  ="p_value", varnames=c("genes1","genes2"))

spearman_Mop_NO_IT_CTX_melt<- reshape2::melt(Mop_NO_IT_CTX_cor,  value.name  ="spearman", varnames=c("genes1","genes2")) 
p_value_Mop_NO_IT_CTX_melt<- reshape2::melt(Mop_NO_IT_CTX_p_value,  value.name  ="p_value", varnames=c("genes1","genes2"))

#HIP
spearman_HIP_DG_melt<- reshape2::melt(HIP_DG_cor,  value.name  ="spearman", varnames=c("genes1","genes2")) 
p_value_HIP_DG_melt<- reshape2::melt(HIP_DG_p_value,  value.name  ="p_value", varnames=c("genes1","genes2"))

spearman_HIP_NO_DG_melt<- reshape2::melt(HIP_NO_DG_cor,  value.name  ="spearman", varnames=c("genes1","genes2")) 
p_value_HIP_NO_DG_melt<- reshape2::melt(HIP_NO_DG_p_value,  value.name  ="p_value", varnames=c("genes1","genes2"))


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#merge
#*************************************************************************************************************************************************

#MOs_FRP
MOs_FRP=  merge(spearman_MOs_FRP_melt, p_value_MOs_FRP_melt, by=c("genes1","genes2"))   
#Remove rows with NA's using na.omit()
MOs_FRP <- na.omit(MOs_FRP)
saveRDS(MOs_FRP, file= "MOs_FRP_p&R.rds")
#*************************************************************************************************************************************************

#Mop
Mop_IT_CTX=  merge(spearman_Mop_IT_CTX_melt, p_value_Mop_IT_CTX_melt, by=c("genes1","genes2"))  
#Remove rows with NA's using na.omit()
Mop_IT_CTX <- na.omit(Mop_IT_CTX)
saveRDS(Mop_IT_CTX, file= "Mop_IT_CTX_p&R.rds")

Mop_NO_IT_CTX=  merge(spearman_Mop_NO_IT_CTX_melt, p_value_Mop_NO_IT_CTX_melt, by=c("genes1","genes2"))   
#Remove rows with NA's using na.omit()
Mop_NO_IT_CTX <- na.omit(Mop_NO_IT_CTX)
saveRDS(Mop_NO_IT_CTX, file= "Mop_NO_IT_CTX_p&R.rds")
#*************************************************************************************************************************************************

#HIP
HIP_DG=  merge(spearman_HIP_DG_melt, p_value_HIP_DG_melt, by=c("genes1","genes2"))   
#Remove rows with NA's using na.omit()
HIP_DG <- na.omit(HIP_DG)
saveRDS(HIP_DG, file= "HIP_DG_p&R.rds")

HIP_NO_DG=  merge(spearman_HIP_NO_DG_melt, p_value_HIP_NO_DG_melt, by=c("genes1","genes2")) 
#Remove rows with NA's using na.omit()
HIP_NO_DG <- na.omit(HIP_NO_DG)
saveRDS(HIP_NO_DG, file= "HIP_NO_DG_p&R.rds")

#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*                           Adjust P values
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************


HIP_DG_PandR <- readRDS("~/cluster_jobs/p&r/HIP_DG_p&R.rds")
HIP_NO_DG_PandR <- readRDS("~/cluster_jobs/p&r/HIP_NO_DG_p&R.rds")
Mop_IT_CTX_PandR <- readRDS("~/cluster_jobs/p&r/Mop_IT_CTX_p&R.rds")
Mop_NO_IT_CTX_PandR <- readRDS("~/cluster_jobs/p&r/Mop_NO_IT_CTX_p&R.rds")
MOs_FRP_PandR <- readRDS("~/cluster_jobs/p&r/MOs_FRP_p&R.rds")


list_melt<- list(HIP_DG_PandR = HIP_DG_PandR, HIP_NO_DG_PandR = HIP_NO_DG_PandR,
                 Mop_IT_CTX_PandR = Mop_IT_CTX_PandR, Mop_NO_IT_CTX_PandR = Mop_NO_IT_CTX_PandR,
                 MOs_FRP_PandR = MOs_FRP_PandR)

#*************************************************************************************************************************************************
#*************************************************************************************************************************************************


#test["bonferroni"]  <- p.adjust(  test["p_value"][,]  , method = "bonferroni", n = length(test["p_value"][,]))


#list[[i]]["bonferroni"]  <- p.adjust(  list[[i]]["p_value"][,]  , method = "bonferroni", n = length(list[[i]]["p_value"][,]))



AdjustPvalues <- function(list) {
  for (i in 1:length(list)) {
    
    
    list[[i]]["bonferroni"]  <- p.adjust(  list[[i]]["p_value"][,]  , method = "bonferroni", n = length(list[[i]]["p_value"][,]))
    
    list[[i]]["Benjamini&Hochberg"]  <- p.adjust(  list[[i]]["p_value"][,]  , method = "BH", n = length(list[[i]]["p_value"][,]))
    
    print( list[[i]])
    saveRDS(list[[i]], file= paste( names(list[i]) ,"AdjustPvalues.rds"))
  }
  
}


#AdjustPvalues(test)
AdjustPvalues(list_melt)



AdjustPvalues(list_melt)





HIP_DG_AdjustPvalues <- readRDS("~/cluster_jobs/HIP_DG_PandR AdjustPvalues.rds")
HIP_NO_DG_AdjustPvalues <- readRDS("~/cluster_jobs/HIP_NO_DG_PandR AdjustPvalues.rds")
Mop_IT_CTX_AdjustPvalues <- readRDS("~/cluster_jobs/Mop_IT_CTX_PandR AdjustPvalues.rds")
Mop_NO_IT_CTX_AdjustPvalues <- readRDS("~/cluster_jobs/Mop_NO_IT_CTX_PandR AdjustPvalues.rds")
MOs_FRP_AdjustPvalues <- readRDS("~/cluster_jobs/MOs_FRP_PandR AdjustPvalues.rds")

options(digits=10)



#filter the data : top and bottem 50 significant correlations

#********************* HIP_PandR ****************************************************************



spearman_cutOff_bh<-list(HIP_DG_AdjustPvalues = HIP_DG_AdjustPvalues,HIP_NO_DG_AdjustPvalues= HIP_NO_DG_AdjustPvalues,  
                         Mop_IT_CTX_AdjustPvalues = Mop_IT_CTX_AdjustPvalues,Mop_NO_IT_CTX_AdjustPvalues= Mop_NO_IT_CTX_AdjustPvalues, 
                         MOs_FRP_AdjustPvalues = MOs_FRP_AdjustPvalues)


spearman_cutOff_BH <- function(list) {
  
  for (i in 1:length(list))  {
    
    # print( list[[i]])
    print(names(list)[i])
    BH<- subset(list[[i]],  `Benjamini&Hochberg` <= 0.05) # FILTER none significant results 
    top_spearman<- top_n(BH, 50, spearman) # filter top spearman coroaltion
    bot_spearman<- top_n(BH, 50, -spearman) #filter bottom spearman coroaltion
    spearman_ranks<- rbind(top_spearman,bot_spearman) # join 2 dataframs into 1
    
    print(spearman_ranks)
    saveRDS(spearman_ranks, file= paste(names(list)[i] ,"BH_spearman_cutOff.rds"))
    
  }
}

spearman_cutOff_BH(spearman_cutOff_bh)



spearman_cutOff_bonf<-list(HIP_DG_AdjustPvalues = HIP_DG_AdjustPvalues,HIP_NO_DG_AdjustPvalues= HIP_NO_DG_AdjustPvalues,  
                           Mop_IT_CTX_AdjustPvalues = Mop_IT_CTX_AdjustPvalues,Mop_NO_IT_CTX_AdjustPvalues= Mop_NO_IT_CTX_AdjustPvalues, 
                           MOs_FRP_AdjustPvalues = MOs_FRP_AdjustPvalues)


spearman_cutOff_bonferroni <- function(list) {
  
  for (i in 1:length(list))  {
    
    # print( list[[i]])
    print(names(list)[i])
    BH<- subset(list[[i]],  bonferroni <= 0.05) # FILTER none significant results 
    top_spearman<- top_n(BH, 50, spearman) # filter top spearman coroaltion
    bot_spearman<- top_n(BH, 50, -spearman) #filter bottom spearman coroaltion
    spearman_ranks<- rbind(top_spearman,bot_spearman) # join 2 dataframs into 1
    
    print(spearman_ranks)
    saveRDS(spearman_ranks, file= paste(names(list)[i] ,"bonf_spearman_cutOff.rds"))
    
  }
}

spearman_cutOff_bonferroni(spearman_cutOff_bonf)



spearman_cutOff_bh

#############################################


spearman_cutOff_BH_CAMKMT <- function(list) {
  
  for (i in 1:length(list))  {
    
    # print( list[[i]])
    print(names(list)[i])
    BH<- subset(list[[i]],  genes1== "Camkmt"  &`Benjamini&Hochberg` <= 0.05) # FILTER none significant results 
    top_spearman<- top_n(BH, 50, spearman) # filter top spearman coroaltion
    bot_spearman<- top_n(BH, 50, -spearman) #filter bottom spearman coroaltion
    spearman_ranks<- rbind(top_spearman,bot_spearman) # join 2 dataframs into 1
    
    print(spearman_ranks)
    saveRDS(spearman_ranks, file= paste("Camkmt",names(list)[i] ,"BH_spearman_cutOff.rds"))
    
  }
}

spearman_cutOff_BH_CAMKMT(spearman_cutOff_bh)


spearman_cutOff_bonf

spearman_cutOff_bonferroni_CAMKMT <- function(list) {
  
  for (i in 1:length(list))  {
    
    # print( list[[i]])
    print(names(list)[i])
    BH<- subset(list[[i]],  genes1== "Camkmt"  & bonferroni <= 0.05) # FILTER none significant results 
    top_spearman<- top_n(BH, 50, spearman) # filter top spearman coroaltion
    bot_spearman<- top_n(BH, 50, -spearman) #filter bottom spearman coroaltion
    spearman_ranks<- rbind(top_spearman,bot_spearman) # join 2 dataframs into 1
    
    print(spearman_ranks)
    saveRDS(spearman_ranks, file= paste("Camkmt",names(list)[i] ,"bonf_spearman_cutOff.rds"))
    
  }
}


spearman_cutOff_bonferroni_CAMKMT(spearman_cutOff_bonf)


################################################

#saveRDS(HIP_DG_filterBh, file= ,"HIP_DG_filterBh.rds")
#saveRDS(HIP_DG_filterbonf , file= ,"HIP_DG_filterbonf.rds")

saveRDS(HIP_NO_DG_filterBh, file= ,"HIP_NO_DG_filterBh.rds")
saveRDS(HIP_NO_DG_filterbonf , file= ,"HIP_NO_DG_filterbonf.rds")
saveRDS(Mop_IT_CTX_filterBh, file= ,"Mop_IT_CTX_filterBh.rds")
saveRDS(Mop_IT_CTX_filterbonf , file= ,"Mop_IT_CTX_filterbonf.rds")

saveRDS(Mop_NO_IT_CTX_filterBh, file= ,"Mop_NO_IT_CTX_filterBh.rds")
saveRDS(Mop_NO_IT_CTX_filterbonf , file= ,"Mop_NO_IT_CTX_filterbonf.rds")

saveRDS(MOs_FRP_filterBh, file= ,"MOs_FRP_filterBh.rds")
saveRDS(MOs_FRP_filterbonf , file= ,"MOs_FRP_filterbonf.rds")
################################################

HIP_DG_filterBh <- readRDS("~/cluster_jobs/HIP_DG_filterBh.rds")
HIP_DG_filterbonf <- readRDS("~/cluster_jobs/HIP_DG_filterbonf.rds")

HIP_NO_DG_filterBh <- readRDS("~/cluster_jobs/HIP_NO_DG_filterBh.rds")
HIP_NO_DG_filterbonf <- readRDS("~/cluster_jobs/HIP_NO_DG_filterbonf.rds")

Mop_IT_CTX_filterBh <- readRDS("~/cluster_jobs/Mop_IT_CTX_filterBh.rds")
Mop_IT_CTX_filterbonf <- readRDS("~/cluster_jobs/Mop_IT_CTX_filterbonf.rds")

Mop_NO_IT_CTX_filterBh <- readRDS("~/cluster_jobs/Mop_NO_IT_CTX_filterBh.rds")
Mop_NO_IT_CTX_filterbonf <- readRDS("~/cluster_jobs/Mop_NO_IT_CTX_filterbonf.rds")

MOs_FRP_filterBh <- readRDS("~/cluster_jobs/MOs_FRP_filterBh.rds")
MOs_FRP_filterbonf <- readRDS("~/cluster_jobs/MOs_FRP_filterbonf.rds")




################################################
########### csv saving 
#####################################
################################################


write.csv(HIP_DG_filterBh, file= ,"HIP_DG_filterBh.csv")
write.csv(HIP_DG_filterbonf , file= ,"HIP_DG_filterbonf.csv")

write.csv(HIP_NO_DG_filterBh, file= ,"HIP_NO_DG_filterBh.csv")
write.csv(HIP_NO_DG_filterbonf , file= ,"HIP_NO_DG_filterbonf.csv")


write.csv(Mop_IT_CTX_filterBh, file= ,"Mop_IT_CTX_filterBh.csv")
write.csv(Mop_IT_CTX_filterbonf , file= ,"Mop_IT_CTX_filterbonf.csv")

write.csv(Mop_NO_IT_CTX_filterBh, file= ,"Mop_NO_IT_CTX_filterBh.csv")
write.csv(Mop_NO_IT_CTX_filterbonf , file= ,"Mop_NO_IT_CTX_filterbonf.csv")

write.csv(MOs_FRP_filterBh, file= ,"MOs_FRP_filterBh.csv")
write.csv(MOs_FRP_filterbonf , file= ,"MOs_FRP_filterbonf.csv")
###

