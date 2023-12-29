library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library("Hmisc")
library(tidyverse)

#Loading the normalized data
load("~/cluster_jobs/UMAP_Objects.RData")

seurat_list_log


my_genes <- c("Camkmt","Calm1","Calm2","Calm3","Actb")

#my_genes <- list(Camkmt="Camkmt",Calm1="Calm1",Calm2="Calm2",Calm3="Calm3",Actb="Actb")


# Focusing on three relevant brain regions: HIP ,MOs_FRP,MOp.  creat a totalof 100 cells per cluster objects


HIP_seurat_list_log<- seurat_list_log[1]


MOs_FRP_seurat_list_log<- seurat_list_log[14]


MOp_seurat_list_log<- seurat_list_log[2]




HIP_seurat_list_log[[1]]@meta.data$number <- (HIP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

HIP_seurat_list_log[[1]]@meta.data$CellNumber<- length(HIP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

HIP_subnumber<- as.data.frame(table(HIP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(HIP_subnumber)<- c("subclass","Count")



MOs_FRP_seurat_list_log[[1]]@meta.data$number <- (MOs_FRP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOs_FRP_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOs_FRP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOs_FRP_subnumber<- as.data.frame(table(MOs_FRP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOs_FRP_subnumber)<- c("subclass","Count")


MOp_seurat_list_log[[1]]@meta.data$number <- (MOp_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOp_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOp_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOp_subnumber<- as.data.frame(table(MOp_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOp_subnumber)<- c("subclass","Count")



#Idents(object = MOs_FRP_seurat_list_log[[1]])

#Idents(object = MOp_seurat_list_log[[1]])

#Idents(object = HIP_seurat_list_log[[1]])



Idents(object = MOs_FRP_seurat_list_log[[1]]) <- "subclass"

Idents(object = MOp_seurat_list_log[[1]]) <- "subclass"

Idents(object = HIP_seurat_list_log[[1]]) <- "subclass"



#MOs_FRP_subnumber
#HIP_subnumber
#MOp_subnumber



# 2 options forsubseting : only one neded:
#1

as.list(subset(MOs_FRP_subnumber, Count >100)["subclass" ] )  


as.list(subset(HIP_subnumber, Count >100)["subclass" ] )  


as.list(subset(MOp_subnumber, Count >100)["subclass" ] )  


MOs_FRPsubsetOver100 <-  subset(x = MOs_FRP_seurat_list_log[[1]], idents =  c(  "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX"  , 
                                                                                "L5/6 NP CTX" ,"L6 CT CTX" ,  "L6 IT CTX" , 
                                                                                "L6b CTX"   ,  "Lamp5"      , "Oligo"  ,   
                                                                                "Pvalb"   ,    "Sncg"    ,  "Sst"    ,
                                                                                "Vip"   ))



MOpsubsetOver100 <-  subset(x = MOp_seurat_list_log[[1]], idents =  c(  "Astro" ,      "Car3"     ,   "Endo"     , 
                                                                        "L2/3 IT CTX" ,"L2/3 IT PPP" ,"L4/5 IT CTX", "L5 IT CTX" , 
                                                                        "L5 PT CTX"  , "L5/6 NP CTX","L6 CT CTX"  , "L6 IT CTX"  ,
                                                                        "L6b CTX"  ,   "Lamp5"     ,  "Micro-PVM" ,  "Oligo"    ,  
                                                                        "Pvalb"    ,   "Sncg"    ,    "Sst"      ,   "Sst Chodl" , 
                                                                        "Vip"   ))


HIP_FRPsubsetOver100 <-  subset(x = HIP_seurat_list_log[[1]], idents =  c(  "Astro","CA1-ProS" , "CA2-IG-FC" ,
                                                                            "CA3"   ,    "DG"    ,"Lamp5"   ,  
                                                                            "Oligo"    , "Sncg"     , "SUB-ProS" ,
                                                                            "Vip"    ))
#..................................................................................................................................................
#2

as.list(subset(MOs_FRP_subnumber, Count <100)["subclass" ] )  


as.list(subset(HIP_subnumber, Count <100)["subclass" ] )  


as.list(subset(MOp_subnumber, Count <100)["subclass" ] )  


MOs_FRPsubsetOver100<- subset(x = MOs_FRP_seurat_list_log[[1]], idents =  c(  "Car3","L2/3 IT PPP", "L5 PT CTX","SMC-Peri","Sst Chodl"  ),invert=TRUE)

MOpsubsetOver100<- subset(x = MOp_seurat_list_log[[1]], idents =  c(  "CR","L5/6 IT TPE-ENT" ,"SMC-Peri", "VLMC"  ),invert=TRUE)


HIP_FRPsubsetOver100<- subset(x = HIP_seurat_list_log[[1]], idents =  c(  "CR","CT SUB" ,"Endo","L2/3 IT ENTl", "L2/3 IT PPP",
                                                                          "L2/3 IT RHP",  "L4/5 IT CTX",  "L5 PT CTX","L6 CT CTX","L6 IT CTX", 
                                                                          "L6b CTX","Micro-PVM" ,   "NP SUB","Pvalb","SMC-Peri","Sst"    
                                                                          ,"Sst Chodl" ,"VLMC"),invert=TRUE)
#..................................................................................................................................................
#..................................................................................................................................................

#cheking levels and idents:
Idents(object = MOs_FRPsubsetOver100) 
levels(x = MOs_FRPsubsetOver100)


Idents(object = MOpsubsetOver100) 
levels(x = MOpsubsetOver100)


Idents(object = HIP_FRPsubsetOver100) 
levels(x = HIP_FRPsubsetOver100)


#list of object with a cut of 100 cells per cluster
listsubsetOver100<-c()

listsubsetOver100 <- append(listsubsetOver100, c( HIP=  HIP_FRPsubsetOver100 ,MOp= MOpsubsetOver100,MOs_FRP=MOs_FRPsubsetOver100))

#********** ploting *****************************************************************************************************

ploting_expression_log_subNumber <- function(list, gene_list) {
  
  for (i in 1:length(list))
  {
    
    
    list[[i]]@meta.data$number <- (list[[i]]@meta.data$subclass)# creating metadata column named number with subclass as rows
    
    list[[i]]@meta.data$CellNumber<- length(list[[i]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells
    
    subnumber<- as.data.frame(table(list[[i]]@meta.data$subclass))   #creating df from subclass name and number of cells 
    
    colnames(subnumber)<- c("number","subclassCount")
    
    
    list[[i]]@meta.data <- merge(list[[i]]@meta.data,subnumber,by="number") #adding mew metadata column of subclass cell number 
    
    
    # merging all relevant column to 1 column
    list[[i]]@meta.data$class.ident <- paste(list[[i]]@meta.data$class, list[[i]]@meta.data$subclass,list[[i]]@meta.data$subclassCount,list[[i]]@meta.data$CellNumber ,sep = "_")
    
    
    
    jpeg( paste("expression_log", " ",names(list)[i]," " ,".jpg"),width= 800, height =600 )
    
    
    Idents(object = list[[i]]) <- "subclass"
    
    
    g<-DotPlot(list[[i]] , features = gene_list ,  col.min = -2.0, col.max = 2.0
               
               ,group.by	= "class.ident", cols = c("blue","yellow") ) + RotatedAxis()+ 
      
      ggtitle( "expression dotplot", subtitle = names(list)[i])+
      
      
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1.0,size = 10),
            axis.text.y = element_text(hjust=1,vjust=0.5,size = 10),
            plot.title = element_text(hjust = 0.5,margin=margin(0,0,30,0)),
            plot.subtitle = element_text(hjust = 0.5,vjust = 9.0,margin=margin(0,0,1,0))) 
    
    
    
    print(g)
    
    x<- as.data.frame(g$data)  
    print(head(x))
    write.csv(x, file = paste("expression_log", names(list)[i], "UMAP",".csv"))
    
    
    dev.off()
    
    
    
    
    
  }
  
}


ploting_expression_log_subNumber(listsubsetOver100,my_genes)


