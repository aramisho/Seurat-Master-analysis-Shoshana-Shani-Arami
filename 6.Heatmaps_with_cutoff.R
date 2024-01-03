library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)


load("~/cluster_jobs/UMAP_Objects.RData")





#        filtering and subseting 



HIP_seurat_list_log<- seurat_list_log[1]

all.genes_HIP <- rownames(HIP_seurat_list_log[[1]])

HIP_seurat_list_log[[1]] <- ScaleData(HIP_seurat_list_log[[1]], features = all.genes_HIP)
################################################################################################

MOs_FRP_seurat_list_log<- seurat_list_log[14]

all.genes_MOs_FRP <- rownames(MOs_FRP_seurat_list_log[[1]])

MOs_FRP_seurat_list_log[[1]] <- ScaleData(MOs_FRP_seurat_list_log[[1]], features = all.genes_MOs_FRP)
################################################################################################


MOp_seurat_list_log<- seurat_list_log[2]

all.genes_MOp <- rownames(MOp_seurat_list_log[[1]])

MOp_seurat_list_log[[1]] <- ScaleData(MOp_seurat_list_log[[1]], features = all.genes_MOp)




#RSP_seurat_list_log[[1]]@meta.data$number <- (RSP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

#RSP_seurat_list_log[[1]]@meta.data$CellNumber<- length(RSP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

#RSP_subnumber<- as.data.frame(table(RSP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

#colnames(RSP_subnumber)<- c("number","subclassCount")

#RSP_seurat_list_log[[1]]@meta.data<- merge(RSP_seurat_list_log[[1]]@meta.data,RSP_subnumber,by="number") 









HIP_seurat_list_log[[1]]@meta.data$number <- (HIP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

HIP_seurat_list_log[[1]]@meta.data$CellNumber<- length(HIP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

HIP_subnumber<- as.data.frame(table(HIP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(HIP_subnumber)<- c("number","subclassCount")





#HIP_seurat_list_log[[1]]@meta.data$number [HIP_seurat_list_log[[1]]@meta.data$subclass == as.character(HIP_subnumber$number[1])] <- HIP_subnumber$subclassCount[1]



NUMBER <- function(seuratobj,numberlist) {
  
  for (i in 1:length(numberlist$number)) {
    
    
    seuratobj[[1]]@meta.data$number [seuratobj[[1]]@meta.data$subclass == as.character(numberlist$number[i])] <- numberlist$subclassCount[i]
    
    
    meta<- seuratobj
    print(meta)
    #saveRDS(meta, file= paste( names(seuratobj[i]) ,"meta.rds"))
    save(meta, file= paste( names(seuratobj[i]) ,"meta.RData"))
    
  }
  
}


NUMBER(HIP_seurat_list_log,HIP_subnumber)













MOs_FRP_seurat_list_log[[1]]@meta.data$number <- (MOs_FRP_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOs_FRP_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOs_FRP_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOs_FRP_subnumber<- as.data.frame(table(MOs_FRP_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOs_FRP_subnumber)<- c("number","subclassCount")




MOp_seurat_list_log[[1]]@meta.data$number <- (MOp_seurat_list_log[[1]]@meta.data$subclass)# creating metadata column named number with subclass as rows

MOp_seurat_list_log[[1]]@meta.data$CellNumber<- length(MOp_seurat_list_log[[1]]@meta.data$orig.ident) #creating metadata coulmn named CellNumber with the total length of cells

MOp_subnumber<- as.data.frame(table(MOp_seurat_list_log[[1]]@meta.data$subclass))   #creating df from subclass name and number of cells 

colnames(MOp_subnumber)<- c("number","subclassCount")



#Idents(object = RSP_seurat_list_log[[1]])

Idents(object = MOs_FRP_seurat_list_log[[1]])

Idents(object = MOp_seurat_list_log[[1]])

Idents(object = HIP_seurat_list_log[[1]])




#Idents(object = RSP_seurat_list_log[[1]]) <- "subclass"



Idents(object = MOs_FRP_seurat_list_log[[1]]) <- "subclass"

Idents(object = MOp_seurat_list_log[[1]]) <- "subclass"

Idents(object = HIP_seurat_list_log[[1]]) <- "subclass"



MOs_FRP_subnumber
HIP_subnumber
MOp_subnumber


#RSP_subnumber




#as.list(subset(RSP_subnumber, subclassCount >100)["number" ] )  




as.list(subset(MOs_FRP_subnumber, subclassCount >100)["number" ] )  


as.list(subset(HIP_subnumber, subclassCount >100)["number" ] )  


as.list(subset(MOp_subnumber, subclassCount >100)["number" ] )  






#RSP_subsetOver100 <-  subset(x = RSP_seurat_list_log[[1]], idents =  c(  "Astro",  "L2/3 IT CTX", "L2/3 IT PPP" ,"L4 RSP-ACA" , 
#                                                                             "L4/5 IT CTX", "L5 IT CTX",   "L5 PT CTX",   "L5/6 NP CTX", "L6 CT CTX",  "L6 IT CTX",
#                                                                           "Lamp5","Oligo","Pvalb","Sncg","Sst",  "Vip" ))


MOs_FRPsubsetOver100 <-  subset(x = MOs_FRP_seurat_list_log[[1]], idents =  c(  "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX"  , 
                                                                                "L5/6 NP CTX" ,"L6 CT CTX" ,  "L6 IT CTX" , 
                                                                                "L6b CTX"   ,  "Lamp5"      , "Oligo"  ,   
                                                                                "Pvalb"   ,    "Sncg"    ,  "Sst"    ,
                                                                                "Vip"   ))

#MOs_FRPsubsetOver100 <- RenameIdents(object = MOs_FRPsubsetOver100, `L2/3 IT CTX` = "L2L3 IT CTX")
#MOs_FRPsubsetOver100 <- RenameIdents(object = MOs_FRPsubsetOver100, `L4/5 IT CTX` = "L4L5 IT CTX")
#MOs_FRPsubsetOver100 <- RenameIdents(object = MOs_FRPsubsetOver100, `L5/6 NP CTX` = "L5L6 NP CTX")


MOpsubsetOver100 <-  subset(x = MOp_seurat_list_log[[1]], idents =  c(  "Astro" ,      "Car3"     ,   "Endo"     , 
                                                                        "L2/3 IT CTX" ,"L2/3 IT PPP" ,"L4/5 IT CTX", "L5 IT CTX" , 
                                                                        "L5 PT CTX"  , "L5/6 NP CTX","L6 CT CTX"  , "L6 IT CTX"  ,
                                                                        "L6b CTX"  ,   "Lamp5"     ,  "Micro-PVM" ,  "Oligo"    ,  
                                                                        "Pvalb"    ,   "Sncg"    ,    "Sst"      ,   "Sst Chodl" , 
                                                                        "Vip"   ))

#MOpsubsetOver100 <- RenameIdents(object = MOpsubsetOver100, `L2/3 IT CTX` = "L2L3 IT CTX")
#MOpsubsetOver100 <- RenameIdents(object = MOpsubsetOver100, `L2/3 IT PPP` = "L2L3 IT PPP")
#MOpsubsetOver100 <- RenameIdents(object = MOpsubsetOver100, `L4/5 IT CTX` = "L4L5 IT CTX")
#MOpsubsetOver100 <- RenameIdents(object = MOpsubsetOver100, `L5/6 NP CTX` = "L5L6 NP CTX")
#MOpsubsetOver100 <- RenameIdents(object = MOpsubsetOver100, `Micro-PVM` = "Micro PVM")


HIP_FRPsubsetOver100 <-  subset(x = HIP_seurat_list_log[[1]], idents =  c(  "Astro","CA1-ProS" , "CA2-IG-FC" ,
                                                                            "CA3"   ,    "DG"    ,"Lamp5"   ,  
                                                                            "Oligo"    , "Sncg"     , "SUB-ProS" ,
                                                                            "Vip"    ))

#HIP_FRPsubsetOver100 <- RenameIdents(object = HIP_FRPsubsetOver100, `CA1-ProS` = "CA1 ProS")
#HIP_FRPsubsetOver100 <- RenameIdents(object = HIP_FRPsubsetOver100, `CA2-IG-FC` = "CA2 IG FC")
#HIP_FRPsubsetOver100 <- RenameIdents(object = HIP_FRPsubsetOver100, `SUB-ProS` = "SUB ProS")




#Idents(object = RSP_subsetOver100) 
#levels(x = RSP_subsetOver100)



Idents(object = MOs_FRPsubsetOver100) 
levels(x = MOs_FRPsubsetOver100)


Idents(object = MOpsubsetOver100) 
levels(x = MOpsubsetOver100)


Idents(object = HIP_FRPsubsetOver100) 
levels(x = HIP_FRPsubsetOver100)


listsubsetOver100<-c()

listsubsetOver100 <- append(listsubsetOver100, c( HIP=  HIP_FRPsubsetOver100 ,MOp= MOpsubsetOver100,MOs_FRP=MOs_FRPsubsetOver100))#, RSP= RSP_subsetOver100))




table(Idents(listsubsetOver100[[3]])) #MOs_FRP
MOs_FRPsubclassFreq<-as.data.frame(table(Idents(listsubsetOver100[[3]])))
colnames(MOs_FRPsubclassFreq)=c("subclass","Freq")
MOs_FRPsubclassFreq

table(Idents(listsubsetOver100[[2]]))#MOp
MOpsubclassFreq<-as.data.frame(table(Idents(listsubsetOver100[[2]])))
colnames(MOpsubclassFreq)=c("subclass","Freq")
MOpsubclassFreq

table(Idents(listsubsetOver100[[1]])) # HIP
HIPsubclassFreq<- as.data.frame(table(Idents(listsubsetOver100[[1]])) )
colnames(HIPsubclassFreq)=c("subclass","Freq")
HIPsubclassFreq



write.csv(HIPsubclassFreq,file = "HIPsubclassFreq.csv", row.names=FALSE)
write.csv(MOpsubclassFreq,file = "MOpsubclassFreq.csv", row.names=FALSE)
write.csv(MOs_FRPsubclassFreq,file = "MOs_FRPsubclassFreq.csv", row.names=FALSE)

#******* just_camkmt *********************************************************

just_camkmt <- function(list) {
  
  
  
  for (i in 1:length(list))
  {
    
    print(dim(list[[i]]))
    
    list[[i]] <- subset(x = list[[i]], subset = Camkmt > 0)
    
    print(dim(list[[i]]))
    
    
    
  }
  return(list)
}




seurat_list_just_camkmt<- just_camkmt(listsubsetOver100)

#********** PLOTING ******************

my_genes <- c("Camkmt","Calm1","Calm2","Calm3","Actb")

ploting_heatmap <- function(list, gene_list,number) {
  
  for (i in 1:length(list))
  {
    
    
      jpeg( paste("heatmap", " ",names(list)[i]," " ,".jpg"),width= 900, height =600 )
    
    
    Idents(object = list[[i]]) <- "subclass"
    
    
    g<- DoHeatmap(subset(list[[i]], downsample = number), 
                  features = my_genes, group.by = "subclass" ,slot='data' ,angle = 90, assay = "RNA",
                  draw.lines= FALSE,size = 3)  + scale_fill_distiller(palette = "YlGnBu") + #scale_fill_gradientn(colors = PurpleAndYellow(), na.value = "white") ,na.value = "white"
      
    
    ggtitle( "log expression heatmap ", subtitle = names(list)[i])+
      
      
      theme( 
        legend.direction="horizontal",legend.position="bottom", legend.title = element_text(size=12),legend.text = element_text(size=8),legend.key.size = unit(1, 'cm'),
        axis.text.y = element_text(hjust=1,vjust=0.5,size = 15),
        plot.title = element_text(hjust = 0.5,margin=margin(0,0,30,0),size = 20),
        plot.subtitle = element_text(hjust = 0.5,vjust = 9.0,margin=margin(0,0,1,0),size = 15)) +guides(color = FALSE, size = FALSE)#+guides(color = FALSE, size = FALSE)
    
    
    
    print( tail(g$data)  )
    
    x<- as.data.frame(g$data) 
    
    write.csv(x, file = paste("heatmap_log", names(list)[i], "heatmap",".csv"))
    
    print(g)
    dev.off()
    
    
    
  }
  
}


ploting_heatmap(listsubsetOver100[1],my_genes,number =3000 ) #HIP
ploting_heatmap(listsubsetOver100[2],my_genes,number= 1500) #MOp
ploting_heatmap(listsubsetOver100[3],my_genes,number=2300) #MOs_FRP

#******* just_camkmt ****************
ploting_heatmap(seurat_list_just_camkmt,my_genes, number = 2000)





