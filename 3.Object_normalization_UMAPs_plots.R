library(data.table)
library("Hmisc")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)

#Loading and repeating the previous step
#No need to run previous steps again

load("~/cluster_jobs/NEWsplits/Object_ACA.RData") 
load("~/cluster_jobs/NEWsplits/Object_AI.RData")
load("~/cluster_jobs/NEWsplits/Object_AUD.RData")
load("~/cluster_jobs/NEWsplits/Object_ENT.RData")
load("~/cluster_jobs/NEWsplits/Object_HIP.RData")
load("~/cluster_jobs/NEWsplits/Object_MOp.RData")
load("~/cluster_jobs/NEWsplits/Object_MOs_FRP.RData")
load("~/cluster_jobs/NEWsplits/Object_PAR_POST_PRE_SUB_ProS.RData")
load("~/cluster_jobs/NEWsplits/Object_PL_ILA_ORB.RData")
load("~/cluster_jobs/NEWsplits/Object_PTLp.RData")
load("~/cluster_jobs/NEWsplits/Object_RSP.RData")
load("~/cluster_jobs/NEWsplits/Object_SSp.RData")
load("~/cluster_jobs/NEWsplits/Object_SSs_GU_VISC_AIp.RData")
load("~/cluster_jobs/NEWsplits/Object_TEa_PERI_ECT.RData")

save(Object_ACA, Object_AI,Object_SSs_GU_VISC_AIp, Object_TEa_PERI_ECT,Object_SSp,
     Object_RSP,Object_PTLp,Object_AUD,Object_ENT,Object_HIP,Object_MOp,Object_MOs_FRP,
     Object_PAR_POST_PRE_SUB_ProS,Object_PL_ILA_ORB, file = "UMAP_Objects.RData")


#creat seurat objects list 
seurat_list <- c()

seurat_list <- append(seurat_list, c( ACA=Object_ACA ,HIP= Object_HIP   ,
                                      
                                      SSs_GU_VISC_AIp= Object_SSs_GU_VISC_AIp ,
                                      
                                      MOp= Object_MOp  ,ENT= Object_ENT  ,AUD= Object_AUD ,
                                      
                                      AI= Object_AI, 
                                      
                                      PAR_POST_PRE_SUB_ProS= Object_PAR_POST_PRE_SUB_ProS ,#VISl = Object_VISl ,
                                      
                                      # VISp = Object_VISp  ,VIS= Object_VIS  ,VISm= Object_VISm  ,
                                      
                                      TEa_PERI_ECT= Object_TEa_PERI_ECT  ,SSp= Object_SSp  ,RSP= Object_RSP ,
                                      
                                      PTLp= Object_PTLp , PL_ILA_ORB = Object_PL_ILA_ORB , MOs_FRP= Object_MOs_FRP))

#or
load("~/cluster_jobs/UMAP_Objects.RData")
#Normalize the Data

seurat_list_log<- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method ="LogNormalize",scale.factor = 1e6 ,verbose = TRUE)
  x <- FindVariableFeatures(x,selection.method= "mean.var.plot",verbose = TRUE)
})

#creat list of interes genes 

my_genes <- c("Camkmt","Calm1","Calm2","Calm3","Actb")

#Select genes as anchors

features_log <- SelectIntegrationFeatures(object.list =seurat_list_log ,  nfeatures = 2000,assay = NULL,verbose = TRUE)
features_log <- append(features_log,my_genes)

# Scaling, dimensionality reduction analysis
seurat_list_log <- lapply(X =seurat_list_log, FUN = function(x) {
  x <- ScaleData(x, features = features_log, verbose = TRUE)
  x <- RunPCA(x, features = features_log, verbose = TRUE)
  x <- RunUMAP(x, dims = 1:50)
})


#################################################  ploting ############################


#.....................subclass...........organized by class.........................

ploting_expression_log_subclass_split <- function(list) {
  
  for (i in 1:length(list))
  {
    
    jpeg( paste(paste( "subclass of Brain region organized by class" ,levels(x =list[[i]]) ) ,".jpg"),width= 800, height =600 )
    
    
    
    
    g<-DimPlot(list[[i]],raster=FALSE, group.by = "subclass", label = T,repel = T,split.by="class")+ ggtitle(paste( "subclass of Brain region organized by class" ,levels(x =list[[i]]) ))
    g  +  FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) + NoLegend()
    
    
    
    print(g)
    dev.off()
    

  }
  
}


ploting_expression_log_subclass_split(seurat_list_log)



#.....................subclass........NO spliting ............................

ploting_expression_log_subclass <- function(list) {
  
  for (i in 1:length(list))
  {
    
    jpeg( paste(paste( "subclass of Brain region" ,levels(x =list[[i]]) ) ,".jpg"),width= 800, height =600 )
    
    
    
    
    g<-DimPlot(list[[i]],raster=FALSE, group.by = "subclass", label = T,repel = T)+ ggtitle(paste( "subclass of Brain region" ,levels(x =list[[i]]) ))
    g  +  FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) + NoLegend()
    
    
    
    print(g)
    dev.off()
    
    
  }
  
}


ploting_expression_log_subclass (seurat_list_log)



#.....................class....................................


ploting_expression_log_class <- function(list) {
  
  for (i in 1:length(list))
  {
    
    jpeg( paste(paste( "class of Brain region" ,levels(x =list[[i]]) ) ,".jpg"),width= 800, height =600 )
    
    
    
    
    g<-DimPlot(list[[i]],raster=FALSE, group.by = "class", label = T,repel = T)+ ggtitle(paste( "class of Brain region" ,levels(x =list[[i]]) ))
    g  +  FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) + NoLegend()
    
    
    
    print(g)
    dev.off()
    

    
  }
  
}


ploting_expression_log_class(seurat_list_log)


















#  ,repel = T  legend.position = "none",


plot <- lapply(X = seurat_list_log, FUN = function(x) {
  
  x <-DimPlot(x,raster=FALSE, group.by = "subclass", label = T,repel = T,split.by="class")  + ggtitle(paste( "subclass of Brain region" ,levels(x) ))
  
  
  x + FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) +theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                                                    legend.key.height = unit(0.3, 'cm'), #change legend key height
                                                                                    legend.key.width = unit(0.5, 'cm'), #change legend key width
                                                                                    legend.title = element_text(size=10), #change legend title font size
                                                                                    legend.text = element_text(size=8)) #change legend text font size
  
  
  print(x)
  
})

######################################################################################theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
legend.key.height = unit(0.3, 'cm'), #change legend key height
legend.key.width = unit(0.5, 'cm'), #change legend key width
legend.title = element_text(size=10), #change legend title font size
legend.text = element_text(size=8)) #change legend text font size
#####################################################################################





# 
jpeg( paste("TESTTTT:",".jpg"),width= 1000, height =600 )



rc_plot_sub <-DimPlot(seurat_list_rc[[3]],raster=FALSE, group.by = "subclass", repel = T,label=TRUE,label.size=4)+ ggtitle(paste( "subclass of Brain region" ,levels(x =seurat_list_rc[[3]]) ))
rc_plot_sub + FontSize(x.title = 15, y.title = 15 , x.text = 13,  y.text = 13,main = 20) +theme(legend.key.size = unit(100, 'cm'), #change legend key size
                                                                                                legend.key.height = unit(1.1, 'cm'), #change legend key height
                                                                                                legend.key.width = unit(0.1, 'cm'), #change legend key width
                                                                                                legend.title = element_text(size=10), #change legend title font size
                                                                                                legend.text = element_text(size=10), #change legend text font size
                                                                                                legend.spacing.y = unit(0.2, "cm"))

rc_plot_sub+ guides(col = guide_legend(nrow = 1))
dev.off()






rc_plot_sub<-DimPlot(seurat_list_rc[[1]],raster=FALSE, group.by = "cluster")
rc_plot_sub + FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) +theme(legend.key.size = unit(0.00000001, 'cm'), #change legend key size
                                                                                            legend.key.height = unit(0.3, 'cm'), #change legend key height
                                                                                            legend.key.width = unit(0.5, 'cm'), #change legend key width
                                                                                            legend.title = element_text(size=10), #change legend title font size
                                                                                            legend.text = element_text(size=8)) #change legend text font size

############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
############## chosennnnnnnnnnnnnnnnnnnnnn
#***************log+mpv*********************
#*
#*

seurat_list_log<- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method ="LogNormalize",scale.factor = 1e6 ,verbose = TRUE)
  x <- FindVariableFeatures(x,selection.method= "mean.var.plot",verbose = TRUE)
})


features_log <- SelectIntegrationFeatures(object.list =seurat_list_log ,  nfeatures = 2000,assay = NULL,verbose = TRUE)

seurat_list_log <- lapply(X =seurat_list_log, FUN = function(x) {
  x <- ScaleData(x, features = features_log, verbose = TRUE)
  x <- RunPCA(x, features = features_log, verbose = TRUE)
})
#


anchors_log <- FindIntegrationAnchors(object.list = seurat_list_log, reference = NULL, reduction = "rpca",dims = 1:50)


#features_log <- SelectIntegrationFeatures(object.list =seurat_list_log ,  nfeatures = 2000,assay = NULL,verbose = TRUE)

#seurat_list_log <- lapply(X =seurat_list_log, FUN = function(x) {
# x <- ScaleData(x, features = features_log, verbose = TRUE)
#  x <- RunPCA(x, features = features_log, verbose = TRUE)
#})
#


#anchors_log <- FindIntegrationAnchors(object.list = seurat_list_log, reference = NULL, reduction = "rpca",dims = 1:50)



#list.integrated_log <- IntegrateData(anchorset = anchors_log, dims = 1:50)

#list.integrated_log <- ScaleData(list.integrated_log, verbose = TRUE)

#list.integrated_log <- RunCCA(list.integrated_log, verbose = TRUE)

list.integrated_log <- RunPCA(list.integrated_log, verbose = TRUE)
list.integrated_log <- RunUMAP(list.integrated_log, dims = 1:50)

#save(list.integrated_rc,file= "list.integrated_log.RData")


# ploting 

log_plot<-DimPlot(list.integrated_log,raster=FALSE, group.by = "region_label")
log_plot + FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) +theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                                                         legend.key.height = unit(0.3, 'cm'), #change legend key height
                                                                                         legend.key.width = unit(0.5, 'cm'), #change legend key width
                                                                                         legend.title = element_text(size=10), #change legend title font size
                                                                                         legend.text = element_text(size=8)) #change legend text font size


log_plot<-DimPlot(list.integrated_log,raster=FALSE, group.by = "subclass")
log_plot + FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) +theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                                                         legend.key.height = unit(0.3, 'cm'), #change legend key height
                                                                                         legend.key.width = unit(0.5, 'cm'), #change legend key width
                                                                                         legend.title = element_text(size=10), #change legend title font size
                                                                                         legend.text = element_text(size=8)) #change legend text font size




log_plot<-DimPlot(list.integrated_log,raster=FALSE, group.by = "class")
log_plot + FontSize(x.title = 7, y.title = 7 , x.text = 7,  y.text = 5,main = 10) +theme(legend.key.size = unit(3, 'cm'), #change legend key size
                                                                                         legend.key.height = unit(0.3, 'cm'), #change legend key height
                                                                                         legend.key.width = unit(0.5, 'cm'), #change legend key width
                                                                                         legend.title = element_text(size=10), #change legend title font size
                                                                                         legend.text = element_text(size=8)) #change legend text font size




