library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Hmisc)


load("~/cluster_jobs/UMAP_Objects.RData")

my_genes <- c("Camkmt","Calm1","Calm2","Calm3","Actb")

seurat_list<- seurat_list_log

##############  Feature Plot  ###################################################


for (i in 1:length(seurat_list)) {
  
  Idents(object = seurat_list[[i]]) <- "region_label"
  
  jpeg( paste("region_label:",levels(seurat_list[[i]]),".jpg"),width= 1500, height =400 )
  
  
  Idents(object = seurat_list[[i]]) <- "subclass"
  
  location_plot<- FeaturePlot(seurat_list[[i]], features =c("Camkmt","Calm1","Calm2","Calm3","Actb")
                              # ,cols = c("aquamarine","blue4")
                              # ,cols = rev(brewer.pal(n =11 , name = "RdYlBu"))
                              , raster=FALSE ,
                              keep.scale="all",ncol = 5, label=TRUE,label.size= 2.5,repel = TRUE,min.cutoff = 0) #,split.by="class") 
  
  
  
  Idents(object = seurat_list[[i]]) <- "region_label"
  print(location_plot)
  
  
  #print(head(as.data.frame(location_plot[[1]]$data)))
  #print(head(as.data.frame(location_plot[[2]]$data)))
  #print(head(as.data.frame(location_plot[[3]]$data)))
  #print(head(as.data.frame(location_plot[[4]]$data)))
  #print(head(as.data.frame(location_plot[[5]]$data)))
  
  Camkmt <-(as.data.frame(location_plot[[1]]$data))
  Calm1 <- (as.data.frame(location_plot[[2]]$data))
  Calm3 <- (as.data.frame(location_plot[[3]]$data))
  Calm2 <- (as.data.frame(location_plot[[4]]$data))
  Actb <-  (as.data.frame(location_plot[[5]]$data))
  
  write.csv(Camkmt, file = paste("location_log", names(seurat_list)[i], "UMAPCamkmt",".csv"))
  write.csv(Calm1, file = paste("location_log", names(seurat_list)[i], "UMAPCalm1",".csv"))
  write.csv(Calm2, file = paste("location_log", names(seurat_list)[i], "UMAPCalm2",".csv"))
  write.csv(Calm3, file = paste("location_log", names(seurat_list)[i], "UMAPCalm3",".csv"))
  write.csv(Actb, file = paste("location_log", names(seurat_list)[i], "UMAPActb",".csv"))
  
  Idents(object = seurat_list[[i]]) <- "subclass"
  print(as.data.frame(AverageExpression(object = seurat_list[[i]], features= c("Camkmt","Calm1","Calm2","Calm3","Actb"), slot = "data")))
  AverageExpression <-(as.data.frame(AverageExpression(object = seurat_list[[i]], features= c("Camkmt","Calm1","Calm2","Calm3","Actb"), slot = "data")))
  write.csv(AverageExpression, file = paste("AverageExpression", names(seurat_list)[i], "UMAPAverageExpression",".csv"))
  
  
  Idents(object = seurat_list[[i]]) <- "orig.ident"
  print(as.data.frame(AverageExpression(object = seurat_list[[i]], features= c("Camkmt","Calm1","Calm2","Calm3","Actb"), slot = "data")))
  All_AverageExpression <- (as.data.frame(AverageExpression(object = seurat_list[[i]], features= c("Camkmt","Calm1","Calm2","Calm3","Actb"), slot = "data")))
  write.csv(All_AverageExpression, file = paste("All_AverageExpression", names(seurat_list)[i], "UMAPAll_AverageExpression",".csv"))
  
  
  dev.off()
}

















