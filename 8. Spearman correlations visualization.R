library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library("tidyverse")
library("Hmisc")
library(ggh4x)
library(hrbrthemes)


# load data 
HIP_DG_BH_spearman_cutOff <- readRDS("~/cluster_jobs/HIP_DG_AdjustPvalues BH_spearman_cutOff.rds")
HIP_NO_DG_BH_spearman_cutOff <- readRDS("~/cluster_jobs/HIP_NO_DG_AdjustPvalues BH_spearman_cutOff.rds")
Mop_IT_CTX_BH_spearman_cutOff <- readRDS("~/cluster_jobs/Mop_IT_CTX_AdjustPvalues BH_spearman_cutOff.rds")
Mop_NO_IT_CTX_BH_spearman_cutOff <- readRDS("~/cluster_jobs/Mop_NO_IT_CTX_AdjustPvalues BH_spearman_cutOff.rds")
MOs_FRP_BH_spearman_cutOff <- readRDS("~/cluster_jobs/MOs_FRP_AdjustPvalues BH_spearman_cutOff.rds")



HIP_DG_BH_spearman_cutOff = mutate(HIP_DG_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 
HIP_NO_DG_BH_spearman_cutOff =mutate(HIP_NO_DG_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman)))
Mop_IT_CTX_BH_spearman_cutOff = mutate(Mop_IT_CTX_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 
Mop_NO_IT_CTX_BH_spearman_cutOff =  mutate(Mop_NO_IT_CTX_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman)))
MOs_FRP_BH_spearman_cutOff =mutate(MOs_FRP_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 


HIP_DG_BH_spearman_cutOff = mutate(HIP_DG_BH_spearman_cutOff,genes2 = fct_relevel(genes2, "Camkmt" , "Calm1" ,"Calm2" ,"Calm3","Actb" ))
HIP_NO_DG_BH_spearman_cutOff =mutate(HIP_NO_DG_BH_spearman_cutOff,genes2 = fct_relevel(genes2, "Camkmt" , "Calm1" ,"Calm2" ,"Calm3","Actb"))
Mop_IT_CTX_BH_spearman_cutOff = mutate(Mop_IT_CTX_BH_spearman_cutOff,genes2 = fct_relevel(genes2, "Camkmt" , "Calm1" ,"Calm2" ,"Calm3","Actb" ))
Mop_NO_IT_CTX_BH_spearman_cutOff = mutate(Mop_NO_IT_CTX_BH_spearman_cutOff,genes2 = fct_relevel(genes2, "Camkmt" , "Calm1" ,"Calm2" ,"Calm3","Actb" ))
MOs_FRP_BH_spearman_cutOff =mutate(MOs_FRP_BH_spearman_cutOff,genes2 = fct_relevel(genes2, "Camkmt" , "Calm1" ,"Calm2" ,"Calm3","Actb" ))



BH <- list("hippocampus DG cluster" = HIP_DG_BH_spearman_cutOff,"hippocampus clusters"= HIP_NO_DG_BH_spearman_cutOff,  
           "Primary motor area IT_CTX cluster " = Mop_IT_CTX_BH_spearman_cutOff, "Primary motor area clusters "= Mop_NO_IT_CTX_BH_spearman_cutOff, 
           "Secondary motor area and the frontal pole" = MOs_FRP_BH_spearman_cutOff)






ploting_coor <- function(list) {
  
  for (i in 1:length(list))
  {
    
    jpeg( paste("filterBH", " ",names(list)[i]," " ,".jpg"),width= 700, height =1500 )
    
    
    P =  ggplot(list[[i]] , aes(x = genes1, y =  genes2, fill =  spearman)) + 
      geom_tile()+
      geom_text(aes( x= genes1,y= genes2 , label =  round(spearman, digits = 3)  ), color = "black", size = 3) +
      scale_fill_distiller(palette = "RdYlBu", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) +
      theme_ipsum() +
      ggtitle( "spearman corralatiom heatmap", subtitle = names(list)[i]) +
      
      theme(
        plot.title = element_text(hjust = 0.5,margin=margin(0,0,50,0),size = 20),
        plot.subtitle = element_text(hjust = 0.5,vjust = 9.0,margin=margin(0,0,1,0),size = 18),
        axis.text.y = element_text(hjust=1,vjust=1.0,size = 10),
        axis.text.x = element_text(hjust=1,vjust=0.5,size = 13),
        legend.position="bottom", 
        
        legend.key.width= unit(2.5, 'cm') ,
        
        legend.title=element_blank(),
        
        axis.title.x = element_blank(),
        axis.title.y = element_blank()   )
    
    print(P)  
    
    dev.off()
  }
  
}


#runing functions


ploting_coor(list=BH )




Camkmt_HIP_DG_BH_spearman_cutOff <- readRDS("~/cluster_jobs/Camkmt HIP_DG_AdjustPvalues BH_spearman_cutOff.rds")
Camkmt_HIP_NO_DG_BH_spearman_cutOff <- readRDS("~/cluster_jobs/Camkmt HIP_NO_DG_AdjustPvalues BH_spearman_cutOff.rds")
Camkmt_Mop_IT_CTX_BH_spearman_cutOff  <- readRDS("~/cluster_jobs/Camkmt Mop_IT_CTX_AdjustPvalues BH_spearman_cutOff.rds")
Camkmt_Mop_NO_IT_CTX_BH_spearman_cutOff  <- readRDS("~/cluster_jobs/Camkmt Mop_NO_IT_CTX_AdjustPvalues BH_spearman_cutOff.rds")
Camkmt_MOs_FRP_BH_spearman_cutOff<- readRDS("~/cluster_jobs/Camkmt MOs_FRP_AdjustPvalues BH_spearman_cutOff.rds")


Camkmt_HIP_DG_BH_spearman_cutOff = mutate(Camkmt_HIP_DG_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 
Camkmt_HIP_NO_DG_BH_spearman_cutOff =mutate(Camkmt_HIP_NO_DG_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman)))
Camkmt_Mop_IT_CTX_BH_spearman_cutOff = mutate(Camkmt_Mop_IT_CTX_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 
Camkmt_Mop_NO_IT_CTX_BH_spearman_cutOff =  mutate(Camkmt_Mop_NO_IT_CTX_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman)))
Camkmt_MOs_FRP_BH_spearman_cutOff =mutate(Camkmt_MOs_FRP_BH_spearman_cutOff, genes2 = fct_reorder(genes2, desc(spearman))) 


Camkmt_BH <- list("hippocampus DG cluster " = Camkmt_HIP_DG_BH_spearman_cutOff,"hippocampus clusters"= Camkmt_HIP_NO_DG_BH_spearman_cutOff,  
                  "Primary motor area IT_CTX cluster " = Camkmt_Mop_IT_CTX_BH_spearman_cutOff, "Primary motor area clusters "= Camkmt_Mop_NO_IT_CTX_BH_spearman_cutOff, 
                  "Secondary motor area and the frontal pole" = Camkmt_MOs_FRP_BH_spearman_cutOff)




Camkmt_ploting_coor <- function(list) {
  
  for (i in 1:length(list))
  {
    
    jpeg( paste("Camkmt_filterBH", " ",names(list)[i]," " ,".jpg"),width= 700, height =1500 )
    
    
    P =  ggplot( list[[i]], aes(x = genes1, y =  genes2, fill =  spearman)) + 
      geom_tile()+
      geom_text(aes( x= genes1,y= genes2 , label =  round(spearman, digits = 3)  ), color = "black", size = 3) +
      scale_fill_distiller(palette = "RdYlBu", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) +
      theme_ipsum() +
      ggtitle( "spearman corralatiom heatmap - Camkmt  ", subtitle = names(list)[i]) +
      
      theme(
        plot.title = element_text(hjust = 0.5,margin=margin(0,0,50,0),size = 20),
        
        plot.subtitle = element_text(hjust = 0.5,vjust = 9.0,margin=margin(0,0,1,0),size = 18),
        
        axis.text.y = element_text(hjust=1,vjust=1.0,size = 10),
        
        axis.text.x = element_text(hjust=1,vjust=0.5,size = 13),
        
        legend.position="bottom", 
        
        legend.key.width= unit(2.5, 'cm') ,
        
        legend.title=element_blank(),
        
        
        axis.title.x = element_blank(),
        
        axis.title.y = element_blank()   ) 
    
    
    
    
    print(P)  
    
    dev.off()
  }
  
}


#runing functions


Camkmt_ploting_coor(list=Camkmt_BH )






