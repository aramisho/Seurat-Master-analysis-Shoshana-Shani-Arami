library(data.table)
library("Hmisc")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)

#read csv files:

#1

ACA<-read.csv(file ="ACA.csv",row.names = NULL , header = TRUE) 

#2
AI<-read.csv(file ="AI.csv" ,row.names = "sample_name", header = TRUE )#

#3
AUD<-read.csv(file ="AUD.csv" ,row.names = "sample_name", header = TRUE)

#4
ENT<-read.csv(file ="ENT.csv" ,row.names = "sample_name", header = TRUE)

#5
HIP<-read.csv(file ="HIP.csv" ,row.names = "sample_name", header = TRUE)

#6
MOp<-read.csv(file ="MOp.csv",row.names = "sample_name", header = TRUE)

#7
MOs_FRP<-read.csv(file ="MOs_FRP.csv",row.names = "sample_name", header = TRUE)

#8
PAR_POST_PRE_SUB_ProS<-read.csv(file ="PAR-POST-PRE-SUB-ProS.csv" ,row.names = "sample_name", header = TRUE)
#9
PL_ILA_ORB<-read.csv(file ="PL-ILA-ORB.csv",row.names = "sample_name", header = TRUE)

#10
PTLp<-read.csv(file ="PTLp.csv" ,row.names = "sample_name", header = TRUE)

#11
RSP<-read.csv(file ="RSP.csv"  ,row.names = "sample_name", header = TRUE)

#12
SSp<-read.csv(file ="SSp.csv" ,row.names = "sample_name", header = TRUE)

#13
SSs_GU_VISC_AIp<-read.csv(file ="SSs-GU-VISC-AIp.csv" ,row.names = "sample_name", header = TRUE)

#14
TEa_PERI_ECT<-read.csv(file ="TEa-PERI-ECT.csv" ,row.names = "sample_name", header = TRUE)

#15
VIS<-read.csv(file ="VIS.csv",row.names = "sample_name", header = TRUE)

#16
VISl<-read.csv(file ="VISl.csv",row.names = "sample_name", header = TRUE)

#17
VISp<-read.csv(file ="VISp.csv",row.names = "sample_name", header = TRUE)

#18
VISm<-read.csv(file ="VISm.csv",row.names = "sample_name", header = TRUE)


#create  matadata or read: cluster_label:Cell type cluster name



matadata_ACA_class_label<- select(ACA,class_label)
matadata_ACA_cluster_label<- select(ACA,cluster_label)
matadata_ACA_region_label<- select(ACA,region_label)
matadata_ACA_subclass_label<- select(ACA,subclass_label)



matadata_AI_class_label<- select(AI,class_label)
matadata_AI_cluster_label<- select(AI,cluster_label)
matadata_AI_region_label<- select(AI,region_label)
matadata_AI_subclass_label<- select(AI,subclass_label)



matadata_AUD_class_label<- select(AUD,class_label)
matadata_AUD_cluster_label<- select(AUD,cluster_label)
matadata_AUD_region_label<- select(AUD,region_label)
matadata_AUD_subclass_label<- select(AUD,subclass_label)



matadata_ENT_class_label<- select(ENT,class_label)
matadata_ENT_cluster_label<- select(ENT,cluster_label)
matadata_ENT_region_label<- select(ENT,region_label)
matadata_ENT_subclass_label<- select(ENT,subclass_label)



matadata_HIP_class_label<- select(HIP,class_label)
matadata_HIP_cluster_label<- select(HIP,cluster_label)
matadata_HIP_region_label<- select(HIP,region_label)
matadata_HIP_subclass_label<- select(HIP,subclass_label)



matadata_MOp_class_label<- select(MOp,class_label)
matadata_MOp_cluster_label<- select(MOp,cluster_label)
matadata_MOp_region_label<- select(MOp,region_label)
matadata_MOp_subclass_label<- select(MOp,subclass_label)



matadata_MOs_FRP_class_label<- select(MOs_FRP,class_label)
matadata_MOs_FRP_cluster_label<- select(MOs_FRP,cluster_label)
matadata_MOs_FRP_region_label<- select(MOs_FRP,region_label)
matadata_MOs_FRP_subclass_label<- select(MOs_FRP,subclass_label)



matadata_PAR_POST_PRE_SUB_ProS_class_label<- select(PAR_POST_PRE_SUB_ProS,class_label)
matadata_PAR_POST_PRE_SUB_ProS_cluster_label<- select(PAR_POST_PRE_SUB_ProS,cluster_label)
matadata_PAR_POST_PRE_SUB_ProS_region_label<- select(PAR_POST_PRE_SUB_ProS,region_label)
matadata_PAR_POST_PRE_SUB_ProS_subclass_label<- select(PAR_POST_PRE_SUB_ProS,subclass_label)



matadata_PL_ILA_ORB_class_label<- select(PL_ILA_ORB,class_label)
matadata_PL_ILA_ORB_cluster_label<- select(PL_ILA_ORB,cluster_label)
matadata_PL_ILA_ORB_region_label<- select(PL_ILA_ORB,region_label)
matadata_PL_ILA_ORB_subclass_label<- select(PL_ILA_ORB,subclass_label)



matadata_PTLp_class_label<- select(PTLp,class_label)
matadata_PTLp_cluster_label<- select(PTLp,cluster_label)
matadata_PTLp_region_label<- select(PTLp,region_label)
matadata_PTLp_subclass_label<- select(PTLp,subclass_label)



matadata_RSP_class_label<- select(RSP,class_label)
matadata_RSP_cluster_label<- select(RSP,cluster_label)
matadata_RSP_region_label<- select(RSP,region_label)
matadata_RSP_subclass_label<- select(RSP,subclass_label)



matadata_SSp_class_label<- select(SSp,class_label)
matadata_SSp_cluster_label<- select(SSp,cluster_label)
matadata_SSp_region_label<- select(SSp,region_label)
matadata_SSp_subclass_label<- select(SSp,subclass_label)



matadata_SSs_GU_VISC_AIp_class_label<- select(SSs_GU_VISC_AIp,class_label)
matadata_SSs_GU_VISC_AIp_cluster_label<- select(SSs_GU_VISC_AIp,cluster_label)
matadata_SSs_GU_VISC_AIp_region_label<- select(SSs_GU_VISC_AIp,region_label)
matadata_SSs_GU_VISC_AIp_subclass_label<- select(SSs_GU_VISC_AIp,subclass_label)



matadata_TEa_PERI_ECT_class_label<- select(TEa_PERI_ECT,class_label)
matadata_TEa_PERI_ECT_cluster_label<- select(TEa_PERI_ECT,cluster_label)
matadata_TEa_PERI_ECT_region_label<- select(TEa_PERI_ECT,region_label)
matadata_TEa_PERI_ECT_subclass_label<- select(TEa_PERI_ECT,subclass_label)



matadata_VIS_class_label<- select(VIS,class_label)
matadata_VIS_cluster_label<- select(VIS,cluster_label)
matadata_VIS_region_label<- select(VIS,region_label)
matadata_VIS_subclass_label<- select(VIS,subclass_label)



matadata_VISl_class_label<- select(VISl,class_label)
matadata_VISl_cluster_label<- select(VISl,cluster_label)
matadata_VISl_region_label<- select(VISl,region_label)
matadata_VISl_subclass_label<- select(VISl,subclass_label)



matadata_VISm_class_label<- select(VISm,class_label)
matadata_VISm_cluster_label<- select(VISm,cluster_label)
matadata_VISm_region_label<- select(VISm,region_label)
matadata_VISm_subclass_label<- select(VISm,subclass_label)



matadata_VISp_class_label<- select(VISp,class_label)
matadata_VISp_cluster_label<- select(VISp,cluster_label)
matadata_VISp_region_label<- select(VISp,region_label)
matadata_VISp_subclass_label<- select(VISp,subclass_label)



#create data frames 

df_ACA <- subset (ACA, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_ACA_t<- t(df_ACA)
#df_ACA<-na.omit(df_ACA)

df_AI <- subset (AI, select=  -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_AI_t<- t(df_AI)
#df_AI<-na.omit(df_AI)

df_AUD <- subset (AUD, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_AUD_t<- t(df_AUD)
#df_AUD <-na.omit(df_AUD )

df_ENT <- subset (ENT, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_ENT_t<- t(df_ENT)
#df_ENT<-na.omit(df_ENT)

df_HIP<- subset (HIP, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_HIP_t<- t(df_HIP)
#df_HIP<-na.omit(df_HIP)

df_MOp <- subset (MOp, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_MOp_t<- t(df_MOp)
#df_MOp<-na.omit(df_MOp)

df_MOs_FRP <- subset (MOs_FRP, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_MOs_FRP_t<- t(df_MOs_FRP)
#df_MOs_FRP<-na.omit(df_MOs_FRP)

df_PAR_POST_PRE_SUB_ProS <- subset (PAR_POST_PRE_SUB_ProS , select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_PAR_POST_PRE_SUB_ProS_t<- t(df_PAR_POST_PRE_SUB_ProS)
#df_PAR_POST_PRE_SUB_ProS<-na.omit(df_PAR_POST_PRE_SUB_ProS)

df_PL_ILA_ORB <- subset (PL_ILA_ORB, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_PL_ILA_ORB_t<- t(df_PL_ILA_ORB)
#df_PL_ILA_ORB<-na.omit(df_PL_ILA_ORB)


df_PTLp <- subset (PTLp, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_PTLp_t<- t(df_PTLp)
#df_PTLp<-na.omit(df_PTLp)


df_RSP <- subset (RSP, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_RSP_t<- t(df_RSP)
#df_RSP<-na.omit(df_RSP)



df_SSp <- subset (SSp, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_SSp_t<- t(df_SSp)
#df_SSp<-na.omit(df_SSp)



df_SSs_GU_VISC_AIp <- subset (SSs_GU_VISC_AIp, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_SSs_GU_VISC_AIp_t<- t(df_SSs_GU_VISC_AIp)
#df_SSs_GU_VISC_AIp<-na.omit(df_SSs_GU_VISC_AIp)


df_TEa_PERI_ECT <- subset (TEa_PERI_ECT, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_TEa_PERI_ECT_t<- t(df_TEa_PERI_ECT)
#df_TEa_PERI_ECT<-na.omit(df_TEa_PERI_ECT)


df_VIS <- subset (VIS, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_VIS_t<- t(df_VIS)
#df_VIS<-na.omit(df_VIS)


df_VISl <- subset (VISl, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_VISl_t<- t(df_VISl)
#df_VISl<-na.omit(df_VISl)


df_VISp <- subset (VISp, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_VISp_t<- t(df_VISp)
#df_VISp<-na.omit(df_VISp)



df_VISm <- subset (VISm, select = -c(class_label,cluster_label,class_label,subclass_label,region_label) )
df_VISm_t<- t(df_VISm)
#df_VISm<-na.omit(df_VISm)


#create seurat objects

Object_ACA<-CreateSeuratObject(counts  = df_ACA_t, project="ACA", 
                               min.cells = 0, min.features = 200)


Object_ACA<-AddMetaData(object=Object_ACA,metadata=matadata_ACA_class_label,col.name="class")
Object_ACA<-AddMetaData(object=Object_ACA,metadata=matadata_ACA_cluster_label,col.name="cluster")
Object_ACA<-AddMetaData(object=Object_ACA,metadata=matadata_ACA_region_label,col.name="region_label")
Object_ACA<-AddMetaData(object=Object_ACA,metadata=matadata_ACA_subclass_label,col.name="subclass")



Object_AI<-CreateSeuratObject(counts  = df_AI_t, project="AI", 
                              min.cells = 0, min.features = 200)

Object_AI<-AddMetaData(object=Object_AI,metadata=matadata_AI_class_label,col.name="class")
Object_AI<-AddMetaData(object=Object_AI,metadata=matadata_AI_cluster_label,col.name="cluster")
Object_AI<-AddMetaData(object=Object_AI,metadata=matadata_AI_region_label,col.name="region_label")
Object_AI<-AddMetaData(object=Object_AI,metadata=matadata_AI_subclass_label,col.name="subclass")


Object_AUD<-CreateSeuratObject(counts  = df_AUD_t, project="AUD", 
                               min.cells = 0, min.features = 200)

Object_AUD<-AddMetaData(object=Object_AUD,metadata=matadata_AUD_class_label,col.name="class")
Object_AUD<-AddMetaData(object=Object_AUD,metadata=matadata_AUD_cluster_label,col.name="cluster")
Object_AUD<-AddMetaData(object=Object_AUD,metadata=matadata_AUD_region_label,col.name="region_label")
Object_AUD<-AddMetaData(object=Object_AUD,metadata=matadata_AUD_subclass_label,col.name="subclass")



Object_ENT<-CreateSeuratObject(counts  = df_ENT_t, project="ENT", 
                               min.cells = 0, min.features = 200)

Object_ENT<-AddMetaData(object=Object_ENT,metadata=matadata_ENT_class_label,col.name="class")
Object_ENT<-AddMetaData(object=Object_ENT,metadata=matadata_ENT_cluster_label,col.name="cluster")
Object_ENT<-AddMetaData(object=Object_ENT,metadata=matadata_ENT_region_label,col.name="region_label")
Object_ENT<-AddMetaData(object=Object_ENT,metadata=matadata_ENT_subclass_label,col.name="subclass")


Object_HIP<-CreateSeuratObject(counts  = df_HIP_t, project="HIP", 
                               min.cells = 0, min.features = 200)

Object_HIP<-AddMetaData(object=Object_HIP,metadata=matadata_HIP_class_label,col.name="class")
Object_HIP<-AddMetaData(object=Object_HIP,metadata=matadata_HIP_cluster_label,col.name="cluster")
Object_HIP<-AddMetaData(object=Object_HIP,metadata=matadata_HIP_region_label,col.name="region_label")
Object_HIP<-AddMetaData(object=Object_HIP,metadata=matadata_HIP_subclass_label,col.name="subclass")



Object_MOp<-CreateSeuratObject(counts  = df_MOp_t, project="MOp", 
                               min.cells = 0, min.features = 200)


Object_MOp<-AddMetaData(object=Object_MOp,metadata=matadata_MOp_class_label,col.name="class")
Object_MOp<-AddMetaData(object=Object_MOp,metadata=matadata_MOp_cluster_label,col.name="cluster")
Object_MOp<-AddMetaData(object=Object_MOp,metadata=matadata_MOp_region_label,col.name="region_label")
Object_MOp<-AddMetaData(object=Object_MOp,metadata=matadata_MOp_subclass_label,col.name="subclass")



Object_MOs_FRP<-CreateSeuratObject(counts  = df_MOs_FRP_t, project="MOs_FRP", 
                                   min.cells = 0, min.features = 200)


Object_MOs_FRP<-AddMetaData(object=Object_MOs_FRP,metadata=matadata_MOs_FRP_class_label,col.name="class")
Object_MOs_FRP<-AddMetaData(object=Object_MOs_FRP,metadata=matadata_MOs_FRP_cluster_label,col.name="cluster")
Object_MOs_FRP<-AddMetaData(object=Object_MOs_FRP,metadata=matadata_MOs_FRP_region_label,col.name="region_label")
Object_MOs_FRP<-AddMetaData(object=Object_MOs_FRP,metadata=matadata_MOs_FRP_subclass_label,col.name="subclass")



Object_PAR_POST_PRE_SUB_ProS<-CreateSeuratObject(counts  = df_PAR_POST_PRE_SUB_ProS_t, project="PAR_POST_PRE_SUB_ProS", 
                                                 min.cells = 0, min.features = 200)


Object_PAR_POST_PRE_SUB_ProS<-AddMetaData(object=Object_PAR_POST_PRE_SUB_ProS,metadata=matadata_PAR_POST_PRE_SUB_ProS_class_label,col.name="class")
Object_PAR_POST_PRE_SUB_ProS<-AddMetaData(object=Object_PAR_POST_PRE_SUB_ProS,metadata=matadata_PAR_POST_PRE_SUB_ProS_cluster_label,col.name="cluster")
Object_PAR_POST_PRE_SUB_ProS<-AddMetaData(object=Object_PAR_POST_PRE_SUB_ProS,metadata=matadata_PAR_POST_PRE_SUB_ProS_region_label,col.name="region_label")
Object_PAR_POST_PRE_SUB_ProS<-AddMetaData(object=Object_PAR_POST_PRE_SUB_ProS,metadata=matadata_PAR_POST_PRE_SUB_ProS_subclass_label,col.name="subclass")



Object_PL_ILA_ORB<-CreateSeuratObject(counts  = df_PL_ILA_ORB_t, project="PL_ILA_ORB", 
                                      min.cells = 0, min.features = 200)


Object_PL_ILA_ORB<-AddMetaData(object=Object_PL_ILA_ORB,metadata=matadata_PL_ILA_ORB_class_label,col.name="class")
Object_PL_ILA_ORB<-AddMetaData(object=Object_PL_ILA_ORB,metadata=matadata_PL_ILA_ORB_cluster_label,col.name="cluster")
Object_PL_ILA_ORB<-AddMetaData(object=Object_PL_ILA_ORB,metadata=matadata_PL_ILA_ORB_region_label,col.name="region_label")
Object_PL_ILA_ORB<-AddMetaData(object=Object_PL_ILA_ORB,metadata=matadata_PL_ILA_ORB_subclass_label,col.name="subclass")




Object_PTLp<-CreateSeuratObject(counts  = df_PTLp_t, project="PTLp", 
                                min.cells = 0, min.features = 200)


Object_PTLp<-AddMetaData(object=Object_PTLp,metadata=matadata_PTLp_class_label,col.name="class")
Object_PTLp<-AddMetaData(object=Object_PTLp,metadata=matadata_PTLp_cluster_label,col.name="cluster")
Object_PTLp<-AddMetaData(object=Object_PTLp,metadata=matadata_PTLp_region_label,col.name="region_label")
Object_PTLp<-AddMetaData(object=Object_PTLp,metadata=matadata_PTLp_subclass_label,col.name="subclass")



Object_RSP<-CreateSeuratObject(counts  = df_RSP_t, project="RSP", 
                               min.cells = 0, min.features = 200)


Object_RSP<-AddMetaData(object=Object_RSP,metadata=matadata_RSP_class_label,col.name="class")
Object_RSP<-AddMetaData(object=Object_RSP,metadata=matadata_RSP_cluster_label,col.name="cluster")
Object_RSP<-AddMetaData(object=Object_RSP,metadata=matadata_RSP_region_label,col.name="region_label")
Object_RSP<-AddMetaData(object=Object_RSP,metadata=matadata_RSP_subclass_label,col.name="subclass")




Object_SSp<-CreateSeuratObject(counts  = df_SSp_t, project="SSp", 
                               min.cells = 0, min.features = 200)


Object_SSp<-AddMetaData(object=Object_SSp,metadata=matadata_SSp_class_label,col.name="class")
Object_SSp<-AddMetaData(object=Object_SSp,metadata=matadata_SSp_cluster_label,col.name="cluster")
Object_SSp<-AddMetaData(object=Object_SSp,metadata=matadata_SSp_region_label,col.name="region_label")
Object_SSp<-AddMetaData(object=Object_SSp,metadata=matadata_SSp_subclass_label,col.name="subclass")


Object_SSs_GU_VISC_AIp<-CreateSeuratObject(counts  = df_SSs_GU_VISC_AIp_t, project="SSs;GU;VISC;AIp", 
                                           min.cells = 0, min.features = 200)


Object_SSs_GU_VISC_AIp<-AddMetaData(object=Object_SSs_GU_VISC_AIp,metadata=matadata_SSs_GU_VISC_AIp_class_label,col.name="class")
Object_SSs_GU_VISC_AIp<-AddMetaData(object=Object_SSs_GU_VISC_AIp,metadata=matadata_SSs_GU_VISC_AIp_cluster_label,col.name="cluster")
Object_SSs_GU_VISC_AIp<-AddMetaData(object=Object_SSs_GU_VISC_AIp,metadata=matadata_SSs_GU_VISC_AIp_region_label,col.name="region_label")
Object_SSs_GU_VISC_AIp<-AddMetaData(object=Object_SSs_GU_VISC_AIp,metadata=matadata_SSs_GU_VISC_AIp_subclass_label,col.name="subclass")



Object_TEa_PERI_ECT<-CreateSeuratObject(counts  = df_TEa_PERI_ECT_t, project="TEa_PERI_ECT", 
                                        min.cells = 0, min.features = 200)


Object_TEa_PERI_ECT<-AddMetaData(object=Object_TEa_PERI_ECT,metadata=matadata_TEa_PERI_ECT_class_label,col.name="class")
Object_TEa_PERI_ECT<-AddMetaData(object=Object_TEa_PERI_ECT,metadata=matadata_TEa_PERI_ECT_cluster_label,col.name="cluster")
Object_TEa_PERI_ECT<-AddMetaData(object=Object_TEa_PERI_ECT,metadata=matadata_TEa_PERI_ECT_region_label,col.name="region_label")
Object_TEa_PERI_ECT<-AddMetaData(object=Object_TEa_PERI_ECT,metadata=matadata_TEa_PERI_ECT_subclass_label,col.name="subclass")

Object_VIS<-CreateSeuratObject(counts  = df_VIS_t, project="VIS", 
                               min.cells = 0, min.features = 200)


Object_VIS<-AddMetaData(object=Object_VIS,metadata=matadata_VIS_class_label,col.name="class")
Object_VIS<-AddMetaData(object=Object_VIS,metadata=matadata_VIS_cluster_label,col.name="cluster")
Object_VIS<-AddMetaData(object=Object_VIS,metadata=matadata_VIS_region_label,col.name="region_label")
Object_VIS<-AddMetaData(object=Object_VIS,metadata=matadata_VIS_subclass_label,col.name="subclass")


Object_VISl<-CreateSeuratObject(counts  = df_VISl_t, project="VISl", 
                                min.cells = 0, min.features = 200)


Object_VISl<-AddMetaData(object=Object_VISl,metadata=matadata_VISl_class_label,col.name="class")
Object_VISl<-AddMetaData(object=Object_VISl,metadata=matadata_VISl_cluster_label,col.name="cluster")
Object_VISl<-AddMetaData(object=Object_VISl,metadata=matadata_VISl_region_label,col.name="region_label")
Object_VISl<-AddMetaData(object=Object_VISl,metadata=matadata_VISl_subclass_label,col.name="subclass")



Object_VISp<-CreateSeuratObject(counts  = df_VISp_t, project="VISp", 
                                min.cells = 0, min.features = 200)


Object_VISp<-AddMetaData(object=Object_VISp,metadata=matadata_VISp_class_label,col.name="class")
Object_VISp<-AddMetaData(object=Object_VISp,metadata=matadata_VISp_cluster_label,col.name="cluster")
Object_VISp<-AddMetaData(object=Object_VISp,metadata=matadata_VISp_region_label,col.name="region_label")
Object_VISp<-AddMetaData(object=Object_VISp,metadata=matadata_VISp_subclass_label,col.name="subclass")




Object_VISm<-CreateSeuratObject(counts  = df_VISm_t, project="VISm", 
                                min.cells = 0, min.features = 200)


Object_VISm<-AddMetaData(object=Object_VISm,metadata=matadata_VISm_class_label,col.name="class")
Object_VISm<-AddMetaData(object=Object_VISm,metadata=matadata_VISm_cluster_label,col.name="cluster")
Object_VISm<-AddMetaData(object=Object_VISm,metadata=matadata_VISm_region_label,col.name="region_label")
Object_VISm<-AddMetaData(object=Object_VISm,metadata=matadata_VISm_subclass_label,col.name="subclass")



#create empty list
seurat_list <- c()
#append seurat objects to list

seurat_list <- append(seurat_list, c(Object_ACA,Object_HIP,Object_SSs_GU_VISC_AIp,
                                     Object_MOp,Object_ENT,Object_AUD,Object_AI,Object_AI,Object_PAR_POST_PRE_SUB_ProS,
                                     Object_VISl,Object_VISp,Object_VIS,Object_VISm, Object_TEa_PERI_ECT,
                                     Object_SSp, Object_RSP, Object_PTLp, Object_PL_ILA_ORB,Object_MOs_FRP))



save(Object_ACA, file = "Object_ACA.RData")
save(Object_AI, file = "Object_AI.RData")
save(Object_AUD, file = "Object_AUD.RData")
save(Object_ENT, file = "Object_ENT.RData")
save(Object_HIP, file = "Object_HIP.RData")
save(Object_MOp, file = "Object_MOp.RData")
save(Object_MOs_FRP, file = "Object_MOs_FRP.RData")
save(Object_PAR_POST_PRE_SUB_ProS, file = "Object_PAR_POST_PRE_SUB_ProS.RData")
save(Object_PL_ILA_ORB, file = "Object_PL_ILA_ORB.RData")
save(Object_PTLp, file = "Object_PTLp.RData")
save(Object_RSP, file = "Object_RSP.RData")
save(Object_SSp, file = "Object_SSp.RData")
save(Object_SSs_GU_VISC_AIp, file = "Object_SSs_GU_VISC_AIp.RData")
save(Object_TEa_PERI_ECT, file = "Object_TEa_PERI_ECT.RData")
save(Object_VIS, file = "Object_VIS.RData")
save(Object_VISl, file = "Object_VISl.RData")
save(Object_VISp, file = "Object_VISp.RData")
save(Object_VISm, file = "Object_VISm.RData")

save(Object_ACA, Object_AI,Object_SSs_GU_VISC_AIp, Object_TEa_PERI_ECT,Object_SSp,
     Object_RSP,Object_PTLp,Object_AUD,Object_ENT,Object_HIP,Object_MOp,Object_MOs_FRP,
     Object_PAR_POST_PRE_SUB_ProS,Object_PL_ILA_ORB, file = "UMAP_Objects.RData")


