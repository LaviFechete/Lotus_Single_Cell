#Data preprocessing and integration 5dpi datasets####

library(stringr)
library(dplyr)
suppressPackageStartupMessages(library(Seurat))
library(Matrix)
library(patchwork)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)
suppressMessages(require(SeuratWrappers))
library(metap)
library(plotly)
library(scater)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(VennDiagram)
library(UpSetR)
library(ComplexUpset)


#Loading the samples from the first round of sequencing, received in January 2022

Cyclops_Control_1_5dpi <- Read10X("cyclops_control_filtered_feature_bc_matrix/")
Cyclops_Mloti_1_5dpi <- Read10X("cyclops_R7A_filtered_feature_bc_matrix/") 
WT_Control_1_dpi <- Read10X("Gifu_control_filtered_feature_bc_matrix/")
WT_Mloti_1_dpi <- Read10X("Gifu_R7A_filtered_feature_bc_matrix/") 


####Quality control and filtering

#Cyclops_Control_1_5dpi

#creating the seurat object
Cyclops_Control_1_5dpi_seurat <- CreateSeuratObject(counts = Cyclops_Control_1_5dpi_seurat, project = "Cyclops_Control_1_5dpi_seurat", min.cells = 3, min.features = 200)

#Adding metadata
Cyclops_Control_1_5dpi_seurat <- AddMetaData(object = Cyclops_Control_1_5dpi_seurat, col.name = "Sample", metadata = "Cyclops_Control_5dpi")
Cyclops_Control_1_5dpi_seurat <- AddMetaData(object = Cyclops_Control_1_5dpi_seurat, col.name = "Genotype", metadata = "Cyclops")
Cyclops_Control_1_5dpi_seurat <- AddMetaData(object = Cyclops_Control_1_5dpi_seurat, col.name = "Rhizobia", metadata = "Mock")

#Finding the percent of mt and chloroplast genes
Cyclops_Control_1_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(Cyclops_Control_1_5dpi_seurat, pattern = "LotjaGiM1v")
Cyclops_Control_1_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(Cyclops_Control_1_5dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(Cyclops_Control_1_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
dim(Cyclops_Control_1_5dpi_seurat)

plot1 <- FeatureScatter(Cyclops_Control_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cyclops_Control_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

#filtering the seurat object
Cyclops_Control_1_5dpi_seurat_filt <- subset(Cyclops_Control_1_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(Cyclops_Control_1_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
FeatureScatter(Cyclops_Control_1_5dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dim(Cyclops_Control_1_5dpi_seurat_filt)

#Saving the object as RDS
saveRDS(Cyclops_Control_1_5dpi_seurat_filt, "Cyclops_Control_1_5dpi_seurat_filt.rds")



####Cyclops_Mloti_1_5dpi
Cyclops_Mloti_1_5dpi_seurat <- CreateSeuratObject(counts = Cyclops_Mloti_1_5dpi_seurat, project = "Cyclops_Mloti_1_5dpi_seurat", min.cells = 3, min.features = 200)

Cyclops_Mloti_1_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_1_5dpi_seurat, col.name = "Sample", metadata = "Cyclops_R7A_5dpi")
Cyclops_Mloti_1_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_1_5dpi_seurat, col.name = "Genotype", metadata = "Cyclops")
Cyclops_Mloti_1_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_1_5dpi_seurat, col.name = "Rhizobia", metadata = "R7A")


Cyclops_Mloti_1_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(Cyclops_Mloti_1_5dpi_seurat, pattern = "LotjaGiM1v")
Cyclops_Mloti_1_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(Cyclops_Mloti_1_5dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(Cyclops_Mloti_1_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
dim(Cyclops_Mloti_1_5dpi_seurat)

plot1 <- FeatureScatter(Cyclops_Mloti_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cyclops_Mloti_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

Cyclops_Mloti_1_5dpi_seurat_filt <- subset(Cyclops_Mloti_1_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(Cyclops_Mloti_1_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

dim(Cyclops_Mloti_1_5dpi_seurat_filt)


saveRDS(Cyclops_Mloti_1_5dpi_seurat_filt, "Cyclops_Mloti_1_5dpi_seurat_filt.rds")


####WT_Control_1_dpi
WT_Control_1_5dpi_seurat <- CreateSeuratObject(counts = WT_Control_1_5dpi_seurat, project = "WT_Control_1_5dpi_seurat", min.cells = 3, min.features = 200)

WT_Control_1_5dpi_seurat <- AddMetaData(object = WT_Control_1_5dpi_seurat, col.name = "Sample", metadata = "Gifu_Control_5dpi")
WT_Control_1_5dpi_seurat <- AddMetaData(object = WT_Control_1_5dpi_seurat, col.name = "Genotype", metadata = "Gifu")
WT_Control_1_5dpi_seurat <- AddMetaData(object = WT_Control_1_5dpi_seurat, col.name = "Rhizobia", metadata = "Mock")


WT_Control_1_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Control_1_5dpi_seurat, pattern = "LotjaGiM1v")
WT_Control_1_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Control_1_5dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(WT_Control_1_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
dim(WT_Control_1_5dpi_seurat)

plot1 <- FeatureScatter(WT_Control_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Control_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

FeatureScatter(WT_Control_1_5dpi_seurat,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
WT_Control_1_5dpi_seurat_filt <- subset(WT_Control_1_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(WT_Control_1_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

saveRDS(WT_Control_1_5dpi_seurat_filt, "WT_Control_1_5dpi_seurat_filt.rds")


####WT_Mloti_1_5dpi
WT_Mloti_1_5dpi_seurat <- CreateSeuratObject(counts = WT_Mloti_1_5dpi_seurat, project = "WT_Mloti_1_5dpi_seurat", min.cells = 3, min.features = 200)

WT_Mloti_1_5dpi_seurat <- AddMetaData(object = WT_Mloti_1_5dpi_seurat, col.name = "Sample", metadata = "Gifu_R7A_5dpi")
WT_Mloti_1_5dpi_seurat <- AddMetaData(object = WT_Mloti_1_5dpi_seurat, col.name = "Genotype", metadata = "Gifu")
WT_Mloti_1_5dpi_seurat <- AddMetaData(object = WT_Mloti_1_5dpi_seurat, col.name = "Rhizobia", metadata = "R7A")


WT_Mloti_1_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Mloti_1_5dpi_seurat, pattern = "LotjaGiM1v")
WT_Mloti_1_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Mloti_1_5dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(WT_Mloti_1_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
dim(WT_Mloti_1_5dpi_seurat)

plot1 <- FeatureScatter(WT_Mloti_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Mloti_1_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

WT_Mloti_1_5dpi_seurat_filt <- subset(WT_Mloti_1_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(WT_Mloti_1_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

saveRDS(WT_Mloti_1_5dpi_seurat_filt, "WT_Mloti_1_5dpi_seurat_filt.rds")


#Loading the samples from the second run, received in June 2022############
Cyclops_Control_2_5dpi <- Read10X("Cyclops_Control_filtered_feature_bc_matrix/")
Cyclops_Mloti_2_5dpi <- Read10X("Cyclops_R7A_filtered_feature_bc_matrix/") 
WT_Control_2_5dpi <- Read10X("Gifu_Control_filtered_feature_bc_matrix/")
WT_Mloti_2_5dpi <- Read10X("Gifu_R7A_filtered_feature_bc_matrix/")


#Cyclops_Control_2
Cyclops_Control_2_5dpi_seurat <- CreateSeuratObject(counts = Cyclops_Control_2_5dpi_seurat, project = "Cyclops_Control_2_5dpi_seurat", min.cells = 3, min.features = 200)
Cyclops_Control_2_5dpi_seurat <- AddMetaData(object = Cyclops_Control_2_5dpi_seurat, col.name = "Sample", metadata = "Cyclops_Control_5dpi")
Cyclops_Control_2_5dpi_seurat <- AddMetaData(object = Cyclops_Control_2_5dpi_seurat, col.name = "Genotype", metadata = "Cyclops")
Cyclops_Control_2_5dpi_seurat <- AddMetaData(object = Cyclops_Control_2_5dpi_seurat, col.name = "Rhizobia", metadata = "Mock")


Cyclops_Control_2_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(Cyclops_Control_2_5dpi_seurat, pattern = "LotjaGiM1v")
Cyclops_Control_2_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(Cyclops_Control_2_5dpi_seurat, pattern = "LotjaGiC1v")

VlnPlot(Cyclops_Control_2_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(Cyclops_Control_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cyclops_Control_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

Cyclops_Control_2_5dpi_seurat_filt <- subset(Cyclops_Control_2_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(Cyclops_Control_2_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
FeatureScatter(Cyclops_Control_2_5dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

saveRDS(Cyclops_Control_2_5dpi_seurat_filt, "Cyclops_Control_2_5dpi_seurat_filt.rds")



#Cyclops_Mloti_2
Cyclops_Mloti_2_5dpi_seurat <- CreateSeuratObject(counts = Cyclops_Mloti_2_5dpi_seurat, project = "Cyclops_Mloti_2_5dpi_seurat", min.cells = 3, min.features = 200)
Cyclops_Mloti_2_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_2_5dpi_seurat, col.name = "Sample", metadata = "Cyclops_R7A_5dpi")
Cyclops_Mloti_2_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_2_5dpi_seurat, col.name = "Genotype", metadata = "Cyclops")
Cyclops_Mloti_2_5dpi_seurat <- AddMetaData(object = Cyclops_Mloti_2_5dpi_seurat, col.name = "Rhizobia", metadata = "R7A")


Cyclops_Mloti_2_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(Cyclops_Mloti_2_5dpi_seurat, pattern = "LotjaGiM1v")
Cyclops_Mloti_2_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(Cyclops_Mloti_2_5dpi_seurat, pattern = "LotjaGiC1v")

VlnPlot(Cyclops_Mloti_2_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(Cyclops_Mloti_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cyclops_Mloti_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

Cyclops_Mloti_2_5dpi_seurat_filt <- subset(Cyclops_Mloti_2_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(Cyclops_Mloti_2_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
FeatureScatter(Cyclops_Mloti_2_5dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

saveRDS(Cyclops_Mloti_2_5dpi_seurat_filt, "Cyclops_Mloti_2_5dpi_seurat_filt.rds")


#WT_Control_2_5dpi
WT_Control_2_5dpi_seurat <- CreateSeuratObject(counts = WT_Control_2_5dpi_seurat, project = "WT_Control_2_5dpi_seurat", min.cells = 3, min.features = 200)
WT_Control_2_5dpi_seurat <- AddMetaData(object = WT_Control_2_5dpi_seurat, col.name = "Sample", metadata = "Gifu_Control_5dpi")
WT_Control_2_5dpi_seurat <- AddMetaData(object = WT_Control_2_5dpi_seurat, col.name = "Genotype", metadata = "Gifu")
WT_Control_2_5dpi_seurat <- AddMetaData(object = WT_Control_2_5dpi_seurat, col.name = "Rhizobia", metadata = "Mock")


WT_Control_2_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Control_2_5dpi_seurat, pattern = "LotjaGiM1v")
WT_Control_2_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Control_2_5dpi_seurat, pattern = "LotjaGiC1v")

VlnPlot(WT_Control_2_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(WT_Control_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Control_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

WT_Control_2_5dpi_seurat_filt <- subset(WT_Control_2_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(WT_Control_2_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
FeatureScatter(WT_Control_2_5dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

saveRDS(WT_Control_2_5dpi_seurat_filt, "WT_Control_2_5dpi_seurat_filt.rds")


#WT_Mloti_2_5dpi
WT_Mloti_2_5dpi_seurat <- CreateSeuratObject(counts = WT_Mloti_2_5dpi_seurat, project = "WT_Mloti_2_5dpi_seurat", min.cells = 3, min.features = 200)
WT_Mloti_2_5dpi_seurat <- AddMetaData(object = WT_Mloti_2_5dpi_seurat, col.name = "Sample", metadata = "Gifu_R7A_5dpi")
WT_Mloti_2_5dpi_seurat <- AddMetaData(object = WT_Mloti_2_5dpi_seurat, col.name = "Genotype", metadata = "Gifu")
WT_Mloti_2_5dpi_seurat <- AddMetaData(object = WT_Mloti_2_5dpi_seurat, col.name = "Rhizobia", metadata = "R7A")


WT_Mloti_2_5dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Mloti_2_5dpi_seurat, pattern = "LotjaGiM1v")
WT_Mloti_2_5dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Mloti_2_5dpi_seurat, pattern = "LotjaGiC1v")

VlnPlot(WT_Mloti_2_5dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(WT_Mloti_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Mloti_2_5dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

WT_Mloti_2_5dpi_seurat_filt <- subset(WT_Mloti_2_5dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 10 & percent.chloroplast < 5)
VlnPlot(WT_Mloti_2_5dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)
FeatureScatter(WT_Mloti_2_5dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

saveRDS(WT_Mloti_2_5dpi_seurat_filt, "WT_Mloti_2_5dpi_seurat_filt.rds")


#Data integration#####


seurat.list <- list(Cyclops_Control_1_5dpi_seurat_filt, Cyclops_Mloti_1_5dpi_seurat_filt, WT_Control_1_5dpi_seurat_filt, WT_Mloti_1_5dpi_seurat_filt, Cyclops_Control_2_5dpi_seurat_filt, Cyclops_Mloti_2_5dpi_seurat_filt,  WT_Control_2_5dpi_seurat_filt,  WT_Mloti_2_5dpi_seurat_filt)


seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("percent.mt",  "percent.chloroplast"), method = "glmGamPoi", verbose = TRUE)
})

seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)


seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = seurat.features)


#Integrating using SCtransform

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
                                       anchor.features = seurat.features)

seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")



DefaultAssay(seurat.integrated) <- "integrated"

seurat.integrated <- RunPCA(object = seurat.integrated, verbose = FALSE)


#Clustering unsing Louvain algorithm 

seurat.integrated <- FindNeighbors(object = seurat.integrated, dims = 1:50)

seurat.integrated <- FindClusters(object = seurat.integrated)


#UMAP
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

#########
DefaultAssay(seurat.integrated) <- "RNA"
seurat.integrated <- NormalizeData(seurat.integrated)


saveRDS(seurat.integrated, "seurat_integrated.rds")


#####Find conserved markers for all clusters######
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(Gifu_Cyclops.integrated,
                       ident.1 = cluster,
                       grouping.var = "Sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}


# Iterate function across desired clusters

DefaultAssay(Gifu_Cyclops.integrated) <- "RNA"

library(multtest)
library(metap)

conserved_markers <- map_dfr(c(0:29), get_conserved)



##Cluster Annotation######

Gifu_Cyclops.integrated <-  readRDS("seurat_integrated.rds")

Idents(Gifu_Cyclops.integrated) <-  "seurat_clusters"

Gifu_Cyclops.integrated<- RenameIdents(object = Gifu_Cyclops.integrated,
                                       "4"="Cortex",
                                       "24"="Atrichoblasts",
                                       "21"="Atrichoblasts",
                                       "1"="Atrichoblasts",
                                       "12"="Cortex",
                                       "8"="Cortex",
                                       "6"="Cortex",
                                       "17"="Cortex",
                                       "13"="Cortex",
                                       "16"="Cortex",
                                       "25"="Endodermis",
                                       "20"="Endodermis",
                                       "28"="Endodermis",
                                       "2"= "Root tip",
                                       "10"= "Root tip",
                                       "19"= "Phloem",
                                       "27"="Phloem",
                                       "29"="Stele",
                                       "18"= "Trichoblasts",
                                       "22"= "Trichoblasts",
                                       "11"="Pericycle",
                                       "0"="Stele",
                                       "5"="Stele",
                                       "3"="Stele",
                                       "9"= "Stele",
                                       "26"="Cortex",
                                       "7"= "Stele",
                                       "15"="Stele",
                                       "14"="Stele",
                                       "23" = "Xylem")


Gifu_Cyclops.integrated@meta.data$Cell_type <- Gifu_Cyclops.integrated@active.ident





######################

#Seurat differential gene expression for WT

table(Gifu_Cyclops.integrated$Sample)

Gifu_Cyclops.integrated$Cluster_Sample <- paste0(Gifu_Cyclops.integrated$seurat_clusters, "_", Gifu_Cyclops.integrated$Sample)
Idents(Gifu_Cyclops.integrated) <- "Cluster_Sample"

clusters <- levels(Gifu_Cyclops.integrated$seurat_clusters)

for (i in clusters){
  tryCatch({
    print(i)
    test <- FindMarkers(Gifu_Cyclops.integrated, ident.1 = paste0(i, "_Gifu_R7A"),
                        ident.2 = paste0(i, "_Gifu_Control"), verbose = TRUE, min.pct = 0.01, test.use = "MAST") %>% filter(p_val_adj <= 0.05)
    test$gene <- rownames(test)
    print(nrow(test))
    
    print(nrow(test))
    
    names(test)[names(test) == "pct.1"] <- "Percent_Inoculated"
    names(test)[names(test) == "pct.2"] <- "Percent_Control"
    writexl::write_xlsx(test, paste0("DE_genes/WT_5dpi_min_0.01/", "WT_5dpi_DE_genes_MAST_Cluster_", i, ".xlsx"))
    
    if (nrow(test)==0|nrow(test)==0) stop("No DE genes found")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}




#List of markers for the cells expression NF-YA1 in Cortex and Root hairs

NFYA1_cells <- WhichCells(Gifu_Cyclops.integrated, expression = LotjaGi5g1v0106700 >0, ident = c("16_Gifu_R7A", "13_Gifu_R7A", "26_Gifu_R7A", "18_Gifu_R7A"))

NFYA1_cells_MAST <- FindMarkers(Gifu_Cyclops.integrated, ident.1 = NFYA1_cells,  verbose = TRUE, min.pct = 0.01, test.use = "MAST") %>% filter(p_val_adj <= 0.05)
NFYA1_cells_MAST$gene <- rownames(test_MAST)

NFYA1_cells_MAST_filtered<- filter(NFYA1_cells_MAST, pct.2 <= 0.02, avg_log2FC >0)



#List of markers for the cells expression NF-YA1 in Cortex

NFYA1_cells_C <- WhichCells(Gifu_Cyclops.integrated, expression = LotjaGi5g1v0106700 >0, ident = c("16_Gifu_R7A", "13_Gifu_R7A", "26_Gifu_R7A"))
DimPlot(Gifu_Cyclops.integrated, label=F, combine= T, pt.size = 2, sizes.highlight=2, cells.highlight= list(NFYA1_cells_C), cols.highlight = c("#7f0023"), cols= "lightgrey", split.by = "Sample")&theme_minimal()&NoGrid()&NoLegend()&NoAxes()&theme(plot.title = element_blank())

NFYA1_cells_C_MAST <- FindMarkers(Gifu_Cyclops.integrated, ident.1 = NFYA1_cells_C, ident.2 = c("13_Cyclops_R7A", "16_Cyclops_R7A", "26_Cyclops_R7A"),  verbose = TRUE, min.pct = 0.01, test.use = "MAST") %>% filter(p_val_adj <= 0.05)
NFYA1_cells_C_MAST$gene <- rownames(NFYA1_cells_C_MAST)


NFYA1_cells_C_MAST_filtered<- filter(NFYA1_cells_C_MAST, pct.2 <= 0.02)



####Subset Root hairs clusters#####

Cluster_18_22 <- subset(Gifu_Cyclops.integrated,  idents = c("18", "22"))
Cluster_18_22$old_clusters <-Cluster_18_22$seurat_clusters

DefaultAssay(Cluster_18_22) <-  "integrated"

Cluster_18_22 <- RunPCA(Cluster_18_22)
ElbowPlot(Cluster_18_22)
Cluster_18_22 <- FindNeighbors(Cluster_18_22, dims = 1:30)
Cluster_18_22 <- RunUMAP(Cluster_18_22, dims = 1:30)
Cluster_18_22 <- RunTSNE(Cluster_18_22, dims = 1:30)

Cluster_18_22 <- FindClusters(Cluster_18_22)

DefaultAssay(Cluster_18_22) <-  "RNA"

DimPlot(object = Cluster_18_22, reduction = "umap",  label = TRUE, repel = TRUE, pt.size=2,  label.size = 7)  + NoAxes()
DimPlot(object = Cluster_18_22, reduction = "umap",  label = TRUE, repel = TRUE, pt.size=2,  split.by="Genotype",  label.size = 7)  + NoAxes()
DimPlot(object = Cluster_18_22, reduction = "umap",  label = TRUE, repel = TRUE, pt.size=2,  group.by="old_clusters", label.size = 7)  + NoAxes()


###Find markers cyclops subcluster RH 5###

Cluster_18_22$Cluster_Sample <- paste0(Cluster_18_22$seurat_clusters, "_", Cluster_18_22$Sample)
Idents(Cluster_18_22) <- "Cluster_Sample"

cyclops_RH <- subset(Cluster_18_22, ident = "Cyclops")
Idents(cyclops_RH) <- "Cluster_Sample"

cyclops_RH_5 <- FindMarkers(Gifu_Cyclops.integrated, ident.1 = WhichCells(cyclops_RH, ident="5_Cyclops_R7A"),  verbose = TRUE, min.pct = 0.01, test.use = "MAST", only.pos = T) %>% filter(p_val_adj <= 0.05)
cyclops_RH_5$gene <- rownames(cyclops_RH_5)


####Find markers Gifu subcluster RH 6 #####

Gifu_RH <- subset(Cluster_18_22, ident = "Gifu")
Idents(Gifu_RH) <- "Cluster_Sample"

Gifu_RH_6 <- FindMarkers(Gifu_Cyclops.integrated, ident.1 =WhichCells(Gifu_RH, ident="6_Gifu_R7A"),  verbose = TRUE, min.pct = 0.01, test.use = "MAST", only.pos = T) %>% filter(p_val_adj <= 0.05)
Gifu_RH_6$gene <- rownames(Gifu_RH_6)

Gifu_RH_6_filter <- filter(Gifu_RH_6, pct.2<= 0.02)




