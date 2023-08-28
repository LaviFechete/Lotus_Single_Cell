
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
library(multtest)
library(plotly)
library(scater)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(VennDiagram)
library(UpSetR)
library(ComplexUpset)


#Data pre-processing and integration 10dpi####


WT_Control1_10dpi <- Read10X("/WT_Control1_10dpi_filtered_feature_bc_matrix/")
WT_Control2_10dpi <- Read10X("/WT_Control2_10dpi_filtered_feature_bc_matrix/") 
WT_Mloti1_10dpi <- Read10X("/WT_Mloti1_10dpi_filtered_feature_bc_matrix/")
WT_Mloti2_10dpi <- Read10X("/WT_Mloti2_10dpi_filtered_feature_bc_matrix/")



#WT_Control1_10dpi
WT_Control1_10dpi_seurat <- CreateSeuratObject(counts = WT_Control1_10dpi, project = "WT_Control1_10dpi_seurat", min.cells = 3, min.features = 200)
WT_Control1_10dpi_seurat <- AddMetaData(object = WT_Control1_10dpi_seurat, metadata = "Control", col.name = "Condition")

WT_Control1_10dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Control1_10dpi_seurat, pattern = "LotjaGiM1v")
WT_Control1_10dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Control1_10dpi_seurat, pattern = "LotjaGiC1v")


VlnPlot(WT_Control1_10dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)


plot1 <- FeatureScatter(WT_Control1_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Control1_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2


FeatureScatter(WT_Control1_10dpi_seurat,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

WT_Control1_10dpi_seurat_filt <- subset(WT_Control1_10dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 5 & percent.chloroplast < 5)
VlnPlot(WT_Control1_10dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

FeatureScatter(WT_Control1_10dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


#WT_Control2_10dpi
WT_Control2_10dpi_seurat <- CreateSeuratObject(counts = WT_Control2_10dpi, project = "WT_Control2_10dpi_seurat", min.cells = 3, min.features = 200)
WT_Control2_10dpi_seurat <- AddMetaData(object = WT_Control2_10dpi_seurat, metadata = "Control", col.name = "Condition")


WT_Control2_10dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Control2_10dpi_seurat, pattern = "LotjaGiM1v")
WT_Control2_10dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Control2_10dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(WT_Control1_10dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(WT_Control2_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Control2_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

FeatureScatter(WT_Control2_10dpi_seurat,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

WT_Control2_10dpi_seurat_filt <- subset(WT_Control2_10dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 5 & percent.chloroplast < 5)
VlnPlot(WT_Control2_10dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

FeatureScatter(WT_Control2_10dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


#WT_Mloti1_10dpi

WT_Mloti1_10dpi_seurat <- CreateSeuratObject(counts = WT_Mloti1_10dpi, project = "WT_Mloti1_10dpi_seurat", min.cells = 3, min.features = 200)

WT_Mloti1_10dpi_seurat <- AddMetaData(object = WT_Mloti1_10dpi_seurat, metadata = "R7A", col.name = "Condition")


WT_Mloti1_10dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Mloti1_10dpi_seurat, pattern = "LotjaGiM1v")
WT_Mloti1_10dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Mloti1_10dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(WT_Mloti1_10dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(WT_Mloti1_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Mloti1_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

FeatureScatter(WT_Mloti1_10dpi_seurat,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

WT_Mloti1_10dpi_seurat_filt <- subset(WT_Mloti1_10dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 5 & percent.chloroplast < 5)
VlnPlot(WT_Mloti1_10dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

FeatureScatter(WT_Mloti1_10dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



##WT_Mloti2_10dpi

WT_Mloti2_10dpi_seurat <- CreateSeuratObject(counts = WT_Mloti2_10dpi, project = "WT_Mloti2_10dpi_seurat", min.cells = 3, min.features = 200)

WT_Mloti2_10dpi_seurat <- AddMetaData(object = WT_Mloti2_10dpi_seurat, metadata = "R7A", col.name = "Condition")

WT_Mloti2_10dpi_seurat[["percent.mt"]] <- PercentageFeatureSet(WT_Mloti2_10dpi_seurat, pattern = "LotjaGiM1v")
WT_Mloti2_10dpi_seurat[["percent.chloroplast"]] <- PercentageFeatureSet(WT_Mloti2_10dpi_seurat, pattern = "LotjaGiC1v")
VlnPlot(WT_Mloti2_10dpi_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

plot1 <- FeatureScatter(WT_Mloti2_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_Mloti2_10dpi_seurat, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2

FeatureScatter(WT_Mloti2_10dpi_seurat,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

WT_Mloti2_10dpi_seurat_filt <- subset(WT_Mloti2_10dpi_seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200 &  nFeature_RNA < 7500 & percent.mt < 5 & percent.chloroplast < 5)
VlnPlot(WT_Mloti2_10dpi_seurat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 4)

FeatureScatter(WT_Mloti2_10dpi_seurat_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



#Data integration####

options(future.globals.maxSize = 8000 * 1024^2)

Gifu.list <- list(WT_Control1_10dpi_seurat, WT_Control2_10dpi_seurat, WT_Mloti1_10dpi_seurat, WT_Mloti2_10dpi_seurat)

Gifu.list <- lapply(X = Gifu.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("percent.mt",  "percent.chloroplast"), verbose = TRUE)
})

Gifu.features <- SelectIntegrationFeatures(object.list = Gifu.list, nfeatures = 3000)


Gifu.list <- PrepSCTIntegration(object.list = Gifu.list, anchor.features = Gifu.features)


#Integrating using SCtransform and reference

Gifu.anchors <- FindIntegrationAnchors(object.list = Gifu.list, normalization.method = "SCT", 
                                       anchor.features = Gifu.features, reference = c(1,2))

Gifu.integrated <- IntegrateData(anchorset = Gifu.anchors, normalization.method = "SCT")


saveRDS(Gifu.integrated, "Gifu.integrated2_SCT.rds")



DefaultAssay(Gifu.integrated) <- "integrated"

Gifu.integrated <- RunPCA(object = Gifu.integrated, verbose = FALSE)


Gifu.integrated <- FindNeighbors(object = Gifu.integrated, dims = 1:50)

#Clustering unsing Louvain algorithm
Gifu.integrated <- FindClusters(object = Gifu.integrated)

#UMAP
Gifu.integrated <- RunUMAP(object = Gifu.integrated, dims = 1:50)

Gifu.integrated <- RunTSNE(object = Gifu.integrated,  check_duplicates = FALSE)


saveRDS(Gifu.integrated, "Gifu.integrated2_SCT_UMAP.rds")



DefaultAssay(Gifu.integrated) <- "RNA"

Gifu.integrated <- NormalizeData(Gifu.integrated, verbose = FALSE)


#######Find conserved markers for all clusters----

get_conserved <- function(cluster){
  FindConservedMarkers(Gifu.integrated,
                       ident.1 = cluster,
                       grouping.var = "Condition",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}



# Iterate function across desired clusters


DefaultAssay(Gifu.integrated) <- "RNA"

conserved_markers <- map_dfr(c(0:31), get_conserved)

conserved_markers$Average_log2FC <- (conserved_markers$Control_avg_log2FC + conserved_markers$Inoculated_avg_log2FC) /2


#######Cluster Annotation#######

Idents(Gifu.integrated) <- "seurat_clusters"

Gifu.integrated<- RenameIdents(object = Gifu.integrated,
                               "0"= "Root tip",
                               "1"="Cortex",
                               "2"= "Trichoblasts",
                               "3"="Stele",
                               "4"="Atrichoblasts",
                               "5"="Cortex",
                               "6"="Pericycle",
                               "7"= "Root tip",
                               "8"="Cortex",
                               "9"= "Trichoblasts",
                               "10"="Cortex",
                               "11"="Atrichoblasts",
                               "12"="Cortex", 
                               "13"= "Root tip",
                               "14"="Cortex",
                               "15"="Cortex",
                               "16"="Cortex",
                               "17"="Stele",
                               "18"="Cortex",
                               "19"= "Root tip",
                               "20"="Endodermis",
                               "21"="Meristem",
                               "22"="Cortex",
                               "23"="Atrichoblasts",
                               "24"="Endodermis",
                               "25"= "Phloem",
                               "26"="Stele",
                               "27"="Endodermis",
                               "28"="Quiescent center",
                               "29" = "Xylem",
                               "30"= "Phloem",
                               "31"="Cortex")

Gifu.integrated$Cell_type <- Idents(Gifu.integrated)

Idents(Gifu.integrated) <- "Cell_type"
DimPlot(object = Gifu.integrated, reduction = "umap",  label = F, repel = TRUE, pt.size = 1, cols = my_cols) & NoAxes()


#####Differential gene expression######################################

Gifu.integrated$Cluster_Condition <- paste0(Gifu.integrated$seurat_clusters, "_", Gifu.integrated$Condition)
Idents(Gifu.integrated) <- "Cluster_Condition"
table(Gifu.integrated$Cluster_Condition)

clusters <- levels(Gifu.integrated$seurat_clusters)

for (i in clusters){
  tryCatch({
    print(i)
    test <- FindMarkers(Gifu.integrated, ident.1 = paste0(i, "_Inoculated"),
                        ident.2 = paste0(i, "_Control"), verbose = TRUE, min.pct = 0.01, test.use = "MAST") %>% filter(p_val_adj <= 0.05)
    test$gene <- rownames(test)
    print(nrow(test))
    
    for (x in 1:nrow(test)){
      test$Percent_Total[x] <- sum(GetAssayData(object = Gifu.integrated, slot = "data")[test$gene[x],]>0)/nrow(Gifu.integrated@meta.data)
    }
    
    test1 <-  filter(test, Percent_Total <= 0.3)
    test1 <- left_join(test1, Annotations, by="gene")
    test1 <- left_join(test1, Nodulation_genes, by="gene")
    print(nrow(test1))
    
    names(test1)[names(test) == "pct.1"] <- "Percent_R7A"
    names(test1)[names(test) == "pct.2"] <- "Percent_Control"
    writexl::write_xlsx(test1, paste0("DE_genes/Seurat_MAST_Mock_vs_Inoculated_10dpi/", "DE_genes_MAST_Cluster", i, ".xlsx"))
    
    if (nrow(test)==0|nrow(test1)==0) stop("No DE genes found")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


###Subsetting the RH Clusters####

RH_clusters <- subset(Gifu.integrated,  idents = c("9", "2"))
RH_clusters$old_clusters <- RH_clusters$seurat_clusters 
DefaultAssay(RH_clusters) <-  "integrated"

RH_clusters <- RunPCA(RH_clusters)
ElbowPlot(RH_clusters)
RH_clusters <- FindNeighbors(RH_clusters, dims = 1:30)
RH_clusters <- RunUMAP(RH_clusters, dims = 1:30)
RH_clusters <- FindClusters(RH_clusters)

Idents(RH_clusters) <- "seurat_clusters"
DimPlot(object = RH_clusters, reduction = "umap",  label = TRUE, repel = TRUE, pt.size=2,  group.by="orig.ident",  label.size = 7)  + NoAxes()

DimPlot(object = RH_clusters, reduction = "umap",  label = TRUE, repel = TRUE, split.by="Condition",  pt.size=2, label.size = 7)  + NoAxes()



####Identifying the infected cells based on the cells in cluster 19 from the Inoculated samples alone#####

#Import the cell IDs for cluster 19

Cluster19_Cortex <- read.delim("Cells_Cluster19_Cortex.txt", sep = " ")
Cluster19_RH <- read.delim("Cells_Cluster19_RH.txt", sep = " ")

Cluster19_ids <- c(Cluster19_Cortex$x, Cluster19_RH$x)


Cluster19_markers <- FindMarkers(Gifu.integrated, ident.1 = Cluster19_ids, verbose = TRUE, test.use = "MAST", min.pct = 0.01) %>% filter(p_val_adj <= 0.05)

Cluster19_markers$gene <- rownames(Cluster19_markers)


#Subset RH and cortex clusters#########

RH_Cortex_10dpi <- subset(Gifu.integrated, idents = c("1", "2", "5","8", "9", "10", "12", "14", "15", "16", "18","22"))
DefaultAssay(RH_Cortex_10dpi) <-  "integrated"

RH_Cortex_10dpi <- RunPCA(RH_Cortex_10dpi)
RH_Cortex_10dpi <- FindNeighbors(RH_Cortex_10dpi, dims = 1:50)
RH_Cortex_10dpi <- RunUMAP(RH_Cortex_10dpi, dims = 1:50)


#Markers for the RH infected cells

Cluster19_markers_RH <- FindMarkers(RH_cortex_10dpi, ident.1 = Cluster19_RH$x, verbose = TRUE, test.use = "MAST", min.pct = 0.01) %>% filter(p_val_adj <= 0.05)
Cluster19_markers_RH$gene <- rownames(Cluster19_markers_RH)

Cluster19_markers_RH_filtered <- filter(Cluster19_markers_RH, pct.2 <= 0.023)

#Markers for the cortical infected cells
Cluster19_markers_cortex <- FindMarkers(RH_cortex_10dpi, ident.1 = Cluster19_Cortex$x, verbose = TRUE, test.use = "MAST", min.pct = 0.01) %>% filter(p_val_adj <= 0.05)

Cluster19_markers_cortex$gene <- rownames(Cluster19_markers_cortex)

Cluster19_markers_cortex_filtered <- filter(Cluster19_markers_cortex, pct.2 <= 0.023)


RH_only <- anti_join(Cluster19_markers_RH_filtered, Cluster19_markers_cortex_filtered, by="gene")

Cortex_only <- anti_join(Cluster19_markers_cortex_filtered, Cluster19_markers_RH_filtered, by="gene")

RH_and_Cortex <- inner_join(Cluster19_markers_RH_filtered, Cluster19_markers_cortex_filtered, by="gene")



####Identifying markers for the bacteroid cells########

Idents(Gifu.integrated) <- "seurat_clusters"

my_cells <- WhichCells(Gifu.integrated, idents= c("8"))

LB1 <- "LotjaGi5g1v0024700"; LB2<- "LotjaGi5g1v0024900"; LB3 <- "LotjaGi5g1v0046500"

LB1_ids <- which(GetAssayData(object = Gifu.integrated, slot = "data")[LB1,my_cells]>0)
LB2_ids <- which(GetAssayData(object = Gifu.integrated, slot = "data")[LB2,my_cells]>0)
LB3_ids <-which(GetAssayData(object = Gifu.integrated, slot = "data")[LB3,my_cells]>0)

Bacterioids_Cells <- unique(c(names(LB1_ids), names(LB2_ids), names(LB3_ids)))

DimPlot(Gifu.integrated, label=T, cells.highlight= list(Bacterioids_Cells), cols.highlight = c("darkblue"), cols= "grey", split.by = "Condition")

Bacterioids_Cells_markers <- FindMarkers(Gifu.integrated, ident.1 = Bacterioids_Cells, verbose = TRUE, test.use = "MAST", min.pct = 0.01) %>% filter(p_val_adj <= 0.05)
Bacterioids_Cells_markers$gene <- rownames(Bacterioids_Cells_markers)

Bacterioids_Cells_markers_filtered <- filter(Bacterioids_Cells_markers, pct.2 <= 0.02)


####Identifying markers for the nodule cells########

CA1 <- "LotjaGi1g1v0180800"

Idents(Gifu.integrated) <- "seurat_clusters"

my_cells <- WhichCells(Gifu.integrated, idents= c("14"))

NOD_ids<- which(GetAssayData(object = Gifu.integrated, slot = "data")[CA1,my_cells]>0)

NOD_ids_markers <- FindMarkers(Gifu.integrated, ident.1 = names(NOD_ids), verbose = TRUE, test.use = "MAST", min.pct = 0.01) %>% filter(p_val_adj <= 0.05)
NOD_ids_markers$gene <- rownames(NOD_ids_markers)

NOD_ids_markers_filtered <- filter(NOD_ids_markers, pct.2 <= 0.02)



