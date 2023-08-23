
library(Seurat)
library(Scissor)
library(preprocessCore)
library(stringr)


#Script for the 10 dpi Data 

bulk_dataset <- as.matrix(read.delim("../Bulk_RNAseq/Genes_Simon/RH_rhizobia_bulk_data.txt"))

R7A_inoculation <- c(0,0,0, 0,0,0, 1, 1, 1, 1, 1, 1)

phenotype <- R7A_inoculation

tag <- c('Control', 'R7A_inoculated')
table(phenotype)


RH_Clusters <- readRDS("RH_Clusters/Cluster2_9_integrated_assay.RDS")

DimPlot(object = RH_Clusters, reduction = "umap",  label = T, repel = TRUE, split.by = "Condition", pt.size = 2)


RH_Clusters@graphs$RNA_snn <- RH_Clusters@graphs$integrated_snn


infos_RH_Clusters <- Scissor(bulk_dataset, RH_Clusters, phenotype,  tag = tag, alpha = NULL, 
                             family = "binomial")

#alpha was selected to 0.8

Scissor_select <- rep(0, ncol(RH_Clusters))
names(Scissor_select) <- colnames(RH_Clusters)


Scissor_select[infos_RH_Clusters$Scissor_pos] <- "Scissor+"
Scissor_select[infos_RH_Clusters$Scissor_neg] <- "Scissor-"


RH_Clusters <- AddMetaData(RH_Clusters, metadata = Scissor_select, col.name = "Scissor_RH_Bulk")


DimPlot(RH_Clusters, reduction = 'umap', group.by = 'Scissor_RH_Bulk', cols = c('grey','indianred1','royalblue'),
        pt.size = 1.2, order = c("Scissor-", "Scissor+"), split.by = "Condition")

saveRDS(RH_Clusters, "Scissor_for_article/RH_Clusters_10days_Scissor.RDS")


evaluate_summary <- evaluate.cell('Scissor_for_article/RH_Clusters_10_days_Bulk.RData', infos_RH_Clusters, FDR = 0.05, bootstrap_n = 100)
evaluate_summary[1:5,1:4]


#Marker genes for the Scissor positive cells - 10 day dataset

Annotations <- read.delim("../../../Reference files/LjGifu_1.2_functional_annotation_for_Single_Cell.tsv", sep = "\t")
Nodulation_genes <- readxl::read_xlsx("../../../Nodulation_curated_genes.xlsx")
colnames(Nodulation_genes)
names(Nodulation_genes)[names(Nodulation_genes) == "ID"] <- "gene"



RH_Clusters <- readRDS("RH_Clusters_10days_Scissor.RDS")

levels(as.factor(RH_Clusters$Scissor_RH_Bulk))


Idents(RH_Clusters) <- "Scissor_RH_Bulk"

markers_RH_positive_cells1 <- FindMarkers(RH_Clusters, ident.1 = "Scissor+",  verbose = TRUE, min.pct = 0.01, test.use = "MAST", only.pos = T) %>% filter(p_val_adj <= 0.05)
markers_RH_positive_cells1$gene <- rownames(markers_RH_positive_cells1)
markers_RH_positive_cells1 <- left_join(markers_RH_positive_cells1, Annotations, by="gene")
markers_RH_positive_cells1 <- left_join(markers_RH_positive_cells1, Nodulation_genes, by="gene")

markers_RH_positive_cells1_filtered <- filter(markers_RH_positive_cells1, pct.2 <= 0.02)


{feature <- "LotjaGi5g1v0269800-LC"
FeaturePlot(RH_Clusters, 
            reduction = "umap", 
            features = feature, 
            order = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            label = F,
            pt.size = 2,
            split.by = "Condition") &NoAxes() & theme(legend.position = "right")& scale_colour_gradientn(colours  = c("grey",brewer.pal(n = 9, name = "YlGnBu")))
}


writexl::write_xlsx(markers_RH_positive_cells1, "Markers_RH_Scissor_positive_Cells.xlsx")

