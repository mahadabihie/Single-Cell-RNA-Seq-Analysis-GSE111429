#GSE111429
  # Identification of a Zeb1 expressing basal stem cell subpopulation in the prostate:https://www.nature.com/articles/s41467-020-14296-y#data-availability
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111429
  # https://github.com/HelenHeZhu/StemCell
  # https://www.youtube.com/watch?v=IjJOTJsd4Mg
  # https://github.com/Lindseynicer/Analysing-10X-ScRNA-data-Wang-et-al-2020-Nat-Commun-

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1. remove low quality cells and contaminated non-epithelial cells
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#library
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(data.table)
library(ggplot2)
library(pathview)
library(biomaRt)
library(fgsea)
library(msigdbr)
library(data.table)
library(ggplot2)

#read in data
data <- Read10X(data.dir = "GSE111429_RAW")

PG_all    <- CreateSeuratObject(count = data, min.cells = 0, min.genes = 0, project = "PG")

#PercentageFeatureSet(): calculate the percentage of all the counts belonging to a subset of the possible features for each cell. This is useful when trying to compute the percentage of transcripts that map to mitochondrial genes
PG_all <- PercentageFeatureSet(PG_all, pattern = "^mt-", col.name = "percent.mt")
before <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))

#Poor-quality cells with less than 1000 genes detected, less than 5000 UMIs or more than 5% UMI mapped to mitochondria genes were removed. 
PG_all <- subset(PG_all, subset = nFeature_RNA > 1400 & nCount_RNA > 5000 & percent.mt < 5)
after <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))

#Supplementary Figure 5e: Violin plots displaying the distribution of detected numbers of gene, uniqumolecular identifiers (UMI) and percentage of UMIs mapped to mitochondrial genes before (e) and after (f) removal of low-quality cells (cells with less than 5000 UMIs or less than 1400 genes detected or more than 5% UMI mapped to mitochondria genes were removed).
#compare the difference before and after filtering
ggarrange(before,after, ncol = 2, nrow = 1)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Step 1: Run the Standard Workflow for Visualization and Clustering
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PG_all <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all <- FindVariableFeatures(PG_all,
                               nfeatures = 1606,
                               mean.cutoff = c(0.0125,3),
                               dispersion.cutoff = c(0.5, Inf))
PG_all <- ScaleData(PG_all, vars.to.regress=c('nCount_RNA','percent.mt'))
PG_all <- RunPCA(PG_all, features = PG_all@assays$RNA@var.features, verbose = TRUE, 
                 ndims.print=1:5, nfeatures.print = 5)
PG_all <- FindNeighbors(PG_all, dims = 1:12)
PG_all <- FindClusters(PG_all, resolution=c(0.2,0.3,0.4,0.5))

#runTSNE used for non-linear reduction of the dimension
PG_all <- RunTSNE(PG_all, dims = 1:11, do.fast=TRUE)

#**Idents():** finds clusters
Idents(object = PG_all) <- PG_all$RNA_snn_res.0.2

#corresponds to figure 6a, non-linear reduction doesn't produce the exact same each ti  me
DimPlot(PG_all, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6)
PCAPlot(PG_all)

#### retrieve non-epithelial cells
  #Cd74, Cd72 and Cd54 are immune cell markers
FeaturePlot(PG_all, features =c('Cd74','Cd72','Cd54') , cols =c('yellow','red'), 
            reduction ='tsne', pt.size = 0.7,ncol = 3, label = T)
  #'Eng','S1pr1', and 'Emcn'are endothelial cell markers (cluster 4)
FeaturePlot(PG_all, features = c('Eng','S1pr1','Emcn'), cols =c('yellow','red'), 
            reduction = 'tsne',pt.size = 0.7,ncol=3, label = T)

#what is this for?
#Immu.Cluster <- FindMarkers(PG_all, ident.1='9', only.pos=TRUE, min.pct=0.25)
#Endo.Cluster <- FindMarkers(PG_all, ident.1='4', only.pos=TRUE, min.pct=0.25)
#Immu.Cluster$Symbol <- rownames(Immu.Cluster)
#Endo.Cluster$Symbol <- rownames(Endo.Cluster)

# violin plots for each class of markers
v1 <- VlnPlot(PG_all, features = c("EpCAM","Cdh1","Krt5", "Krt14"), ncol = 4) #Epithelial cellmarkers
v2 <- VlnPlot(PG_all, features = c("Igfbp4","Fn1","Gng11"), ncol = 3) # stromal markers
v3 <- VlnPlot(PG_all, features = c("Eng","S1pr1","Emcn"), ncol = 3) # endothelial markers
v4 <- VlnPlot(PG_all, features = c("Cd74","Cd72"), ncol = 2) #immune markers

#this is numbered 0-9, but matches with 1-10
ggarrange(v1, v2, v3, v4, ncol = 1)

# standardized heatmap
  # seurat function
  ### heatmap for non-epithelial cells markers
geneSets <- c('Cd74','Cd72','Cd54','Eng','S1pr1','Emcn')
DoHeatmap(PG_all, features = geneSets)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Step 2: Clustering Using Seurat
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #re-do standard workflow for a subset of the data 
PG_all <- subset(PG_all, idents = c("4","8"), invert = TRUE)
DefaultAssay(PG_all) <- "RNA"
PG_all <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all <- FindVariableFeatures(PG_all, mean.function = ExpMean, 
                               dispersion.function= LogVMR,
                               mean.cutoff = c(0.0125, 3),
                               dispersion.cutoff = c(0.5, Inf))
length(PG_all@assays$RNA@var.features)
#
PG_all <- ScaleData(PG_all, vars.to.regress=c('nCount_RNA','percent.mt'))
PG_all <- RunPCA(PG_all, features = PG_all@assays$RNA@var.features, verbose = TRUE, 
                 ndims.print=1:5, nfeatures.print = 5)
PG_all <- FindNeighbors(PG_all, dims = 1:12)
PG_all <- FindClusters(PG_all, resolution=c(0.2,0.3,0.4,0.5))
PG_all <- RunTSNE(PG_all, dims = 1:12, do.fast=TRUE)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#save the work from this point on (standard work flow + subset reclustering)
saveRDS(PG_all, file = "PG_sub.RDS")
PG_sub <- readRDS("PG_all.RDS")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#at this point, PG_all has a totla of 9316 cells 
  #higher than in the paper -> bc we chose lower resolution

#click on PG_all under environment
  #view metadata -> can only cluster of 7 or cluster of 10, no cluster of 9
  #go with 10, use resolution 0.5 (standard)
Idents(object = PG_all) <- PG_all$RNA_snn_res.0.5
DimPlot(PG_all, label = T, pt.size = 1, label.size = 6)

#shows supp figure 6n-o
plot <- VlnPlot(PG_all, features = c("Krt5","Krt14"), combine = FALSE, fill.by = c("feature","ident")) 
plot <-  lapply(
  X = plot,
  FUN = function(p) p + aes(color= PG_all$RNA_snn_res.0.5)
)
CombinePlots(plots = plot, legend = 'none')

################################################
# Doheatmap
listsA <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1','S100a4','Zeb1','Zeb2',
            'Twist1','Twist2','Snai1','Snai2','Prrx1','Prrx2','Foxc2','Cdh2','Acta2')
# Option 1
DoHeatmap(PG_all, features = listsA)
#one cluster with enrichment of the EMT genes


# Option 2: customise
DoHeatmap(PG_all, features = listsA, disp.min = -1,
          slot = 'scale.data', 
          group.colors = rainbow(9)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 1)
 #cluster 7 showed both epithelial and stromal markers

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Step 03 DEGs and GSEA analysis
  #normally should look into all clusters (FindAllMarkers), but cluster 7 is special?
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PG.markers.7 <- FindMarkers(PG_all, ident.1 = "7", min.pct = 0.20)
PG.markers.7$gene <- rownames(PG.markers.7)
gene_id <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")))
PG.markers.7 <- merge(PG.markers.7, gene_id[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
PG.markers.7.filt <- PG.markers.7
PG.markers.7.filt <- PG.markers.7.filt[!duplicated(PG.markers.7.filt$entrezgene_id),]
#PG.markers.7.filt <- PG.markers.7[which(PG.markers.7$p_val_adj < 0.05),]
genelist <- PG.markers.7$avg_log2FC
names(genelist) <-  PG.markers.7$entrezgene_id
genelist <- sort(genelist, decreasing = TRUE)

#fsgea package is used for gsea
  #check the enrichment in accordance to the Molecular Signature Database (MSigDB) version 6.2
mdf <- msigdbr(species="Mus musculus", 
                category="C2") 
mdf <- mdf %>% dplyr::select(gs_name, entrez_gene)
mdf_fgsea <- split(x=mdf$entrez_gene, 
                    f=mdf$gs_name)

#including more cell markers from the paper
  #this data was taken from Supplementary Data 1 Genes upregulated by ≥2 fold and FDR<0.05 in human benign prostatic basal epithelial cells compared to the corresponding luminal cells in the RNA-Seq analysis. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4773505/)
cellmarkers <- read.csv("../GSE111429/cellmarkers.csv", header = FALSE, fileEncoding="UTF-8-BOM")
