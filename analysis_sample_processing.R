# CMV scRNAseq nextseq2000 data


setwd("/Users/yumingcao/Dropbox (UMass Medical School)/CMV_scRNAseq/nextseq2000_06_30_23/")


library(Seurat)


cmv.names <- read.table("../merlin.gene.names")
cmv.names <- c(cmv.names$V1, "EGFP")
cmv.names <- gsub("_","-",cmv.names)


dirs <- Sys.glob("cellranger_nointron/Exp*/outs/filtered_feature_bc_matrix/")
for (i in 1:length(dirs)) {
  files <- dirs[i]
  print(files)
  counts <- Read10X(files)
  sampname <- strsplit(files, split = "/")[[1]][2]
  colnames(counts) <- gsub("1",sampname, colnames(counts))
  sobj <- CreateSeuratObject(counts = counts, 
                             project = sampname,
                             min.cells = 0, 
                             min.features = 0)
  if (i == 1) {
    cmv.srt <- sobj
  } else {
    cmv.srt <- merge(cmv.srt, sobj)
  }
}


table(cmv.srt$orig.ident)




#### data QC
cmv.srt[["percent.mt"]] <- PercentageFeatureSet(cmv.srt, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(cmv.srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.05)


plot1 <- FeatureScatter(cmv.srt, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1)
plot2 <- FeatureScatter(cmv.srt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
plot1 + plot2


## filter out cells with more than 25% MT
cmv.srt.filt <- subset(cmv.srt, subset = nFeature_RNA > 2500  & percent.mt < 25)
dim(cmv.srt.filt)


VlnPlot(cmv.srt.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(cmv.srt.filt, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1)
plot2 <- FeatureScatter(cmv.srt.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
plot1 + plot2


# add metadata
samp <- sapply(strsplit(colnames(cmv.srt.filt), split = "-"),'[[',2)
cmv.srt.filt[["strain"]] <- sapply(strsplit(samp, split = "_"),'[[',3)
cmv.srt.filt[["dpi"]] <- sapply(strsplit(samp, split = "_"),'[[',4)
cmv.srt.filt[["exp"]] <- sapply(strsplit(samp, split = "_"),'[[',2)




cmv.g <- rownames(cmv.srt.filt)[which(rownames(cmv.srt.filt) %in% cmv.names)]
cmv.srt.filt[["percent.cmv"]] <- PercentageFeatureSet(cmv.srt.filt, features = cmv.g)
VlnPlot(cmv.srt.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.cmv"), 
        ncol = 3, pt.size = 0.05)
cmv.srt.filt[["cmv_exp"]] <- cmv.srt.filt[['percent.cmv']] * cmv.srt.filt[['nCount_RNA']] / 100


halfmin_virus <- min(cmv.srt.filt$percent.cmv[cmv.srt.filt$percent.cmv>0])
halfmin_virus


cmv.srt.filt$percent.cmv.log10 <- log10(cmv.srt.filt$percent.cmv+halfmin_virus)
range(cmv.srt.filt$percent.cmv.log10)




ggplot(cmv.srt.filt@meta.data, aes(log10(cmv_exp+1), fill = orig.ident)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV UMI counts")


ggplot(cmv.srt.filt@meta.data, aes(percent.cmv.log10, fill = orig.ident)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV expresion percentage log10")


saveRDS(cmv.srt.filt, file = "analysis_Robjects/cmv_seurat_filtered_counts.rds")
cmv.srt.filt <- readRDS("analysis_Robjects/cmv_seurat_filtered_counts.rds")


## normalize data
cmv.srt.filt <- NormalizeData(cmv.srt.filt, normalization.method = "LogNormalize", 
                              scale.factor = 10000)


### dimenstion reduction on raw data before integration
# remove cmv genes in identifying variable genes
cmv.srt.filt <- FindVariableFeatures(cmv.srt.filt,
                                     selection.method = "vst", nfeatures = 2000)
VariableFeatures(cmv.srt.filt) <-  VariableFeatures(cmv.srt.filt)[-which(VariableFeatures(cmv.srt.filt) %in% cmv.names)]
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cmv.srt.filt), 10)


# plot variable features with and without labels
pdf("plots/vargenes.pdf", height = 10, width = 12)
plot1 <- VariableFeaturePlot(cmv.srt.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()


# scale data
all.genes <- rownames(cmv.srt.filt)
cmv.srt.filt <- ScaleData(cmv.srt.filt, features = all.genes)


### PCA
cmv.srt.filt <- RunPCA(cmv.srt.filt, features = VariableFeatures(object = cmv.srt.filt))


ElbowPlot(cmv.srt.filt, ndims = 50) 


### cluster cells
nPC = 25
cmv.srt.filt <- FindNeighbors(cmv.srt.filt, dims = 1:nPC)
cmv.srt.filt <- FindClusters(cmv.srt.filt, resolution = 0.5)


### Run UMAP
cmv.srt.filt <- RunUMAP(cmv.srt.filt, dims = 1:nPC)
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi"), shuffle = T)
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi"))
FeaturePlot(cmv.srt.filt, reduction = "umap", features = "percent.cmv")


cmv.srt.filt[["virus.presence"]] <- "no"
cmv.srt.filt[["virus.presence"]][which(cmv.srt.filt[["percent.cmv"]] > 0),] <- "yes"
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi"))


DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("virus.presence"), cols = c("#808080","#FF0000"))


cmv.srt.filt[["percent.EGFP"]] <- PercentageFeatureSet(cmv.srt.filt, pattern = "^EGFP$")
cmv.srt.filt$GFP.presence <- "no"
cmv.srt.filt$GFP.presence[which(cmv.srt.filt$percent.EGFP> 0)] <- "yes"
table(cmv.srt.filt$orig.ident, cmv.srt.filt$GFP.presence)
saveRDS(cmv.srt.filt, file = "analysis_Robjects/cmv_seurat_filtered_counts.rds")


#########
cmv.srt.filt <- readRDS("analysis_Robjects/cmv_seurat_filtered_counts.rds")
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("virus.presence"), cols = c("#808080","#FF0000"))
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi"))


# plot gene expression
features <- c("SPRR1A", "SPRR1B","SPRR2A","SPRR2E", "A2ML1", "KRT14","KRT5")


FeaturePlot(cmv.srt.filt, reduction = "umap", features = features)
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi",
                                                       "seurat_clusters"))
DimPlot(cmv.srt.filt, reduction = "umap", group.by = c("seurat_clusters"))




#### identify marker genes for each cluster
BiocManager::install("MAST")
library(MAST)
for ( i in 0:(length(unique(cmv.srt.filt$seurat_clusters))-1)) {
  clname <- paste("cluster",i, sep = "")
  marker <- FindMarkers(cmv.srt.filt, ident.1 = i, test.use = "MAST")
  outname <- paste("cluster_marker/", clname, "_marker_genes.tsv", sep = "")
  write.table(marker,file = outname, col.names = T, row.names = T, quote = F, sep = "\t")
}




for ( i in 0:(length(unique(cmv.srt.filt$seurat_clusters))-1)) {
  print(i)
  clname <- paste("cluster",i, sep = "")
  outname <-  paste("cluster_marker/", clname, "_marker_genes.tsv", sep = "")
  marker <- read.table(outname, header = T, row.names = 1, sep = "\t")
  marker <- marker[order(marker$p_val_adj, -marker$pct.1, -marker$avg_log2FC),]
  len <- length(which(marker$avg_log2FC >0))
  if (len > 10) {
    gene <- head(rownames(marker)[which(marker$avg_log2FC >0)], 10)
  } else {
    gene <- rownames(marker)[which(marker$avg_log2FC >0)]
  }
  if (i ==0) {
    marker_list <- list(gene)
  } else {
    marker_list <- c(marker_list, list(gene))
  }
  
}
marker_list
names(marker_list)


FeaturePlot(cmv.srt.filt, reduction = "umap", features = marker_list[[10]])
FeaturePlot(cmv.srt.filt, reduction = "umap", features = marker_list[[11]])
FeaturePlot(cmv.srt.filt, reduction = "umap", features = marker_list[[12]])
