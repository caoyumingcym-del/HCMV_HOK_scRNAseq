if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Seurat")
BiocManager::install("decontX")
BiocManager::install("SoupX")


library(SoupX)
library(Seurat)
library(decontX)
library(ggplot2)


cmv.names <- read.table("../merlin.gene.names")
cmv.names <- c(cmv.names$V1, "EGFP")
cmv.names <- gsub("_","-",cmv.names)
cmv.names <- c(cmv.names, "RNA5.0","RNA4.9","RNA1.2","RNA2.7")


cmv.srt.filt <- readRDS("analysis_Robjects/cmv_seurat_filtered_counts_redo.rds")


samp <- unique(cmv.srt.filt$orig.ident)
samp <- samp[-grep("Mock", samp)]


###### run soupx
dirs <- Sys.glob("cellranger_nointron/Exp*/outs/")


for (i in 1:length(dirs)) {
  files <- dirs[i]
  print(files)
  counts <- load10X(files)
  if (i %in% c(9,12)) {
    counts = autoEstCont(counts,tfidfMin=0.7,soupQuantile=0.70)
  } else {
    counts = autoEstCont(counts,tfidfMin=0.8,soupQuantile=0.80)
  }
  
  counts.soupx = adjustCounts(counts)
  sampname <- strsplit(files, split = "/")[[1]][2]
  colnames(counts.soupx) <- gsub("1",sampname, colnames(counts.soupx))
  sobj <- CreateSeuratObject(counts = counts.soupx, 
                             project = sampname,
                             min.cells = 0, 
                             min.features = 0)
  if (i == 1) {
    cmv.srt.soupx <- sobj
  } else {
    cmv.srt.soupx <- merge(cmv.srt.soupx, sobj)
  }
}


table(cmv.srt.soupx$orig.ident)
table(cmv.srt.filt$orig.ident)


saveRDS(cmv.srt.soupx, file ="analysis_Robjects/cmv.srt.soupx.rds")


#### data QC
cmv.srt.soupx[["percent.mt"]] <- PercentageFeatureSet(cmv.srt.soupx, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(cmv.srt.soupx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.05, log = T)


plot1 <- FeatureScatter(cmv.srt.soupx, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1)
plot2 <- FeatureScatter(cmv.srt.soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
plot1 + plot2


cmv.srt.soupx.filt <- subset(cmv.srt.soupx, subset = nFeature_RNA > 2500  & percent.mt < 25)
table(cmv.srt.soupx.filt$orig.ident)
table(cmv.srt.filt$orig.ident)
saveRDS(cmv.srt.soupx.filt, file ="analysis_Robjects/cmv.srt.soupx.filt.rds")




cmv.g <- rownames(cmv.srt.soupx.filt)[which(rownames(cmv.srt.soupx.filt) %in% cmv.names)]
cmv.srt.soupx.filt[["percent.cmv"]] <- PercentageFeatureSet(cmv.srt.soupx.filt, features = cmv.g)
VlnPlot(cmv.srt.soupx.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.cmv"), 
        ncol = 3, pt.size = 0.05)
cmv.srt.soupx.filt[["cmv_exp"]] <- cmv.srt.soupx.filt[['percent.cmv']] * cmv.srt.soupx.filt[['nCount_RNA']] / 100
cmv.srt.soupx.filt[["cmv_exp_log"]] <-  log10(cmv.srt.soupx.filt[["cmv_exp"]] + 1)


ggplot(cmv.srt.soupx.filt@meta.data, aes(log10(cmv_exp+1), fill = orig.ident)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV UMI counts")


# add metadata
samp <- sapply(strsplit(colnames(cmv.srt.soupx.filt), split = "-"),'[[',2)
cmv.srt.soupx.filt[["strain"]] <- sapply(strsplit(samp, split = "_"),'[[',3)
cmv.srt.soupx.filt[["dpi"]] <- sapply(strsplit(samp, split = "_"),'[[',4)
cmv.srt.soupx.filt[["exp"]] <- sapply(strsplit(samp, split = "_"),'[[',2)




cmv.g <- rownames(cmv.srt.soupx.filt)[which(rownames(cmv.srt.soupx.filt) %in% cmv.names)]
cmv.srt.soupx.filt[["percent.cmv"]] <- PercentageFeatureSet(cmv.srt.soupx.filt, features = cmv.g)
VlnPlot(cmv.srt.soupx.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.cmv"), 
        ncol = 3, pt.size = 0.05)
cmv.srt.soupx.filt[["cmv_exp"]] <- cmv.srt.soupx.filt[['percent.cmv']] * cmv.srt.soupx.filt[['nCount_RNA']] / 100


halfmin_virus <- min(cmv.srt.soupx.filt$percent.cmv[cmv.srt.soupx.filt$percent.cmv>0])
halfmin_virus


cmv.srt.soupx.filt$percent.cmv.log10 <- log10(cmv.srt.soupx.filt$percent.cmv+halfmin_virus)
range(cmv.srt.soupx.filt$percent.cmv.log10)


cmv.srt.soupx.filt[["cmv_exp_log"]] <-  log10(cmv.srt.soupx.filt[["cmv_exp"]] + 1)
ggplot(cmv.srt.soupx.filt@meta.data, aes(log10(cmv_exp+1), fill = orig.ident)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV UMI counts")


ggplot(cmv.srt.soupx.filt@meta.data, aes(percent.cmv.log10, fill = orig.ident)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV expresion percentage log10")


saveRDS(cmv.srt.soupx.filt, file = "analysis_Robjects/cmv.srt.soupx.filt.rds")




## normalize data
cmv.srt.soupx.filt <- NormalizeData(cmv.srt.soupx.filt, normalization.method = "LogNormalize", 
                              scale.factor = 10000)


### dimenstion reduction on raw data before integration
# remove cmv genes in identifying variable genes
cmv.srt.soupx.filt <- FindVariableFeatures(cmv.srt.soupx.filt,
                                     selection.method = "vst", nfeatures = 2000)
VariableFeatures(cmv.srt.soupx.filt) <-  VariableFeatures(cmv.srt.soupx.filt)[-which(VariableFeatures(cmv.srt.soupx.filt) %in% cmv.names)]
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cmv.srt.soupx.filt), 10)


# plot variable features with and without labels
pdf("plots/vargenes.pdf", height = 10, width = 12)
plot1 <- VariableFeaturePlot(cmv.srt.soupx.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()


# scale data
all.genes <- rownames(cmv.srt.soupx.filt)
cmv.srt.soupx.filt <- ScaleData(cmv.srt.soupx.filt, features = all.genes)


### PCA
cmv.srt.soupx.filt <- RunPCA(cmv.srt.soupx.filt, features = VariableFeatures(object = cmv.srt.soupx.filt))


ElbowPlot(cmv.srt.soupx.filt, ndims = 50) 


### cluster cells
nPC = 20
cmv.srt.soupx.filt <- FindNeighbors(cmv.srt.soupx.filt, dims = 1:nPC)
cmv.srt.soupx.filt <- FindClusters(cmv.srt.soupx.filt, resolution = 0.5)


### Run UMAP
cmv.srt.soupx.filt <- RunUMAP(cmv.srt.soupx.filt, dims = 1:nPC)
DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by = c("orig.ident","strain","exp","dpi"), shuffle = T)


FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = "percent.cmv", max.cutoff = "q90", cols = c('grey',"blue","red"))
FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = "cmv_exp_log",  cols = c('grey',"blue","red","yellow"))


cmv.srt.soupx.filt[["virus.presence"]] <- "no"
cmv.srt.soupx.filt[["virus.presence"]][which(cmv.srt.soupx.filt[["percent.cmv"]] > 0),] <- "yes"
DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by = c("virus.presence"), cols = c("#808080","#FF0000"))


cmv.srt.soupx.filt[["percent.EGFP"]] <- PercentageFeatureSet(cmv.srt.soupx.filt, pattern = "^EGFP$")
cmv.srt.soupx.filt$GFP.presence <- "no"
cmv.srt.soupx.filt$GFP.presence[which(cmv.srt.soupx.filt$percent.EGFP> 0)] <- "yes"
DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by = c("GFP.presence"), cols = c("#808080","#FF0000"))
table(cmv.srt.soupx.filt$orig.ident, cmv.srt.soupx.filt$GFP.presence)
saveRDS(cmv.srt.soupx.filt, file = "analysis_Robjects/cmv.srt.soupx.filt.rds")




DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by = c("seurat_clusters"), shuffle = T)


FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = c("KRT16","KRT6B","TGM1","SPRR1B"))


####### markers for each cluster ######
p <- DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by = c("seurat_clusters"))
LabelClusters(plot = p, id = "seurat_clusters")


BiocManager::install("MAST")
BiocManager::install("GenomeInfoDb", force = T)
library(MAST)
cmv.srt.soupx.filt.new <- JoinLayers(cmv.srt.soupx.filt.new)
for ( i in 0:(length(unique(cmv.srt.soupx.filt.new$seurat_clusters))-1)) {
  clname <- paste("cluster",i, sep = "")
  marker <- FindMarkers(cmv.srt.soupx.filt.new, ident.1 = i, test.use = "MAST")
  outname <- paste("cluster_marker_soupx/", clname, "_marker_genes.tsv", sep = "")
  write.table(marker,file = outname, col.names = T, row.names = T, quote = F, sep = "\t")
}


for ( i in 0:(length(unique(cmv.srt.soupx.filt.new$seurat_clusters))-1)) {
  print(i)
  clname <- paste("cluster",i, sep = "")
  outname <-  paste("cluster_marker_soupx/", clname, "_marker_genes.tsv", sep = "")
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


FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = c("ISG15"))
FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = c("IFNB1"))
FeaturePlot(cmv.srt.soupx.filt, reduction = "umap", features = c("ISG15"), split.by = "Infection_state_bkgd", pt.size = 0.3)
DimPlot(cmv.srt.soupx.filt, reduction = "umap", group.by  = "Infection_state_bkgd")


cmv.srt.soupx.filt$Infection_state_bkgd <- cmv.srt.soupx.filt$Infection_state
cmv.srt.soupx.filt$Infection_state_bkgd[which(cmv.srt.soupx.filt$Infection_state == "Bystander" &
                                                cmv.srt.soupx.filt$virus.presence == "yes")] <- "Bystander_bkdg"
DimPlot(cmv.srt.soupx.filt, reduction = "umap", split.by = "Infection_state_bkgd")


VlnPlot(cmv.srt.soupx.filt, features = "ISG15", split.by = "Infection_state_bkgd", group.by = "orig.ident")
