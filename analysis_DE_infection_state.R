## =============================================================================
## Reproduce Fig. 4 edgeR DE results (*_edgeRDEresults.tsv) from cmv.srt.soupx.filt.rds
##
## Methods: edgeR (glmFit / glmLRT) via SignallingSingleCell findDEgenes() + edgeRDE().
##   HCMV genes removed; genes expressed (>0) in >10% of the tested cells kept;
##   significant = |log2FC| > 1 & FDR < 0.05.
##
## Fig. 4 comparisons (groups from Infection_state_bkgd):
##   Mock vs highly infected       -> 6_Infected_edgeRDEresults.tsv       (Infected vs Mock)
##   Mock vs marginally infected   -> 5_Bystander_bkgd_edgeRDEresults.tsv (Marginal vs Mock)
##   Mock vs bystander             -> 4_Bystander_edgeRDEresults.tsv      (Bystander vs Mock)
##   Bystander vs highly infected  -> 1_Infected_edgeRDEresults.tsv       (Infected vs Bystander)
##
## Each test: cells in `cells` are tested; contrastID cells vs all other tested
## cells (logFC = log2(contrastID / rest)), no batch term.
##
## Usage: Rscript analysis_DE_infection_state.R      (writes DE_repro/results/*.tsv + report)
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(Biobase)
  library(edgeR)
})
setwd("~/hcmv")
OUTDIR <- "DE_repro/results/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- SignallingSingleCell (ee0dac6), verbatim -------------------------------
construct_ex_sc <- function(input){
  HSMM_expr_matrix <- input
  num.columns <- ncol(HSMM_expr_matrix)
  name.columns <- colnames(HSMM_expr_matrix)
  num.rows <- nrow(HSMM_expr_matrix)
  name.rows <- rownames(HSMM_expr_matrix)
  HSMM_sample_sheet <- matrix(nrow = num.columns, ncol= 0)
  rownames(HSMM_sample_sheet) <- c(name.columns)
  HSMM_sample_sheet <- data.frame(HSMM_sample_sheet)
  HSMM_gene_annotation <- matrix(nrow = num.rows, ncol= 0)
  rownames(HSMM_gene_annotation) <- c(name.rows)
  HSMM_gene_annotation <- data.frame(HSMM_gene_annotation)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  ex_sc <- ExpressionSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  return(ex_sc)
}

edgeRDE <- function(input, groups, batch, sizefactor, lib_size, minCells,
                    pVal = 1, contrast = list(c(-1,1))) {
  ## remove any zero-variance genes
  rvar <- apply(exprs(input),1,var)
  idx.keep <- which(rvar>0)
  input <- input[idx.keep,]
  ## remove genes expressed in less than minCells fraction of the cells
  gfrac <- apply(exprs(input),1,function(a) length(a[a>0])/length(a))
  idx.keep <- which(gfrac>minCells)
  input <- input[idx.keep,]
  ## make DGElist
  y <- edgeR::DGEList(counts=exprs(input), group=groups)
  y$samples$lib.size <- lib_size          # library size BEFORE gene filtering
  if (is.null(sizefactor)) {
    message(sprintf('Calculating norm factors...'));flush.console()
    y <- edgeR::calcNormFactors(y)
  } else {
    y$samples$norm.factors = pData(input)[,sizefactor]/lib_size
  }
  if (length(unique(batch))>1) {
    design <- model.matrix(~0+batch+group, data=y$samples)
  } else {
    design <- model.matrix(~0+group, data=y$samples)
  }
  message('Estimating dispersion...');flush.console()
  y <- edgeR::estimateDisp(y, design)
  message('Doing likelihood ratio fit...');flush.console()
  fit <- edgeR::glmFit(y, design)
  tab <- list()
  if (length(unique(batch))>1) {
    lrt <- edgeR::glmLRT(fit)
    tab[['contrast_1']] <- as.data.frame(edgeR::topTags(lrt, p.value=pVal, n=Inf, sort.by='logFC'))
  } else {
    for (j in 1:length(contrast)) {
      cStr <- sprintf('contrast_%d', j)
      lrt <- edgeR::glmLRT(fit, contrast=contrast[[j]])
      tab[[cStr]] <- as.data.frame(edgeR::topTags(lrt, p.value=pVal, n=Inf, sort.by='logFC'))
    }
  }
  return(tab)
}

findDEgenes = function(input, pd = pData(input), DEgroup, contrastID, sizefactor = NULL,
                       lib_size = NULL, facet_by, minCells = 0.1, batchID = NULL,
                       outdir = "DEresults/", outsuffix = "DEresults.tsv", pVal = 1,
                       contrast = list(c(-1,1))) {
  for (i in 1:length(unique(pd[,facet_by]))) {
    name = sort(unique(as.character(pd[,facet_by])))[i]
    print(paste0("Performing DE for ", name))
    idx = rownames(pd)[which(pd[,facet_by]==name)]
    cntr = rownames(pd)[which(pd[,facet_by]==name & pd[,DEgroup]==contrastID)]
    z = as.matrix(exprs(input))[,idx]
    if(ncol(z) > 2){
      if (is.null(lib_size)) {
        lib_size = colSums(exprs(input))  # input must contain only the tested cells
      } else {
        lib_size = pd[idx,lib_size]
      }
      if (is.null(batchID)) {
        batch = as.factor(rep(1,ncol(input)))
      } else {
        batch = as.factor(pd[idx,batchID])
      }
      groupList = rep(0, times=ncol(z))
      groupList[which(colnames(z) %in% cntr)] = 1
      if(length(unique(groupList)) > 1){
        group = factor(groupList)
        z = construct_ex_sc(z)
        pData(z) = pd[idx,,drop=F]
        tab = edgeRDE(input = z, groups = group, batch = batch, sizefactor = sizefactor, lib_size = lib_size,
                      minCells = minCells, pVal = pVal, contrast = contrast)
        DEtablecntr = tab[['contrast_1']]
        outfile = paste(name, contrastID, outsuffix, sep = "_")
        write.table(DEtablecntr, paste(outdir, outfile, sep = ""), sep = "\t")
      } else {
        stop("Only one group in comparison")
      }
    }
  }
}

## ---- data: SoupX counts, HCMV genes removed ---------------------------------
obj <- readRDS("cmv.srt.soupx.filt.rds")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
md <- obj@meta.data
rm(obj); invisible(gc())

## same HCMV gene list construction as analysis_soupx_sample_processing.R
cmv.names <- read.table("merlin.gene.names")
cmv.names <- c(cmv.names$V1, "EGFP")
cmv.names <- gsub("_","-",cmv.names)
cmv.names <- c(cmv.names, "RNA5.0","RNA4.9","RNA1.2","RNA2.7")
counts <- counts[!rownames(counts) %in% cmv.names, ]

st <- as.character(md$Infection_state_bkgd)   # Mock / Bystander / Bystander_bkgd / Infected

## ---- the 8 tests ------------------------------------------------------------
## name = file prefix; cells = tested cells; grp = labels; id = contrast label
TESTS <- list(
  `1` = list(cells = st %in% c("Bystander", "Infected"),                 # Infected vs Bystander (no HCMV reads)
             grp = st, id = "Infected"),
  `2` = list(cells = st %in% c("Mock", "Bystander", "Infected"),         # Infected vs Mock + Bystander (no HCMV reads)
             grp = ifelse(st == "Infected", "Yes_CMV", "No_CMV"), id = "Yes_CMV"),  # NOT virus.presence: 1 Mock cell has HCMV reads
  # `3` = High_CMV: definition still unknown -- cells = virus.presence == "yes" (3326/3327 cells, 11114 genes)
  `4` = list(cells = st %in% c("Mock", "Bystander", "Infected"),         # Bystander vs Mock + Infected
             grp = st, id = "Bystander"),
  `5` = list(cells = st %in% c("Mock", "Bystander_bkgd"),                # Bystander_bkgd vs Mock
             grp = st, id = "Bystander_bkgd"),
  `6` = list(cells = st %in% c("Mock", "Infected"),                      # Infected vs Mock
             grp = st, id = "Infected"),
  `7` = list(cells = st %in% c("Bystander", "Infected"),                 # identical to test 1 (original files identical)
             grp = st, id = "Infected"),
  `8` = list(cells = md$Infection_state == "Bystander",                  # Bystander_bkgd vs Bystander (no HCMV reads)
             grp = st, id = "Bystander_bkgd")
)

for (nm in names(TESTS)) {
  t <- TESTS[[nm]]
  pd <- md[t$cells, ]
  pd$test <- nm
  pd$DEgrp <- t$grp[t$cells]
  input <- construct_ex_sc(as.matrix(counts[, rownames(pd)]))
  findDEgenes(input, pd = pd, DEgroup = "DEgrp", contrastID = t$id, facet_by = "test",
              minCells = 0.1, outdir = OUTDIR, outsuffix = "edgeRDEresults.tsv")
}

## ---- compare with original files --------------------------------------------
sig <- function(d) rownames(d)[abs(d$logFC) > 1 & d$FDR < 0.05]
report <- do.call(rbind, lapply(list.files(OUTDIR, "_edgeRDEresults.tsv$"), function(f) {
  a <- read.delim(file.path(OUTDIR, f)); b <- read.delim(f); g <- rownames(b)
  data.frame(file = f, same_genes = setequal(rownames(a), g),
             max_abs_dlogFC = signif(max(abs(a[g, "logFC"] - b$logFC)), 3),
             max_abs_dlogCPM = signif(max(abs(a[g, "logCPM"] - b$logCPM)), 3),
             sig_orig = length(sig(b)), sig_new = length(sig(a)),
             sig_shared = length(intersect(sig(a), sig(b))))
}))
print(report, row.names = FALSE)
write.table(report, "DE_repro/reproduction_report.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "DE_repro/sessionInfo.txt")
