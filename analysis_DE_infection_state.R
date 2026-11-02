## =============================================================================
## Differential expression across HCMV infection states in primary HOKs
## Implements the edgeR DE step in the Methods of:
##   "Human Cytomegalovirus Infection of Primary Human Oral Keratinocytes Induces
##    Intermediate Keratinocyte Differentiation and an Altered Innate Immune
##    Response" (PMID 41000686; bioRxiv 10.1101/2025.09.16.676604)
##
## Input : analysis_Robjects/cmv.srt.soupx.filt.rds
##         (output of analysis_soupx_sample_processing.R)
##
## Methods being reproduced:
##   - Infection state per cell (Fig 3A/B):
##       Mock      : cells from Mock samples
##       Bystander : cells in infected samples with 0 viral UMIs
##       Marginal / High : virus+ cells split by a per-sample threshold at the
##                         first local minimum of the density of
##                         log10(% viral transcripts + min non-zero %)
##                         (BROOKELAB/SingleCell Run_filter_and_count_virus_cal07.R)
##   - DE (Fig 4A-D): edgeR; HCMV genes removed; genes detected in >=10% of the
##     tested cells; DEG = |log2FC| > 1 and FDR < 0.05
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
set.seed(1)

outdir <- "DE_infection_state"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

## ---- parameters (from Methods) ---------------------------------------------
min.pct  <- 0.10   # gene detected (count > 0) in >= 10% of tested cells
lfc.cut  <- 1      # |log2FC| > 1  (two-fold)
fdr.cut  <- 0.05

## Set TRUE to reuse an existing obj$Infection_state column
## (values must be Mock / Bystander / Marginal / High)
use.existing.states <- FALSE

## Optional manual per-sample thresholds on the log10(% viral + fudge) scale,
## to override the automatic local-minimum call, e.g.
## manual.thr <- c("Exp213_TB_3dpi" = -1.2)
manual.thr <- c()

## ---- load ------------------------------------------------------------------
obj <- readRDS("analysis_Robjects/cmv.srt.soupx.filt.rds")
obj <- JoinLayers(obj)

cmv.names <- read.table("../merlin.gene.names")$V1
cmv.names <- gsub("_", "-", c(cmv.names, "EGFP"))
cmv.names <- c(cmv.names, "RNA5.0", "RNA4.9", "RNA1.2", "RNA2.7")

md <- obj@meta.data
md$is_mock <- grepl("Mock", md$orig.ident)

## ---- infection state -------------------------------------------------------
if (use.existing.states && "Infection_state" %in% colnames(md)) {
  message("Using existing Infection_state column")
} else {
  # Threshold method from BROOKELAB/SingleCell Run_filter_and_count_virus_cal07.R:
  #   halfmin_virus <- min(virus_pct[virus_pct > 0])
  #   virus_pct_log10 <- log10(virus_pct + halfmin_virus)
  #   des <- density(virus_pct_log10)                       # default bandwidth
  #   mins <- des$x[which(diff(sign(diff(des$y))) == 2) + 1]
  #   infected = virus_pct_log10 > mins[1]
  # Applied per infected sample to virus+ cells (paper Fig 3A).
  fudge <- min(md$percent.cmv[md$percent.cmv > 0])   # = halfmin_virus
  md$pcmv_log10 <- log10(md$percent.cmv + fudge)

  find_local_min <- function(x) {
    des  <- density(x)
    mins <- des$x[which(diff(sign(diff(des$y))) == 2) + 1]
    if (length(mins) == 0) NA_real_ else mins[1]
  }

  inf.samples <- sort(unique(md$orig.ident[!md$is_mock]))
  all.mins <- list()
  thr <- sapply(inf.samples, function(s) {
    x <- md$pcmv_log10[md$orig.ident == s & md$percent.cmv > 0]
    des <- density(x)
    all.mins[[s]] <<- des$x[which(diff(sign(diff(des$y))) == 2) + 1]
    find_local_min(x)
  })
  print(all.mins)   # inspect: if mins[1] is a small dip, override in manual.thr
  if (length(manual.thr)) thr[names(manual.thr)] <- manual.thr
  if (any(is.na(thr))) {
    stop("No local minimum found for: ",
         paste(names(thr)[is.na(thr)], collapse = ", "),
         " -> set these in manual.thr")
  }

  thr.tab <- data.frame(sample = names(thr),
                        threshold_log10 = thr,
                        threshold_pct = 10^thr - fudge)
  print(thr.tab)
  write.table(thr.tab, file.path(outdir, "infection_thresholds.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  md$Infection_state <- "Mock"
  inf <- !md$is_mock
  md$Infection_state[inf & md$percent.cmv == 0] <- "Bystander"
  pos <- inf & md$percent.cmv > 0
  md$Infection_state[pos] <- ifelse(
    md$pcmv_log10[pos] > thr[md$orig.ident[pos]], "High", "Marginal")

  # QC plot: Fig 3A-style histograms with thresholds
  p <- ggplot(md[pos, ], aes(pcmv_log10)) +
    geom_histogram(bins = 50, fill = "grey40") +
    geom_vline(data = thr.tab, aes(xintercept = threshold_log10),
               colour = "red", linetype = 2) +
    facet_wrap(~orig.ident, scales = "free_y", nrow = 2,
               labeller = labeller(orig.ident = function(x) x)) +
    labs(x = "log10(% viral transcripts + fudge)", y = "cells") +
    theme_bw()
  ggsave(file.path(outdir, "plots", "infection_threshold_hist.pdf"),
         p, width = 12, height = 5)
}

md$Infection_state <- factor(md$Infection_state,
                             levels = c("Mock", "Bystander", "Marginal", "High"))
print(table(md$orig.ident, md$Infection_state))
obj$Infection_state <- md$Infection_state
saveRDS(obj, "analysis_Robjects/cmv.srt.soupx.filt.infstate.rds")

## ---- host count matrix (HCMV genes removed) ---------------------------------
counts.all  <- GetAssayData(obj, assay = "RNA", layer = "counts")
host.genes  <- setdiff(rownames(counts.all), cmv.names)
counts.host <- counts.all[host.genes, ]
lib.host    <- colSums(counts.host)   # library size on host genes only
rm(counts.all); gc()

## ---- edgeR DE ----------------------------------------------------------------
run_edger <- function(case, ref) {
  cells <- rownames(md)[md$Infection_state %in% c(case, ref)]
  grp   <- factor(as.character(md[cells, "Infection_state"]), levels = c(ref, case))
  cnt   <- counts.host[, cells, drop = FALSE]

  keep <- rowSums(cnt > 0) >= min.pct * length(cells)
  cnt  <- cnt[keep, ]

  pct.case <- rowSums(cnt[, grp == case] > 0) / sum(grp == case)
  pct.ref  <- rowSums(cnt[, grp == ref]  > 0) / sum(grp == ref)

  y <- DGEList(counts = as.matrix(cnt), group = grp,
               lib.size = lib.host[cells])
  y <- calcNormFactors(y, method = "TMMwsp")   # TMMwsp copes with sparse cells
  design <- model.matrix(~ grp)
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)

  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$gene     <- rownames(tab)
  tab$pct.case <- pct.case[tab$gene]
  tab$pct.ref  <- pct.ref[tab$gene]
  tab$DE <- "NS"
  tab$DE[tab$logFC >  lfc.cut & tab$FDR < fdr.cut] <- "Up"
  tab$DE[tab$logFC < -lfc.cut & tab$FDR < fdr.cut] <- "Down"
  tab[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR",
          "pct.case", "pct.ref", "DE")]
}

volcano <- function(tab, title, n.lab = 20) {
  tab$nlog10FDR <- -log10(pmax(tab$FDR, .Machine$double.xmin))
  lab <- bind_rows(
    tab %>% filter(DE == "Up")   %>% arrange(FDR, desc(logFC)) %>% head(n.lab),
    tab %>% filter(DE == "Down") %>% arrange(FDR, logFC)       %>% head(n.lab))
  ggplot(tab, aes(logFC, nlog10FDR, colour = DE)) +
    geom_point(size = 0.6) +
    geom_vline(xintercept = c(-lfc.cut, lfc.cut), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(fdr.cut), linetype = 2, colour = "grey50") +
    geom_text_repel(data = lab, aes(label = gene), size = 2.5,
                    max.overlaps = Inf, show.legend = FALSE) +
    scale_colour_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    labs(title = title, x = "log2 fold change", y = "-log10 FDR") +
    theme_bw()
}

## case vs reference (positive logFC = higher in case), as in Fig 4A-D
comps <- list(
  High_vs_Mock          = c("High",      "Mock"),
  Marginal_vs_Mock      = c("Marginal",  "Mock"),
  Bystander_vs_Mock     = c("Bystander", "Mock"),
  High_vs_Bystander     = c("High",      "Bystander"),
  Marginal_vs_Bystander = c("Marginal",  "Bystander")   # reported as no DEGs
)

de.res <- list()
for (nm in names(comps)) {
  message("edgeR: ", nm)
  tab <- run_edger(comps[[nm]][1], comps[[nm]][2])
  de.res[[nm]] <- tab
  write.table(tab, file.path(outdir, paste0("DE_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(outdir, "plots", paste0("volcano_", nm, ".pdf")),
         volcano(tab, nm), width = 7, height = 6)
  message(sprintf("  %s: %d up, %d down", nm,
                  sum(tab$DE == "Up"), sum(tab$DE == "Down")))
}
saveRDS(de.res, file.path(outdir, "DE_results_list.rds"))

## ---- overlaps among the three vs-Mock comparisons (Fig 4I/J) ---------------
vs.mock <- c("High_vs_Mock", "Marginal_vs_Mock", "Bystander_vs_Mock")
for (dir in c("Up", "Down")) {
  sets <- lapply(de.res[vs.mock], function(t) t$gene[t$DE == dir])
  all.g <- unique(unlist(sets))
  memb <- data.frame(gene = all.g,
                     sapply(sets, function(s) all.g %in% s), check.names = FALSE)
  memb$pattern <- apply(memb[, vs.mock], 1,
                        function(r) paste(vs.mock[r], collapse = "&"))
  write.table(memb, file.path(outdir, paste0("overlap_vsMock_", dir, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(table(memb$pattern))
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
