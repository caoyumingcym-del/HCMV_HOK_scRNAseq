# HCMV infection of primary human oral keratinocytes: scRNA-seq analysis

Analysis code for:

> Cojohari O, Cao Y, Temirbek A, Garber M, Colubri A, Kowalik T. **Human Cytomegalovirus Infection of Primary Human Oral Keratinocytes Induces Intermediate Keratinocyte Differentiation and an Altered Innate Immune Response.** *bioRxiv* (2025). https://doi.org/10.1101/2025.09.16.676604

This repository reproduces the single-cell analyses behind **Figures 2, 3 and 4**.

## Data

- Raw and processed data: GEO [GSE348416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348416)
- Cells: primary human oral keratinocytes (HOK), mock-infected or infected with HCMV strains TB40/E-GFP (`TB`) or MOLD (`MOLD`), sampled at several days post-infection (dpi)
- Reads were aligned with Cell Ranger (introns excluded) to a combined human + HCMV Merlin + EGFP reference, so viral transcripts are counted alongside host genes

## Repository contents

| Script | Purpose |
|---|---|
| `analysis_sample_processing.R` | Loads the Cell Ranger filtered matrices, merges samples, runs QC and filtering, then normalization, PCA, clustering, UMAP and MAST marker genes per cluster. Calculates viral load per cell (`percent.cmv`, `cmv_exp`) and flags virus and GFP detection. |
| `analysis_soupx_sample_processing.R` | Runs the same pipeline after removing ambient RNA with **SoupX**. SoupX removes viral reads that leak into uninfected cells, so these are the counts used for infection calls. Also produces cluster markers and ISG15/IFNB1 expression by infection state. |
| `analysis_cmv_gene_exp_soupx.R` | Classifies cells as infected or not infected and models how viral genes change with viral load (see below). |
| `analysis_DE_infection_state.R` | Runs the edgeR differential expression between infection states for Figure 4 and compares the results with the published tables. |

## Analysis overview

1. **QC**: keep cells with more than 2,500 detected genes and less than 25% mitochondrial reads.
2. **Ambient RNA correction**: SoupX `autoEstCont` (`tfidfMin`/`soupQuantile` = 0.8 in most samples, 0.7 in two), then `adjustCounts`.
3. **Host clustering**: LogNormalize, then 2,000 variable genes with viral genes removed, then PCA (20–25 PCs), Louvain clustering (resolution 0.5) and UMAP.
4. **Viral load**: the viral share of each cell's UMIs (`percent.cmv`), log10-transformed with a pseudocount of half the smallest non-zero value.
5. **Infection calls**: within each infected sample, a cell counts as `Infected` if its log10 viral fraction exceeds the second local minimum of the density curve (`Infection_localminima`). An alternative method, based on the distribution of UL123 (IE1) expression, is also included.
6. **Viral gene trajectory**: infected cells are grouped into 5% viral-load bins and averaged within each bin. Viral genes in the top 30% for both mean and SD are kept, each is scaled to its maximum, and k-means (k = 5) groups them into expression programs. The results are plotted as a heatmap and as one trajectory plot per cluster.
7. **Differential expression (Figure 4)**: edgeR likelihood-ratio tests (`glmFit`/`glmLRT`) on the SoupX counts, using the `findDEgenes`/`edgeRDE` functions from [SignallingSingleCell](https://github.com/garber-lab/SignallingSingleCell). The script includes the pre-June-2020 version of these functions, in which the library size is each cell's total count before gene filtering. HCMV genes are removed, only genes detected in more than 10% of the tested cells are kept, and differentially expressed genes are those with |log2FC| > 1 and FDR < 0.05. Cell groups come from `Infection_state_bkgd`: `Mock`, `Bystander` (infected sample, no viral reads), `Bystander_bkgd` (marginally infected: viral reads below the infection cutoff) and `Infected` (highly infected).

   | Comparison | Output table | Cells tested |
   |---|---|---|
   | Mock vs highly infected | `6_Infected_edgeRDEresults.tsv` | Infected vs Mock |
   | Mock vs marginally infected | `5_Bystander_bkgd_edgeRDEresults.tsv` | Bystander_bkgd vs Mock |
   | Mock vs bystander (Fig. 4D) | `4_Bystander_edgeRDEresults.tsv` | Bystander vs Mock and Infected pooled |
   | Bystander vs highly infected | `1_Infected_edgeRDEresults.tsv` | Infected vs Bystander |

   In each table, `logFC` is the log2 ratio of the first group listed under "Cells tested" over the rest.

## Requirements

R (≥ 4.3) with:

- **Bioconductor/CRAN:** `Seurat` (v5), `SoupX`, `decontX`, `MAST`, `ggplot2`, `dplyr`, `tidyr`, `scales`, `viridis`, `RColorBrewer`
- **Distribution fitting (UL123-based calls):** `fitdistrplus`, `MASS`, `actuar`, `distrMod`
- **Differential expression:** `edgeR`, `Biobase`

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Seurat", "SoupX", "decontX", "MAST", "edgeR", "Biobase"))
install.packages(c("ggplot2", "dplyr", "tidyr", "scales", "viridis",
                   "RColorBrewer", "fitdistrplus", "actuar", "distrMod"))
```

## Expected directory layout

The scripts use relative paths and expect to run from a project directory structured like this:

```
project/
├── cellranger_nointron/
│   └── Exp<N>_<exp>_<strain>_<dpi>/outs/   # one Cell Ranger output per sample
├── analysis_Robjects/                       # intermediate .rds files
├── plots/
├── cluster_marker/  cluster_marker_soupx/
../merlin.gene.names                          # HCMV Merlin gene list, one per line
```

Sample metadata (`exp`, `strain`, `dpi`) comes from the sample folder names, so keep the `Exp_<exp>_<strain>_<dpi>` naming.

## Running

Run the scripts in this order:

```
1. analysis_sample_processing.R
2. analysis_soupx_sample_processing.R
3. analysis_cmv_gene_exp_soupx.R
4. analysis_DE_infection_state.R
```

Scripts 1–3 are interactive analysis scripts, so run them section by section in RStudio rather than with `Rscript`. Update the `setwd()` path in script 1 before you start. Script 4 runs from start to finish with `Rscript analysis_DE_infection_state.R`, which takes about an hour. Before running it, update its `setwd()` path to the directory that holds `cmv.srt.soupx.filt.rds`, `merlin.gene.names` and the published `*_edgeRDEresults.tsv` tables.

## Main outputs

- `analysis_Robjects/cmv.srt.soupx.filt.rds`: SoupX-corrected, filtered Seurat object with viral load and infection-state metadata
- `plots/cmv_exp_histogram.pdf`: per-sample viral-fraction distributions with infection calls
- `plots/heatmap_bin_cmv_exp_5_soupx.pdf`: viral gene expression across viral-load bins
- `plots/cluster_trajectory_cmv.pdf`: viral gene expression programs (k-means clusters) plotted against viral load
- `cluster_marker*/`: MAST marker genes for each cluster
- `DE_repro/results/*_edgeRDEresults.tsv`: edgeR results for each comparison
- `DE_repro/reproduction_report.tsv`: comparison with the published tables (same tested genes and significant genes; logFC differs by about 2e-4 because of the edgeR version)

## Citation

If you use this code, please cite the bioRxiv preprint above.
