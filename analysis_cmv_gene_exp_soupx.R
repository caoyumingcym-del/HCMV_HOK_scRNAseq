library(ggplot2)
library(Seurat)


cmv.srt.soupx.filt <- readRDS("analysis_Robjects/cmv.srt.soupx.filt.rds")


cmv.names <- read.table("../merlin.gene.names")$V1
cmv.names <- gsub("_","-", cmv.names)
cmv.names <- c(cmv.names, "EGFP", "RNA5.0","RNA4.9","RNA1.2","RNA2.7")


idx <- which(rownames(cmv.srt.soupx.filt) %in% cmv.names)


cmv.srt.soupx.filt.new <- JoinLayers(cmv.srt.soupx.filt)


cmv.count <- GetAssayData(cmv.srt.soupx.filt.new, layer = "RNA", slot = "counts")
cmv.count <- cmv.count[idx,]
cmv.count <- cmv.count[,which(cmv.srt.soupx.filt$strain %in% c("TB","MOLD") & cmv.srt.soupx.filt$cmv_exp > 0)]


dim(cmv.count)
dim(cmv.srt.soupx.filt)
head(cmv.srt.soupx.filt@meta.data)


rsum <- apply(cmv.count, 1, sum)
r_idx <- which(rsum == 0)


cmv.count <- cmv.count[-r_idx,]
dim(cmv.count)


# convert to log10 percentage
vload <-  cmv.srt.soupx.filt$cmv_exp[match(colnames(cmv.count), colnames(cmv.srt.soupx.filt))]
vload_fudge <- vload + 100
vpercent <- cmv.srt.soupx.filt$percent.cmv[match(colnames(cmv.count), colnames(cmv.srt.soupx.filt))]
#cmv.count.pct <- sweep(cmv.count, 2, vload, FUN = '/') * 100
cmv.count.pct <- sweep(cmv.count, 2, vload_fudge, FUN = '/') * 100 ## use fudge factor


csum <- colSums(cmv.count.pct)
cmv.count.pct.df <- as.data.frame(t(as.matrix(cmv.count.pct)))


colnames(cmv.count.pct.df) <- gsub("-\\(","_",colnames(cmv.count.pct.df)  )
colnames(cmv.count.pct.df) <- gsub("\\)","",colnames(cmv.count.pct.df)  )
colnames(cmv.count.pct.df) <- gsub("-","_",colnames(cmv.count.pct.df)  )
cmv.count.pct.df$cell <- rownames(cmv.count.pct.df)


## order cells 
cmv.count.pct.df$vload <- vload
cmv.count.pct.df$vpercent <- vpercent
cmv.count.pct.df <- cmv.count.pct.df[order(cmv.count.pct.df$vpercent, decreasing = F),] # order by viral percentage
cmv.count.pct.df$Infection_localminima <- cmv.srt.soupx.filt$Infection_localminima[rownames(cmv.count.pct.df)]
hist(log10(cmv.count.pct.df$vpercent), breaks = 100)


cmv.count.pct.df.sub <- cmv.count.pct.df
hist(log10(cmv.count.pct.df.sub$vpercent), breaks = 100)


# try with all cells first
cmv.count.pct.df$cell <- factor(cmv.count.pct.df$cell, levels =cmv.count.pct.df$cell )
cmv.count.pct.df.sub$cell <- factor(cmv.count.pct.df.sub$cell, levels =cmv.count.pct.df.sub$cell )


# scale by gene between 0 and 1 and # cluster by gene
library(scales)
#cmv.count.scale.sub <- as.data.frame(rescale(as.matrix(cmv.count.pct.df.sub[,1:344]))) # EGFP removed
cmv.count.scale.sub <- as.data.frame(rescale(as.matrix(cmv.count.pct.df.sub[,1:345]))) # EGFP included
# remove lowly expressed cmv genes
range.df <- data.frame("gene" = colnames(cmv.count.scale.sub))
for(i in 1:ncol(cmv.count.scale.sub)) {
  rg = range(cmv.count.scale.sub[,i])
  range.df[i,2] <- rg[1]
  range.df[i,3] <- rg[2]
}


ggplot(range.df, aes(V3)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = mean(range.df$V3)) +
  geom_vline(xintercept = median(range.df$V3),linetype = "dashed")


keep.cmvgene <- range.df$gene[which(range.df$V3 > mean(range.df$V3))]
# write.table(keep.cmvgene, file = "cmv_exp_pattern/cmv_gene_high_expression.tsv", 
#             quote = F, col.names = F, row.names = F)
cmv.count.scale.sub <- cmv.count.scale.sub[,which(colnames(cmv.count.scale.sub)%in% keep.cmvgene)]
which(rowSums(cmv.count.scale.sub) == 0 )


d <- dist(t(cmv.count.scale.sub))
hc <- hclust(d, method = "average")


# k means
#kmeans_res <- kmeans(t(cmv.count.scale.sub), centers = 3,nstart = 25)
#kmeans_order <- order(kmeans_res$cluster)


cmv_order <- colnames(cmv.count.scale.sub)[hc$order]
cmv_out <- factor(cmv_order, levels = cmv_order)
# saveRDS(cmv_out, file = "cmv_exp_pattern/cmv_gene_high_expression_factor.rds")
cmv.count.pct.df.scale <- cmv.count.scale.sub
cmv.count.pct.df.scale$rank <- seq(1,nrow(cmv.count.pct.df.scale),1)
cmv.count.pct.df.scale$vload <- cmv.count.pct.df.sub$vload
cmv.count.pct.df.scale$cell <- cmv.count.pct.df.sub$cell
cmv.count.pct.df.scale$vpercent <- cmv.count.pct.df.sub$vpercent
cmv.count.pct.df.scale$Infection_localminima <- cmv.count.pct.df.sub$Infection_localminima
cmv.count.pct.df.scale.lng <- tidyr::gather(cmv.count.pct.df.scale, cmv_gene, exp,
                                            -cell, - vload, -vpercent, -rank, -Infection_localminima)
cmv.count.pct.df.scale.lng$cmv_gene <- factor(cmv.count.pct.df.scale.lng$cmv_gene , levels = cmv_order)
cmv.count.pct.df.scale.lng$dpi <- gsub(".*_","",cmv.count.pct.df.scale.lng$cell)
cmv.count.pct.df.scale.lng$sample <- gsub(".*-","", cmv.count.pct.df.scale.lng$cell)
cmv.count.pct.df.scale.lng.tmp <- cmv.count.pct.df.scale.lng
cmv.count.pct.df.scale.lng.tmp$exp[which(cmv.count.pct.df.scale.lng.tmp$exp > 0.1)] <- 0.1
cmv.count.pct.df.scale.lng.tmp$height <- 1


library(ggplot2)
library(viridis)
ggplot(cmv.count.pct.df.scale.lng.tmp, aes(cell, cmv_gene, fill = exp)) +
  geom_tile() +
#  facet_wrap(~Infection_localminima, scales = "free_x") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_viridis()
ggsave("plots/cmv_fraction_hm.pdf", height = 10, width = 8)


early_genes <- c("ORFL264C_UL123", "ORFL253W_UL112","ORFL268C")
g <- "ORFL264C_UL123"
g2 <- "RNA2.7"
g3<- "ORFL253W_UL112"
g4 <- "ORFL36W_UL6"
cmv.count.pct.df.scale.lng.g <- cmv.count.pct.df.scale.lng[which(cmv.count.pct.df.scale.lng$cmv_gene %in% c(g,g2,g3,g4)),]
cmv.count.pct.df.scale.lng.g$order <- rep(seq(1, nrow(cmv.count.pct.df.sub), 1),4)
ggplot(cmv.count.pct.df.scale.lng.g, aes(order, exp, color = cmv_gene)) +
  geom_point(size = 0.5)+
  geom_smooth() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle(paste(g, g3))


#########


cmv.count.pct.df.scale.lng$order <- rep(seq(1, nrow(cmv.count.pct.df.sub), 1),length(keep.cmvgene))
ggplot(cmv.count.pct.df.scale.lng, aes(rank, exp)) +
  geom_point(size = 0.5)+
  geom_smooth() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_vline(xintercept = 2073)


##### fit a distribution ########


library(fitdistrplus)
library(stats4)
library(MASS)
# for other necessary test or graphical tools
library(survival)
BiocManager::install('actuar')
BiocManager::install('distrMod')
library(actuar)
library(distrMod)


plotdist(cmv.count.pct.df.scale.lng$exp, histo = TRUE, demp = TRUE)




cmv.count.pct.df.scale.lng$rank_group <- cut(cmv.count.pct.df.scale.lng$rank, breaks=10)
ggplot(cmv.count.pct.df.scale.lng, aes(x=rank_group, y=exp)) +
  geom_boxplot() +
  labs(title="Boxplot of y by Rank Group", x="Rank Group", y="y")


cmv.count.pct.df.scale.lng.g <- cmv.count.pct.df.scale.lng[which(cmv.count.pct.df.scale.lng$cmv_gene %in% c(g)),]
cmv.count.pct.df.scale.lng.g$rank_group <- cut(cmv.count.pct.df.scale.lng.g$rank, breaks=10)


ggplot(cmv.count.pct.df.scale.lng.g, aes(x=rank_group, y=exp)) +
  geom_boxplot() +
  labs(title="Boxplot of y by Rank Group", x="Rank Group", y="y")


ggplot(cmv.count.pct.df.scale.lng.g, aes(exp)) +
  geom_histogram(bins = 100)


library(dplyr)


early_g.mean <- cmv.count.pct.df.scale.lng.g %>%
  group_by(cell) %>%
  summarise_at(vars(exp), funs(mean(., na.rm=TRUE)))
early_g.mean <- as.data.frame(early_g.mean)
early_g.mean <- cmv.count.pct.df.scale.lng.g
early_g.mean <- merge(early_g.mean, cmv.count.pct.df.scale[,c(64,66)], by= "cell", all = 1)
head(early_g.mean)


range(early_g.mean$exp)
ggplot(early_g.mean, aes(x=rank, y=exp)) +
  geom_point(size = 0.5) 


early_g.mean$exp_est <- round(early_g.mean$exp*100)


vec <- c()
for (i in 1:nrow(early_g.mean)) {
  new <- rep(early_g.mean$rank[i], early_g.mean$exp_est[i])
  vec <- c(vec, new)
}


hist(vec, breaks = 100)
vec_scaled <- vec / 10000
library(fitdistrplus)
# Fit different distributions
fit_norm <- fitdist(vec, "norm")
fit_lognorm <- fitdist(vec, "lnorm")
fit_exp <- fitdist(vec, "exp")
fit_gamma <- fitdist(vec, "gamma")
fit_negbin <- fitdist(vec,"nbinom")
fit_binom <- fitdist(vec, "binom")
fit_wb <- fitdist(vec,"weibull")
fit_poisson <- fitdist(vec, "pois", method = "mle")
fit_beta <- fitdist(vec_scaled, "beta")


# Summary of fits
summary(fit_norm)
summary(fit_lognorm)
summary(fit_gamma)
summary(fit_negbin)
summary(fit_poisson)
summary(fit_beta)


# Compare distributions
plot.legend <- c("Normal",  "Log-normal",  "Gamma", "Negative Binomial", "Poisson", "Weibull")
denscomp(list(fit_norm,  fit_lognorm, fit_gamma,fit_negbin, fit_poisson, fit_wb), legendtext = plot.legend)
qqcomp(list(fit_norm, fit_lognorm,  fit_gamma,fit_negbin, fit_poisson, fit_wb), legendtext = plot.legend)
cdfcomp(list(fit_norm, fit_lognorm,  fit_gamma,fit_negbin, fit_poisson, fit_wb), legendtext = plot.legend)
ppcomp(list(fit_norm, fit_lognorm,  fit_gamma,fit_negbin, fit_poisson, fit_wb), legendtext = plot.legend)




plot.legend <- c( "Beta")
denscomp(list(fit_beta), legendtext = plot.legend)
qqcomp(list(fit_beta), legendtext = plot.legend)
cdfcomp(list(fit_beta), legendtext = plot.legend)
ppcomp(list(fit_beta), legendtext = plot.legend)






# Goodness-of-Fit statistics
gofstat(list(fit_norm, fit_lognorm,  fit_gamma, fit_negbin))


# Validate the best fit (e.g., Normal distribution)
best_fit <- fit_norm
plot(best_fit)
ks.test(data$y, "pnorm", mean=mean(data$y), sd=sd(data$y))


fit_mean <- best_fit$estimate['mean']
fit_sd <- best_fit$estimate['sd']


bottom_5th_quantile <- qnorm(0.05, mean = fit_mean, sd = fit_sd)
bottom_5th_quantile


## define infection by UL123 distribution
cmv.count.pct.df.scale.lng.g$UL123_define_infection <- "Infected"
cmv.count.pct.df.scale.lng.g$UL123_define_infection[which(cmv.count.pct.df.scale.lng.g$rank < bottom_5th_quantile)] <- "NotInfected"


cmv.srt.soupx.filt$UL123_define_infection <- "NotInfected"
infected_Cell <- as.character(cmv.count.pct.df.scale.lng.g$cell[which(cmv.count.pct.df.scale.lng.g$UL123_define_infection == "Infected")])
cmv.srt.soupx.filt$UL123_define_infection[which(colnames(cmv.srt.soupx.filt) %in% infected_Cell)] <- "Infected"


table(cmv.srt.soupx.filt$UL123_define_infection, cmv.srt.soupx.filt$Infection_localminima)
ggplot(cmv.srt.soupx.filt@meta.data, aes(percent.cmv.log10, fill = UL123_define_infection)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV expresion percentage log10")
saveRDS(cmv.srt.soupx.filt, file = "analysis_Robjects/cmv.srt.soupx.filt.rds")


### define infection by local minima ######
samples <- unique(cmv.srt.soupx.filt$orig.ident)
samples <- samples[-grep("Mock", samples)]
cmv.srt.soupx.filt$Infection_localminima <- "NotInfected"
for(i in 1:length(samples)) {
  print(samples[i])
  dat <- cmv.srt.soupx.filt@meta.data[which(cmv.srt.soupx.filt$orig.ident == samples[i]),]
  des.all <- density(dat$percent.cmv.log10)
  min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
  cutoff = min.all[2]
  cmv.srt.soupx.filt$Infection_localminima[which(cmv.srt.soupx.filt$orig.ident == samples[i] & 
                                                   cmv.srt.soupx.filt$percent.cmv.log10 > cutoff )] <- "Infected"
}


ggplot(cmv.srt.soupx.filt@meta.data, aes(percent.cmv.log10, fill = Infection_localminima)) +
  geom_histogram( bins = 50) +
  facet_wrap(~orig.ident, scales = "free_y", nrow=2) +
  ggtitle("CMV expresion percentage log10") +
  theme_classic()
ggsave("plots/cmv_exp_histogram.pdf", height = 5, width = 12)
saveRDS(cmv.srt.soupx.filt, file = "analysis_Robjects/cmv.srt.soupx.filt.rds")


cmv.srt.soupx.filt$cell <- rownames(cmv.srt.soupx.filt@meta.data)
cmv.count.pct.df.scale.lng <- merge(cmv.count.pct.df.scale.lng, cmv.srt.soupx.filt@meta.data[,c(17,18)], by = "cell")


ggplot(cmv.count.pct.df.scale.lng, aes(rank, exp, color = Infection_localminima)) +
  geom_point(size = 0.5)+
  geom_smooth() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 




#########
head(cmv.count.pct.df.scale.lng)
cmv.count.pct.df.scale.lng.uninfect <- cmv.count.pct.df.scale.lng[which(cmv.count.pct.df.scale.lng$Infection_localminima == "NotInfected" & 
                                                                          cmv.count.pct.df.scale.lng$exp > 0),]


dim(cmv.count.pct.df.scale.lng.uninfect)


#########
#### bin vpercent by 5% increment, create meta cell by averaging the cmv count #####
cmv.names <- cmv.names[-grep("EGFP",cmv.names)]
idx <- which(rownames(cmv.srt.soupx.filt) %in% cmv.names)


cmv.srt.soupx.filt.new <- JoinLayers(cmv.srt.soupx.filt)
cmv.count <- GetAssayData(cmv.srt.soupx.filt.new, layer = "RNA", slot = "counts")
cmv.count <- cmv.count[idx,]
cmv.count <- cmv.count[,which(cmv.srt.soupx.filt$strain %in% c("TB","MOLD") & cmv.srt.soupx.filt$cmv_exp > 0)]


dim(cmv.count)
dim(cmv.srt.soupx.filt)
head(cmv.srt.soupx.filt@meta.data)


rsum <- apply(cmv.count, 1, sum)
r_idx <- which(rsum == 0)
cmv.count <- cmv.count[-r_idx,]
dim(cmv.count)
csum.new <- apply(cmv.count, 2, sum)


# convert to log10 percentage
vload <-  cmv.srt.soupx.filt$cmv_exp[match(colnames(cmv.count), colnames(cmv.srt.soupx.filt))]
vpercent <- cmv.srt.soupx.filt$percent.cmv[match(colnames(cmv.count), colnames(cmv.srt.soupx.filt))]
cmv.count.df <- as.data.frame(t(as.matrix(cmv.count)))
cmv.count.df$vpercent <- vpercent


infected_Cells <- colnames(cmv.srt.soupx.filt)[which(cmv.srt.soupx.filt$Infection_localminima == "Infected")]
cmv.count.df <- cmv.count.df[which(rownames(cmv.count.df) %in% infected_Cells),]


cuts <-cmv.count.df$vpercent/5
cuts.floor <- floor(cuts)
range(cuts.floor)
cmv.count.df$bin_2percent <- cuts.floor * 5
cmv.count.df <- cmv.count.df[,-c(345)]


## remove lowly expressed genes 
head(cmv.count.df)
cmv.mean <- apply(cmv.count.df[,-345],2,mean)
cmv.sd <- apply(cmv.count.df[,-345],2,sd)


plot(cmv.mean, cmv.sd, bty='n', pch=19, log="xy",  col = alpha("black", 0.1))


filt.mean <- which(cmv.mean > quantile(cmv.mean, 0.7))
filt.sd <- which(cmv.sd > quantile(cmv.sd, 0.7))
keep.idx <- intersect(filt.mean, filt.sd)
cmv.count.df <- cmv.count.df[,c(keep.idx, 345)]
ggplot(cmv.count.df, aes(bin_2percent)) +
  geom_histogram()


library(dplyr)
cmv.count.df.bins <- cmv.count.df %>% 
  group_by(bin_2percent) %>% 
  summarize(across(everything(), mean))


hist_table <- as.data.frame(table(cmv.count.df$bin_2percent))
colnames(hist_table)[1] <- "bin"


ggplot(hist_table, aes(bin, Freq)) + geom_bar(stat = "identity") +
  scale_y_log10() +
  theme_classic()
ggsave("plots/hist_bin_cmv_exp_5_soupx.pdf", height = 2, width = 6)


## relative gene exp by max exp per gene
max_val <- apply(cmv.count.df.bins[,-1], 2, max)
cmv.count.frac.bins <- sweep(cmv.count.df.bins[,-1], 2, max_val, FUN = "/")


cmv.count.frac.bins <- t(cmv.count.frac.bins)


colnames(cmv.count.frac.bins) <- cmv.count.df.bins$bin_2percent
cmv.count.frac.bins <- as.data.frame(cmv.count.frac.bins)
rownames(cmv.count.frac.bins) <- gsub("-\\(","_",rownames(cmv.count.frac.bins)  )
rownames(cmv.count.frac.bins) <- gsub("\\)","",rownames(cmv.count.frac.bins)  )
rownames(cmv.count.frac.bins) <- gsub("-","_",rownames(cmv.count.frac.bins)  )


### make heatmap
cmv.count.frac.bins.scale <- cmv.count.frac.bins
#d <- dist(cmv.count.frac.bins.scale)
#hc <- hclust(d, method = "average")
#cl <- cutree(hc, 7)


kmeans_res <- kmeans(cmv.count.frac.bins.scale, centers = 5, nstart = 25)
kmeans_order <- order(kmeans_res$cluster)
cluster <- kmeans_res$cluster


cmv.count.frac.bins.scale$gene <- rownames(cmv.count.frac.bins.scale)
cmv.count.frac.bins.scale$cluster <- cluster


cmv_order <- rownames(cmv.count.frac.bins.scale)[kmeans_order]
cmv_out <- factor(cmv_order, levels = cmv_order)


cmv.count.frac.bins.scale.lng <- tidyr::gather(cmv.count.frac.bins.scale, vload_percent_bin,  exp, -gene, -cluster)
cmv.count.frac.bins.scale.lng$gene <- factor(cmv.count.frac.bins.scale.lng$gene , levels = rev(cmv_order))
cmv.count.frac.bins.scale.lng$vload_percent_bin <- as.numeric(cmv.count.frac.bins.scale.lng$vload_percent_bin )
cmv.count.frac.bins.scale.lng$vload_percent_bin <- as.factor(cmv.count.frac.bins.scale.lng$vload_percent_bin)


library(viridis)
ggplot(cmv.count.frac.bins.scale.lng, aes(vload_percent_bin, gene, fill = exp)) +
  geom_tile() +
#  facet_grid(rows = cmv.count.frac.bins.scale.lng$cluster, scales = "free_y") +
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  scale_fill_viridis() +
  theme_classic()


ggsave("plots/heatmap_bin_cmv_exp_5_soupx.jpg", height = 13, width = 6)
ggsave("plots/heatmap_bin_cmv_exp_5_soupx.pdf", height = 13, width = 6)
ggsave("plots/heatmap_bin_cmv_exp_5_soupx.eps", height = 13, width = 6, device = "eps")






#### plot cluster in uninfected cells #####
ggplot(cmv.count.pct.df.scale.lng.tmp[which(cmv.count.pct.df.scale.lng.tmp$cmv_gene %in% g.plot),], aes(cell, cmv_gene, fill = exp)) +
  geom_tile() +
  facet_wrap(~Infection_localminima, scales = "free_x") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_viridis()


ggplot(cmv.count.pct.df.scale.lng.uninfect, aes(cell, cmv_gene, fill = exp)) +
  geom_tile() +
  facet_wrap(~Infection_localminima, scales = "free_x") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_viridis()


g.keep <- as.character(unique(cmv.count.frac.bins.scale.lng$gene))
g.orig <- as.character(unique(cmv.count.pct.df.scale.lng.uninfect$cmv_gene))


g.plot <- g.orig[which(g.orig %in% g.keep)]
g.orig[-which(g.orig %in% g.keep)]


cluster.g.orig <- cluster[g.plot]
cluster.g.orig<- as.data.frame(cluster.g.orig)
cluster.g.orig$cmv_gene <- rownames(cluster.g.orig)
colnames(cluster.g.orig)[1] <- "cluster"
cmv.count.pct.df.scale.lng.uninfect <- merge(cmv.count.pct.df.scale.lng.uninfect, cluster.g.orig, by = "cmv_gene")


ggplot(cmv.count.pct.df.scale.lng.uninfect, aes(cell, cmv_gene, fill = exp)) +
  geom_tile() +
  facet_grid(rows = cmv.count.pct.df.scale.lng.uninfect$cluster, scales = "free_y") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_viridis()


chisq.test(cmv.count.pct.df.scale.lng.uninfect[,c(7,10)])


######## plot smooth plot for each cluter  #######
library(ggbeeswarm)
library(RColorBrewer)


library(wesanderson)


head(cmv.count.frac.bins.scale.lng)
cmv.count.frac.bins.scale.lng$vload_percent_bin_num <- as.numeric(cmv.count.frac.bins.scale.lng$vload_percent_bin)




ggplot(cmv.count.frac.bins.scale.lng, aes(vload_percent_bin_num, exp)) +
  geom_jitter(aes(color= as.factor(cluster)), alpha = 0.5, size = 0.5)+
  geom_line(aes(group = gene, color = as.factor(cluster)),alpha = 0.5) + 
  geom_smooth()+
  facet_wrap(~cluster, ncol = 1) +
 # scale_color_manual(values=wes_palette( name="FantasticFox1")) +
  scale_color_brewer(palette = "Set1") +
  theme_classic()+
  ylim(0,1)


ggsave("plots/cluster_trajectory_cmv.pdf", height = 10, width = 7)
