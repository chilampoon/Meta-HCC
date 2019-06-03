# Meta p-values for DEGs gotten from limma/limma-voom of each study

list.tt <- list()
load("~/hcc/data/TCGA/DE.TCGA.Rdata")
list.tt[['TCGA']] <- tt

load("~/hcc/data/ICGC/DE.ICGC.Rdata")
list.tt[['ICGC']] <- tt

load("~/hcc/data/GSE14520/DE_GSE14520_1.Rdata")
list.tt[['GSE14520_1']] <- ttable1

load("~/hcc/data/GSE14520/DE_GSE14520_2.Rdata")
list.tt[['GSE14520_2']] <- ttable2

load("~/hcc/data/GSE22058/DE_GSE22058.Rdata")
list.tt[['GSE22058']] <- ttable

load("~/hcc/data/GSE25097/DE_GSE25097.Rdata")
list.tt[['GSE25097']] <- ttable

# Intersect common genes
cases <- names(list.tt)
sapply(list.tt, nrow)
com.genes <- Reduce(intersect, list(rownames(list.tt[['TCGA']]), rownames(list.tt[['ICGC']]),
                                    rownames(list.tt[['GSE14520_1']]), rownames(list.tt[['GSE14520_2']]),
                                    rownames(list.tt[['GSE22058']]), rownames(list.tt[['GSE25097']])))

list.com <- lapply(list.tt, function(x) x[which(rownames(x) %in% com.genes),])

# Combine multiple p-values
library(metap)

pval.df <- setNames(data.frame(matrix(nrow = 0, ncol = length(list.tt)+2)), c('gene', names(list.tt), 'fisherp'))
for (g in com.genes) {
  ps <- c()
  for (df in list.tt) {
    ps <- c(ps, unlist(df[which(rownames(df)==g),]$P.Value))
  }
  pval.df[nrow(pval.df)+1,]$gene <- g
  pval.df[nrow(pval.df), 2:7] <- ps
  pval.df[nrow(pval.df),]$fisherp <-  sumlog(ps)$p
}

pval.df$adj.Pval <- p.adjust(pval.df$fisherp, method = 'BH', n = nrow(pval.df))
pval.df <- pval.df[order(pval.df$adj.Pval),]
rownames(pval.df) <- NULL
nrow(pval.df[pval.df$adj.Pval < 0.01,])
save(list.tt, pval.df, file='~/hcc/results/pval.Rdata')

