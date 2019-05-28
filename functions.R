library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

# 1. Get color vectors
getColors <- function(n) {
  col <- brewer.pal.info[brewer.pal.info$category=='qual', ] # get max. 74 colours
  col_vector <- unlist(mapply(brewer.pal, col$maxcolors, rownames(col)))
  ifelse (n > length(col_vector), 
          vec <- sample(col_vector, n, replace=T),
          vec <- sample(col_vector, n, replace=F)
  )
  vec
}


# 2. Draw the heatmaps
draw_heatmap <- function(voomObj, topTable, phenoDF, list) {
  hm_cdr <- phenoDF %>% select(Sample.Type, tumor_stage)
  rownames(hm_cdr) <- colnames(voomObj)
  colnames(hm_cdr)  <- c('tumor', 'stage')
  tumor <- c("#99a599", "#37637f")
  names(tumor) <- unique(hm_cdr$tumor)
  stage <- c("#d7191c","#fdae61","#a1d99b","#2b83ba","#bababa")
  names(stage) <- c("I","II","III","IV","NA")
  anno_colors <- list(tumor = tumor,
                      stage = stage) # the name must be consistent
  h <- pheatmap(voomObj$E[list, ], annotation_col=hm_cdr, annotation_colors=anno_colors, 
                labels_row = list, show_colnames = F)
  h
}

# 3. dens.plot
dens.plot <- function(table, colVec, yrange) {
  d <- plot(density(table[, 1]), col=colVec[1], 
            lwd=2, las=2, ylim=yrange, main="", xlab="") +
    abline(v=0, lty=3) + title(xlab="expr values") +
    for (i in 2:ncol(table)) {
      den <- density(table[, i])
      lines(den$x, den$y, col=colVec[i], lwd=2)
    } 
  d
}

# 4. Function to draw the boxplot for a single gene
single.box <- function(v, phenoDF, id, tt){
  t_pdata <- phenoDF %>% select(Sample.ID, Sample.Type) 
  exp_list <- as.data.frame(v$E[rownames(v$E)==id, ])
  exp_list$Sample.ID <- rownames(exp_list)
  colnames(exp_list) <- c("counts", "Sample.ID")
  mdf <- merge(exp_list, t_pdata, by="Sample.ID")
  mdf$Sample.Type <- factor(mdf$Sample.Type)
  symbol <- id
  q_val <- tt[id, ]$adj.P.Val
  
  ggboxplot(mdf, x="Sample.Type", y="counts",
            color="Sample.Type", palette="jco", main=paste0(symbol, "  q-val = ", formatC(q_val, format="e", digits=2)),
            xlab="Tissue", ylab="logCPM", add="jitter", ggtheme = theme_bw())
  
}
