setwd("~/Documents/Low_pass/")
library(biomaRt)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

dt <- read.table("CNV_ONCDRS_output_data_Taylor_Data.txt", sep = "\t",
                 header = T)
seg <- read.table("Tyalor_CPDM_request_Oncopanel.seg", header = T)
cr <- read.table("Tyalor_CPDM_request_Oncopanel.log2ratio", header = T)
cr$gene_symbol <- sapply(strsplit(cr$gene_exon, "_"), "[", 1)

#### look at all genes ####
top100names <- names(sort(table(cr$gene_symbol), decreasing = T))[1:100]
top100 <- cr[which(cr$gene_symbol %in% top100names),]

meanTop100 <- top100 %>% 
  group_by(Sample, gene_symbol) %>% 
  summarise(mean = mean(log2.ratio))

ggplot(cr[which(cr$gene_symbol %in% top100names),], 
       aes(gene_symbol, log2.ratio)) +
  geom_point(aes(col = Sample)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Log2 Ratio for top 100 most common genes in samples")

ggplot(meanTop100, aes(gene_symbol, mean)) +
  geom_point(aes(col = Sample)) +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Mean Log2 Ratio for top 100 most common genes")

# HA genes to admire
#haGenes <- unique(dt[which(dt$CNV_TYPE_CD == "HA"),2])
haGenes <- unique(dt[which(dt$CNV_TYPE_CD == "HA"),1:2])

#### amplifications ####
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_gene_ensembl")
haTable <- getBM(attributes = c("hgnc_symbol","chromosome_name", 
                                "start_position", "end_position"), 
                 filters = c("hgnc_symbol"), 
                 values = list(haGenes[,2]), mart = mart)
haMerge <- merge.data.frame(haGenes, haTable, 
                            by.x = "GENE", by.y = "hgnc_symbol")
haSegments <- data.frame(ID=NULL, chrom=NULL, start=NULL, end=NULL,
                         seg.mean=NULL, gene.symbol=NULL)
for (row in 1:dim(haMerge)[1]) {
  rawSegment <- seg[which(seg$chrom == haMerge[row,3] & 
                            seg$loc.start <= haMerge[row,4] & 
                            seg$loc.end >= haMerge[row,5] &
                            seg$ID == haMerge[row, 2]),]
  df <- cbind(rawSegment, 
              gene.symbol = rep(haMerge[row,1], dim(rawSegment)[1]))
  haSegments <- rbind(haSegments, df)
}

haSegments$copyNum <- sapply(haSegments$seg.mean, function(x) {return((2^x)*2)})

# plot segment values
getPalette <- colorRampPalette(brewer.pal(6, "Set3"))
getPalette2 <- colorRampPalette(brewer.pal(6, "Set1"))
# png("amplified.png", width = 1024, height = 768, units = "px")
ggplot(haSegments, aes(gene.symbol, copyNum)) +
  geom_point(aes(col = factor(ID)), size = 4) +
  geom_hline(yintercept = 2, col = "red", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual(values=c(getPalette2(3),
                              getPalette(3)))
# dev.off()
#### regions ####
seg[which(seg$chrom == haMerge[row,3] &
            seg$loc.start <= haMerge[row,4] &
            seg$loc.end >= haMerge[row,5]),]

findOverlap <- function(segmented, chrom, start, end, margin) {
  start <- start/margin
  end <- end*margin
  # contain
  positions <- which(segmented$loc.start >= start &
        segmented$loc.end <= end &
        segmented$chrom == chrom &
        segmented$seg.mean >= 0.5)
  return(positions)
}


plotRegions <- function(row) {
  # tempPlot <- seg[which(seg$chrom == haMerge[row,3] &
  #                         seg$loc.start <= haMerge[row,4] &
  #                         seg$loc.end >= haMerge[row,5]),]
  tempPlot <- seg[findOverlap(seg, haMerge[row,3], 
                              haMerge[row,4], haMerge[row,5],
                              1.05),]
  p <- ggplot(tempPlot, aes(x = factor(ID), y = loc.start, 
                       xend = factor(ID), yend = loc.end)) + 
    geom_segment(aes(col = seg.mean), size = 4) + 
    geom_segment(aes(x = factor(haMerge[row,1]), y = haMerge[row,4],
                     xend = factor(haMerge[row,1]), yend=haMerge[row,5]),
                 size = 4, col = "red") +
    # geom_segment(aes(x = factor("mean"), y=mean(loc.start),
    #                  xend = factor("mean"), yend=mean(loc.end)),
    #              size = 4, col = "green") +
    xlab("CPDM ID") + ylab ("Position") +
    theme_minimal() +
    ggtitle(haMerge[row,1]) +
    coord_flip()
  return(p)
}

# png("avail_regions.png", width = 2048, height = 768, units = "px")
# grid.arrange(plotRegions(1), plotRegions(2), plotRegions(5), plotRegions(17),
#              plotRegions(29), plotRegions(32), plotRegions(36), plotRegions(40),
#              plotRegions(44), plotRegions(47), plotRegions(52), plotRegions(60),
#              ncol = 4, nrow = 3)
grid.arrange(plotRegions(1), plotRegions(2), plotRegions(5),
             plotRegions(10), plotRegions(13), plotRegions(17),
             plotRegions(22), plotRegions(24), ncol = 4, nrow = 2)
# dev.off()
# png("detail.png", width = 1024, height = 768, units = "px")
plotRegions(1)
# dev.off()

#### save regions ####
regions <- data.frame(ID=NULL, chrom=NULL, loc.start=NULL, loc.end=NULL,
                      seg.mean=NULL, gene=NULL)
for (i in c(1, 2, 13)) {
  df <- seg[findOverlap(seg, haMerge[i,3], haMerge[i,4], haMerge[i,5],
                                   1.05),]
  df$gene <- haMerge[i, 1]
  regions <- rbind(regions, df[which(df$seg.mean == max(df$seg.mean)),])
}
for (i in c(5, 10)) {
  df <- seg[findOverlap(seg, haMerge[i,3], haMerge[i,4], haMerge[i,5],
                        1),]
  df$gene <- haMerge[i, 1]
  regions <- rbind(regions, df[which(df$ID == "CPDM_1883X"),])
}

df <- seg[findOverlap(seg, haMerge[17,3], haMerge[17,4], haMerge[17,5],1),]
df$gene <- haMerge[17, 1]
regions <- rbind(regions, df[which(df$ID == "CPDM_1717X"),])

df <- seg[findOverlap(seg, haMerge[22,3], haMerge[22,4], haMerge[22,5],1),]
df$gene <- haMerge[22, 1]
regions <- rbind(regions, df)

df <- seg[findOverlap(seg, haMerge[24,3], haMerge[24,4], haMerge[24,5],1.01),]
df$gene <- haMerge[24, 1]
regions <- rbind(regions, df[which(df$ID == "CPDM_0890X"),])

write.csv(x = regions, file = "gene_regions")
