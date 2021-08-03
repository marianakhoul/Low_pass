setwd("~/Documents/Low_pass/")
library(biomaRt)
library(dplyr)
library(ggplot2)
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
ggplot(haSegments, aes(gene.symbol, copyNum)) +
  geom_point(aes(col = factor(ID)), size = 4) +
  geom_hline(yintercept = 2, col = "red", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual(values=c(getPalette2(3),
                              getPalette(3)))


#### regions ####


# Extract HGNC gene symbol
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
positions <- seg[1, 2:4]
positions <- data.frame(chr = 4, pos1 = 54244160, pos2 = 54292087)
results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = c("chromosome_name", "start", "end"), 
                 values = list(positions[,1], positions[,2], positions[,3]), mart = mart)
