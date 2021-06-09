# file: QDNAseq.R

#library(QDNAseq, lib="/ifs/rcgroups/dxbx-cb/tools/miniconda/envs/ichorCNA/lib/R/library/")
library(QDNAseq, lib="/mnt/centers/bxdx/dxbx-cb/tools/miniconda/envs/ichorCNA/lib/R/library/")
setwd("F:/DFCI/Low_pass/")

#library(QDNAseq)

# 1. Get Bin Annotations

#bins <- readRDS("/ifs/rcgroups/dxbx-cb/shared/panels/QDNAseq/QDNAseq.hg19.100kbp.SR50.rds")
#bins <- readRDS("/mnt/centers/bxdx/dxbx-cb/shared/panels/QDNAseq/QDNAseq.hg19.100kbp.SR50.rds")
bins <- getBinAnnotations(binSize=100, genome="hg19")

#setwd(/directory with bams)

readCounts <- binReadCounts(bins, pairedEnds = T)


# 2. Calculate and plot Raw Counts of all .bam files in working directory. Plot a raw copy number profile (read counts across the genome), and highlight bins that will be removed with default filtering

#readCounts <- binReadCounts(bins,isPaired=TRUE)
exportBins(readCounts, file="raw_readCounts.txt")

pdf("Raw_read_counts.pdf")

plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))

highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)

dev.off()

# 3. Apply filters to the read counts and plot median read counts as a function of GC content and mappability

readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

pdf("isobarPlot.pdf")

isobarPlot(readCountsFiltered)

dev.off()

readCountsFiltered <- estimateCorrection(readCountsFiltered)

pdf("noisePlot.pdf")

noisePlot(readCountsFiltered)

dev.off()

# 4. Correct bin counts for GC content and mappability

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

pdf("Norm_CN_plots.pdf")

plot(copyNumbersSmooth)

dev.off()

exportBins(copyNumbersSmooth, file="copyNumbers.100kb.txt")
exportBins(copyNumbersSmooth, file="copyNumbers.100kb.igv", format="igv")

# 5. Segmentation with the CBS algorithm from DNAcopy, and calling copy number aberrations with CGHcall or cutoffs


copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)


pdf("Segment_plots_log2.pdf",width=22*0.8,height=8*0.8)

plot(copyNumbersSegmented)

dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)

pdf("Segment_plots_v2.pdf",width=22*0.8,height=8*0.8)

plot(copyNumbersCalled)

dev.off()

save.image(file="QDNASeq_session.RData")
exportBins(copyNumbersCalled, file="segments_called.100kb.txt", format="tsv", type = "segments")
exportBins(copyNumbersCalled, format = "seg")

save.image(file="QDNASeq_session.RData")
