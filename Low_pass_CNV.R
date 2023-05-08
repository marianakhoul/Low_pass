#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
# BiocManager::install("QDNAseq", force = T)
# BiocManager::install("QDNAseq.hg19", force = T)
# BiocManager::install("BSgenome", force = T)
# BiocManager::install("Biobase", force = T)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", force = T)

library(QDNAseq)
library(QDNAseq.hg19)
library(Biobase)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(future)
library(dplyr)
options(scipen=999)
setwd("F:/DFCI/Low_pass/")

exportSEG <- function(obj, fnames=NULL) {
  print("Overwriting segment calling function.")
  
  calls <- assayDataElement(obj, "calls")
  # segments <- log2adhoc(assayDataElement(obj, "segmented"))
  segments <- log2(assayDataElement(obj, "segmented"))
  
  fd <- fData(obj)
  pd <- pData(obj)
  
  if (is.null(fnames)) 
    fnames <- pd$name
  
  if (length(fnames) != length(pd$name)) {
    stop("Length of 'fnames' is too short: ", length(fnames), " != ", length(pd$name))
  }
  
  oopts2 <- options(scipen=100)
  on.exit(options(scipen=oopts2), add=TRUE)
  
  for (i in 1:ncol(calls)) {	
    d <- cbind(fd[,1:3],calls[,i], segments[,i])
    # sel <- d[,4] != 0 & !is.na(d[,4])
    sel <- !is.na(d[,4])
    
    dsel <- d[sel,]
    
    rleD <- rle(paste(d[sel,1], d[sel,4], sep=":"))
    
    endI <- cumsum(rleD$lengths)
    posI <- c(1, endI[-length(endI)] + 1)
    
    chr <- dsel[posI,1]
    pos <- dsel[posI,2]
    end <- dsel[endI,3]
    score <- dsel[posI,4]
    segVal <- round(dsel[posI,5],digits=2)
    bins <- rleD$lengths
    
    out <- cbind(fnames[i], chr, pos, end, bins, segVal)
    colnames(out) <- c("SAMPLE_NAME", "CHROMOSOME", "START", "STOP", "DATAPOINTS", "LOG2_RATIO_MEAN")
    
    fname <- paste(fnames[i], ".seg", sep="")
    
    write.table(out, fname, quote=FALSE, sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

assignInNamespace("exportSEG",exportSEG,ns="QDNAseq")

#### standard workflow ####

binsAndReads <- function(binSize, bamfile, threads){
  if(threads > 1){
    future::plan("multicore", workers=threads)
  }
  bins <- getBinAnnotations(binSize, genome = "hg19")
  readCounts <- binReadCounts(bins, pairedEnds = T, bamfiles=bamfile)
  if(threads > 1){
    future::plan("sequential")
  }
  return(readCounts)
}

filterAndBin <- function(readCounts){
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  return(list("readCounts" = readCounts, 
              "readCountsFiltered" = readCountsFiltered,
              "copyNumbersSmooth" = copyNumbersSmooth,
              "copyNumbersSegmented" = copyNumbersSegmented,
              "copyNumbersCalled" = copyNumbersCalled))
}

exportFiles <- function(QDNAseqObj, name) {
  exportBins(QDNAseqObj$copyNumbersSmooth,
             file=sprintf("copyNumbers.%s.igv", name),
             format="igv")
  QDNAseqObj$copyNumbersCalled@phenoData@data$name <- name
  exportBins(QDNAseqObj$copyNumbersCalled, 
             file = sprintf("%s.segments.igv", name), 
             format="igv", type = "segments")
  QDNAseqObj$copyNumbersCalled@phenoData@data$name <- 
    sprintf("%s.segments", name)
  exportBins(QDNAseqObj$copyNumbersCalled, 
             format = "seg", type="segments")
}

plotAll <- function(QDNAseqObj, title) {
  isobarPlot(QDNAseqObj$readCountsFiltered)
  noisePlot(QDNAseqObj$readCountsFiltered)
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=c(2,2))
  plot(QDNAseqObj$readCounts, logTransform=FALSE, 
       main = paste("Read counts,", title))
  highlightFilters(QDNAseqObj$readCounts, logTransform=FALSE, 
                   residual=TRUE, blacklist=TRUE)
  plot(QDNAseqObj$copyNumbersSmooth, 
       main = paste("Smoothed segments, ", title))
  plot(QDNAseqObj$copyNumbersSegmented, 
       main = paste("Segments,", title))
  plot(QDNAseqObj$copyNumbersCalled, 
       main = paste("Called CNVs,", title), xlab=NULL)
  mtext(QDNAseqObj$readCounts@phenoData@data$name,
        side=3, line=-20, outer = T)
  par(oldpar)
}

bamfile <- "X:/rs26/Files_For_Taylor/Samples_With_EGFR_AMP/KHC-14_SY-7661_S14_disambiguate.sorted.bam"

title <- "50_kbp_14"
readcounts50 <- binsAndReads(50, bamfile, 4) 
filteredReads50 <- filterAndBin(readcounts50)
exportFiles(filteredReads50, name = title)
plotAll(filteredReads50, title = title)

title <- "500_kbp_14"
readcounts500 <- binsAndReads(500, bamfile, 4) 
filteredReads500 <- filterAndBin(readcounts500)
exportFiles(filteredReads500, name = title)
plotAll(filteredReads500, title = title)

#### counting reads ####
egfrStart <- 55086725
egfrEnd   <- 55275031

segOverlap <- function(segmented, chrom, start, end) {
  # returns if the segs on chrom between start and end
  positions <- which(
    ((segmented$START >= start & # does the start fall within region?
        segmented$START <= end) |
        (segmented$STOP >= start & # or, does end fall within region?
           segmented$STOP <= end)) & # and, are those positions in chrom?
      segmented$CHROMOSOME == chrom)
  # if nothing is smaller than start + end, find the bin that captures
  # both
  if (length(positions) == 0) {
    positions <- which(segmented$START <= start & 
                         segmented$STOP >= end &
                         segmented$CHROMOSOME == chrom)
  }
  return(positions)
}

#### substitution ####
# NOTE: only works for singular segment replacements?
seg500 <- read.table("500_kbp_14.segments.seg", header = T)
seg50  <- read.table("50_kbp_14.segments.seg", header = T)

smallPos <- segOverlap(seg50, 7, egfrStart, egfrEnd)
largePos <- segOverlap(seg500, 7, egfrStart, egfrEnd)

smallSpan <- c(seg50$START[smallPos], seg50$STOP[smallPos])

insertFrame <- rbind(seg500[largePos,], seg500[largePos,])

insertFrame$STOP[1] <- (smallSpan[1]-1)
insertFrame$START[2] <- (smallSpan[2]+1)

segOut <- rbind(seg500[1:(largePos-1),],
      rbind(insertFrame[1,], seg50[smallPos,], insertFrame[2,]),
      seg500[(largePos+1):nrow(seg500),])
segOut[,1] <- rep("500+50kbp_hybrid", dim(segOut)[1])

write.table(segOut, file = "500_insert50.seg", quote = F, row.names = F,
            sep = "\t")

## sub igv bins ##
igvSegs500 <- read.table("copyNumbers.500_kbp_14.igv", header = T)
igvSegs50 <- read.table("copyNumbers.50_kbp_14.igv", header = T)

binOverlap <- function(bins, chrom, start, end) {
  # same as above but the names are different
  positions <- which(
    ((bins$start >= start & 
        bins$start <= end) |
       (bins$end >= start & 
          bins$end <= end)) & 
      bins$chromosome == chrom)
  if (length(positions) == 0) {
    positions <- which(bins$start <= start & 
                         bins$end >= end &
                         bins$chromosome == chrom)
  }
  return(positions)
}

largePos <- binOverlap(igvSegs500, 7, egfrStart, egfrEnd)
smallPos <- binOverlap(igvSegs50, 7, igvSegs500$start[largePos], 
                       igvSegs500$end[largePos])

outbins <- rbind(igvSegs500[1:(largePos-1),],
                 igvSegs50[smallPos,],
                 igvSegs500[(largePos+1):nrow(igvSegs500),])

write.table(outbins, file = "500_insert50.CN.igv", quote = F, 
            row.names = F, sep = "\t")





