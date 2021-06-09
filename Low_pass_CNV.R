#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("QDNAseq", force = T)
# BiocManager::install("QDNAseq.hg19", force = T)
# BiocManager::install("BSgenome")
# BiocManager::install("Biobase", force = T)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(QDNAseq)
library(QDNAseq.hg19)
library(Biobase)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(future)
library(dplyr)
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
    future::plan("multiprocess", workers=threads)
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

readOverlap <- function(readcounts, chrom, start, end) {
  # returns if the bins on chrom between start and end
  df <- readcounts@featureData@data
  positions <- which(
    ((df$start >= start & # does the start fall within region?
        df$start <= end) |
        (df$end >= start & # or, does end fall within region?
           df$end <= end)) & # and, are those positions in chrom?
      df$chromosome == chrom)
  # if nothing is smaller than start + end, find the bin that captures
  # both
  if (length(positions) == 0) {
    positions <- which(df$start <= start & 
                         df$end >= end &
                         df$chromosome == chrom)
  }
  return(positions)
}

readOverlap(readcounts50, 7, egfrStart, egfrEnd)

#### substitution ####
smallSegs <- filteredReads50$copyNumbersCalled
largeSegs <- filteredReads500$copyNumbersCalled

largePos <- readOverlap(largeSegs, 7, egfrStart, egfrEnd)
largeStart <- largeSegs@featureData@data$start[largePos]
largeEnd <- largeSegs@featureData@data$end[largePos]

smallPos <- readOverlap(smallSegs, 7, largeStart, largeEnd)





