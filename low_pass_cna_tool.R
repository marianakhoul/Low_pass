#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
current_time <- format(Sys.time(), "%y-%m-%d_%H:%M")
validBins <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)
#### Set up command line arguments ####
option_list = list(
  make_option(c("-b", "--bamfile"), type="character", default=NULL, 
              help="bamfile location", metavar="file path"),
  make_option(c("-c", "--chromosome"), type="integer", default=NULL, 
              help="chromosome where focal cnv is location", 
              metavar="integer"),
  make_option(c("-s", "--start"), type="integer", default=NULL, 
              help="starting position of focal cnv", metavar="coordinate"),
  make_option(c("-e", "--end"), type="integer", default=NULL, 
              help="ending position of focal cnv", metavar="coordinate"),
  make_option(c("-S", "--smallbins"), type="integer", default=50, 
              help=sprintf("size of the small bins to insert in between start and end, in kilobases.\n\t\tValid sizes: %s", paste(validBins, collapse = ", ")), 
              metavar="integer"),
  make_option(c("-L", "--largebins"), type="integer", default=500, 
              help=sprintf("size of the larger bins that cover the remainder of the genome, in kilobase pairs.\n\t\tMUST be a multiple of small bin size.\n\t\tValid sizes: %s", paste(validBins, collapse = ", ")), 
              metavar="integer"),
  make_option(c("-t", "--threads"), type="integer", default=1, 
              help="determine how many threads QDNAseq will use. Default=%default", 
              metavar="integer"),
  make_option(c("-z", "--zerosegments"), type="logical", default=TRUE, 
              help="use custom segment combining function to return more segments. Default=%default", 
              metavar="logical"),
  make_option(c("-p", "--saveplots"), type="logical", default=FALSE, 
              help="Save plots of QDNAseq outputs using --out as labels. Default=%default", 
              metavar="logical"),
  make_option(c("-o", "--out"), type="character", 
              default=sprintf("ulpCNA_%s", current_time), 
              help="output file name. Default=\"%default\"", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

# catch some issues #
if (opt$largebins %% opt$smallbins != 0) {
  stop("Large bin size (-L) is not a multiple of small bin size (-S).")
}
if (!opt$smallbins %in% validBins) {
  stop("%i is not a valid bin size for small bins (-S). Valid bins are: ", opt$smallbins, 
       paste(validBins, collapse = ", "))
}
if (!opt$largebins %in% validBins) {
  stop("%i is not a valid bin size for large bins (-L). Valid bins are: ", opt$smallbins, 
       paste(validBins, collapse = ", "))
}
if (!file.exists(opt$bamfile)) {
  stop(sprintf("Bam file %s not found.", opt$bamfile))
}
for (mandatoryOption in c("bamfile", "chromosome", "start", "end")) {
  if (!exists(mandatoryOption, where=opt)) {
    stop(sprintf("'%s' is a mandatory field with no default. Please specify.",
                 mandatoryOption))
  }
}

# load noisy packages quietly
suppressPackageStartupMessages(require(Biobase))
suppressPackageStartupMessages(require(BSgenome))
suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(future))
suppressPackageStartupMessages(require(QDNAseq))
suppressPackageStartupMessages(require(QDNAseq.hg19))
options(scipen=999999999)

# load custom segments function if flagged
if (opt$zerosegments) {
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
}

#### run QDNAseq two times ####

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

filterAndBin <- function(readCounts, threads){
  # perform the main QDNAseq functions
  if(threads > 1){
    future::plan("multiprocess", workers=threads)
  }
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  if(threads > 1){
    future::plan("sequential")
  }
  return(list("readCounts" = readCounts, 
              "readCountsFiltered" = readCountsFiltered,
              "copyNumbersSmooth" = copyNumbersSmooth,
              "copyNumbersSegmented" = copyNumbersSegmented,
              "copyNumbersCalled" = copyNumbersCalled))
}

exportFiles <- function(QDNAseqObj, name) {
  # write the two files for bin copy num and seg and return names
  names <- list(copynum = sprintf("%s.CN.igv", name),
                seg     = sprintf("%s.SEG", name))
  exportBins(QDNAseqObj$copyNumbersSmooth,
             file=names$copynum,
             format="igv")
  QDNAseqObj$copyNumbersCalled@phenoData@data$name <- names$seg
  exportBins(QDNAseqObj$copyNumbersCalled, format = "seg", 
             type="segments")
  names$seg <- sprintf("%s.SEG.seg", name)
  return(names)
}

plotAll <- function(QDNAseqObj, title) {
  # save a few plots as pngs
  png(filename = sprintf("Isobar_plot_%s.png", title), width = 720)
  isobarPlot(QDNAseqObj$readCountsFiltered)
  dev.off()
  png(filename = sprintf("Noise_plot_%s.png", title))
  noisePlot(QDNAseqObj$readCountsFiltered)
  dev.off()
  png(filename = sprintf("Bins-segments_plot_%s.png", title),
      width = 960, height = 960, units = "px") 
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=c(2,2))
  plot(QDNAseqObj$readCounts, logTransform=FALSE, 
       main = paste("Read counts,", title))
  highlightFilters(QDNAseqObj$readCounts, logTransform=FALSE, 
                   residual=TRUE, blacklist=TRUE)
  plot(QDNAseqObj$copyNumbersSmooth, 
       main = paste("Smoothed segments,", title))
  plot(QDNAseqObj$copyNumbersSegmented, 
       main = paste("Segments,", title))
  plot(QDNAseqObj$copyNumbersCalled, 
       main = paste("Called CNVs,", title), xlab=NULL)
  mtext(QDNAseqObj$readCounts@phenoData@data$name,
        side=3, line=-40, outer = T)
  par(oldpar)
  dev.off()
}

# small bins first, takes longer
smtitle <- sprintf("%ikbp_%s", opt$smallbins, opt$out)
smallRC <- binsAndReads(opt$smallbins, opt$bamfile, opt$threads) 
filteredSmallRC <- filterAndBin(smallRC, opt$threads)
smallNames <- exportFiles(filteredSmallRC, name = smtitle)
if (opt$saveplots) {
  plotAll(filteredSmallRC, title = smtitle)
}
print(smallNames)

# then create the large bin files
lgtitle <- sprintf("%ikbp_%s", opt$largebins, opt$out)
largeRC <- binsAndReads(opt$largebins, opt$bamfile, opt$threads) 
filteredLargeRC <- filterAndBin(largeRC, opt$threads)
largeNames <- exportFiles(filteredLargeRC, name = lgtitle)
if (opt$saveplots) {
  plotAll(filteredLargeRC, title = lgtitle)
}

#### combine the bins and segs into new files ####
findOverlap <- function(segmented, chrom, start, end) {
  if (length(names(segmented)) == 5) { # change col headers to match
    names(segmented)[2:3] <- c("START", "STOP")
  }
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

# combine segments
# largeNames <- list(seg = "testing/1000kbp_test1.SEG.seg")
# smallNames <- list(seg = "testing/50kbp_test1.SEG.seg")
# opt <- list(chromosome = 7, start = 55086725, end = 55275031)
#
largeSegs <- read.table(largeNames$seg, header = T)
smallSegs <- read.table(smallNames$seg, header = T)

# record the positions that overlap the region of interest
largePos <- findOverlap(largeSegs, opt$chromosome, opt$start, opt$end)
smallPos <- findOverlap(smallSegs, opt$chromosome, opt$start, opt$end)

# find the start and stop of the small segment 
smallSpan <- c(smallSegs$START[smallPos], smallSegs$STOP[smallPos])
# the data.frame we will modify and reinsert
insertFrame <- rbind(largeSegs[largePos,], largeSegs[largePos,])
insertFrame$STOP[1] <- (smallSpan[1]-1) # truncate large segment lengths
insertFrame$START[2] <- (smallSpan[2]+1)

# create the full frame with the new segments for large and small
segOut <- rbind(largeSegs[1:(largePos-1),],
                rbind(insertFrame[1,], smallSegs[smallPos,], insertFrame[2,]),
                largeSegs[(largePos+1):nrow(largeSegs),])
segOut[,1] <- rep(sprintf("%ikbp_insert%ikbp_%s", opt$largebins, 
                          opt$smallbins, opt$out),
                  dim(segOut)[1]) # change the name for igv

write.table(segOut, 
            file = sprintf("%ikbp_insert%ikbp_%s.seg", opt$largebins, 
                           opt$smallbins, opt$out), 
            quote = F, row.names = F, sep = "\t")

## sub igv bins ##
# igvSegs500 <- read.table("copyNumbers.500_kbp_14.igv", header = T)
# igvSegs50 <- read.table("copyNumbers.50_kbp_14.igv", header = T)
# 
# 
# largePos <- binOverlap(igvSegs500, 7, egfrStart, egfrEnd)
# smallPos <- binOverlap(igvSegs50, 7, igvSegs500$start[largePos], 
#                        igvSegs500$end[largePos])
# 
# outbins <- rbind(igvSegs500[1:(largePos-1),],
#                  igvSegs50[smallPos,],
#                  igvSegs500[(largePos+1):nrow(igvSegs500),])
# 
# write.table(outbins, file = "500_insert50.CN.igv", quote = F, 
#             row.names = F, sep = "\t")
# 
