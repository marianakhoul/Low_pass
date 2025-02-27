#!/usr/bin/env Rscript

# Authors: 
#Taylor Falk, taylora_falk@dfci.harvard.edu
#Maria Nakhoul, maria_nakhoul@dfci.harvard.edu

Sys.setenv("DISPLAY"=":0.0")
suppressPackageStartupMessages(require(optparse))
current_time <- format(Sys.time(), "%y-%m-%d_%H:%M")
validBins <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)
# generated with cpdm_sample_processing.R
# keep "ALL" as last element
defaultGenes <- 
  data.frame(gene=c("CCND2", "CDK4", "MET", "KDR", "KIT", "PDGFRA",
                    "RHEB", "XRCC2", "ALL"),
             chrom=c(12, 12, 7, 4, 4, 4, 7, 7, NA),
             start=c(4383459, 58143168, 116335903, 55946350, 55524310, 
                     55125049, 151164346, 151164346, NA),
             end=c(4409232, 58145511, 116436179, 55991522, 55604790, 
                   55161497, 151216655, 151216655, NA))
# add some extra space around the area
defaultGenes$start <- defaultGenes$start - (75*1000)
defaultGenes$end <- defaultGenes$end + (75*1000)
#### Set up command line arguments ####
option_list = list(
  make_option(c("-b", "--bamfile"), type="character", default=NULL, 
              help="bamfile location", metavar="file path"),
  make_option(c("-g", "--gene"), type="character", default=NULL,
              help=paste0("specify gene names instead of chrom, start, end",
                          " to use as region of interest.\n\t\t",
                          "Must be comma separated. Specify 'ALL' to select all available genes.\n\tCurrent options:\n\t",
                          paste0(defaultGenes$gene, collapse = ", "))),
  make_option(c("-c", "--chromosome"), type="integer", default=NULL, 
              help="chromosome where focal cnv is location", 
              metavar="integer"),
  make_option(c("-s", "--start"), type="integer", default=NULL, 
              help="starting position of focal cnv", metavar="coordinate"),
  make_option(c("-e", "--end"), type="integer", default=NULL, 
              help="ending position of focal cnv", metavar="coordinate"),
  make_option(c("-S", "--smallbins"), type="integer", default=50, 
              help=sprintf("size of the small bins to insert in between start and end, in kilobases.\n\t\tValid sizes: %s", 
                           paste(validBins, collapse = ", ")), 
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
  make_option(c("-w", "--warnings"), type="logical", default=TRUE, 
              help="Enable warnings. FALSE skips prompts that warn of negative behavior. Default=%default", 
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
if(!exists("bamfile", where=opt)){
  stop(sprintf("--bamfile is a mandatory field with no default. Please specify."))
}
if (!file.exists(opt$bamfile)) {
  stop(sprintf("Bam file %s not found.", opt$bamfile))
}
if(exists("gene", where=opt)){
  genes <- strsplit(opt$gene, ",")[[1]]
  if ("ALL" %in% genes) {
    genes <- defaultGenes$gene[1:length(defaultGenes$gene)-1]
  }
  # place coords into opt if gene is valid
  for (gene in genes) {
    if (!(gene %in% defaultGenes$gene)) {
      stop(sprintf("'%s' not in valid gene list. Valid genes:\n %s",
                   gene, paste0(defaultGenes$gene, collapse = ", ")))
    }
  }
} else {
  for (mandatoryOption in c("chromosome", "start", "end")) {
    if (!exists(mandatoryOption, where=opt)) {
      cat("--gene was not specified, chromosome, start, and end are mandatory.\n")
      stop(sprintf("'%s' is a mandatory field with no default. Please specify.",
                   mandatoryOption))
    }
  }
}

# extend positions
## since we want to compare our area of interest to the surrounding 
## cn distribution, we extend the bin replacement up and downstream by
## the length of the input start - stop

geneDT <- data.frame(gene = NULL, chrom = NULL, span = NULL, 
                     start = NULL, stop = NULL,  statStart = NULL, 
                     statStop = NULL)
if(exists("gene", where=opt)) {
  for (gene in genes) {
    geneRow <- defaultGenes[which(defaultGenes$gene == gene),]
    posSpan <- geneRow$end - geneRow$start
    tmp <- data.frame(gene = gene, span = posSpan,
                      chrom = geneRow$chrom,
                      start = geneRow$start, stop = geneRow$end, 
                      statStart = geneRow$start - posSpan, 
                      statStop = geneRow$end + posSpan)
    geneDT <- rbind(geneDT, tmp)
  }  
} else {
  posSpan <- opt$end - opt$start
  geneDT <- data.frame(gene = "ulp", chrom=opt$chromosome, span = posSpan,
                       start = opt$start, stop = opt$end, 
                       statStart = opt$start - posSpan, 
                       startStop = opt$end + posSpan)
}
# stop if span negative
if (any(geneDT$span < 1)) {
  stop(sprintf("The span of a region of interest cannot be negative. Span: %s",
               paste(geneDT$span, collapse = ", ")))
} 


# detail inputs
cat(sprintf(paste0("Combining QDNAseq outputs for %s.\n",
                   "\nLarge bin size: %i, Small bin size: %i\n",
                   "Running on %i thread(s), custom segmenting function: %s\n",
                   "Files saved with tag %s\n"),
            opt$bamfile, opt$largebins, opt$smallbins, opt$threads, 
            opt$zerosegments, opt$out))
for (row in 1:length(geneDT$gene)) {
  cat(sprintf(
    "For gene %s, replacing bins and segments in chromosome %i pos %i-%i\n",
    geneDT$gene[row], geneDT$chrom[row], geneDT$start[row], geneDT$stop[row]
  ))
}

# load noisy packages quietly
suppressPackageStartupMessages(require(Biobase))
suppressPackageStartupMessages(require(BSgenome))
suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(future))
suppressPackageStartupMessages(require(QDNAseq))
suppressPackageStartupMessages(require(QDNAseq.hg19))
suppressPackageStartupMessages(require(ggplot2))
options(scipen=999999999)

# load custom segments function if flagged
if (opt$zerosegments) {
  exportSEG <- function(obj, fnames=NULL) {
    cat("Overwriting segment calling function.\nd")
    
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
    future::plan("multicore", workers=threads)
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
    future::plan("multicore", workers=threads)
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
  exportBins(QDNAseqObj$copyNumbersCalled, format = "seg", file=names$seg,
             type="segments")
  names$seg <- sprintf("%s.SEG.seg", name)
  return(names)
}

plotAll <- function(QDNAseqObj, title) {
  # save a few plots as pngs
  png(filename = sprintf("Isobar_plot_%s.png", title), width = 720, type="cairo")
  isobarPlot(QDNAseqObj$readCountsFiltered)
  dev.off()
  png(filename = sprintf("Noise_plot_%s.png", title), type="cairo")
  noisePlot(QDNAseqObj$readCountsFiltered)
  dev.off()
  png(filename = sprintf("Bins-segments_plot_%s.png", title),
      width = 960, height = 960, units = "px", type="cairo") 
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
cat("\n##### Creating small bins and segments #####\n")
smtitle <- sprintf("%ikbp_%s", opt$smallbins, opt$out)
smallRC <- binsAndReads(opt$smallbins, opt$bamfile, opt$threads) 
filteredSmallRC <- filterAndBin(smallRC, opt$threads)
smallNames <- exportFiles(filteredSmallRC, name = smtitle)
if (opt$saveplots) {
  plotAll(filteredSmallRC, title = smtitle)
}

# then create the large bin files
cat("\n##### Creating large bins and segments #####\n")
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
    names(segmented)[1] <- "CHROMOSOME"
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

#### combine segments ####
if (F) {
  largeNames <- list(seg = "testing/500kbp_t3.SEG.seg")
  smallNames <- list(seg = "testing/50kbp_t3.SEG.seg")
  opt <- list(chromosome = 2, start = 16078672, end = 16089125)
}

combineSegs <- function(gene, chrom, start, stop) {
  largeSegs <- read.table(largeNames$seg, header = T)
  smallSegs <- read.table(smallNames$seg, header = T)
  
  # record the positions that overlap the region of interest
  largePos <- findOverlap(largeSegs, chrom, start, stop)
  smallPos <- findOverlap(smallSegs, chrom, start, stop)
  
  # find the start and stop of the small segment 
  smallSpan <- c(min(smallSegs$START[smallPos]), max(smallSegs$STOP[smallPos]))
  # the data.frame we will modify and reinsert
  insertFrame <- rbind(largeSegs[largePos,], largeSegs[largePos,])
  insertFrame$STOP[1] <- (smallSpan[1]-1) # truncate large segment lengths
  insertFrame$START[2] <- (smallSpan[2]+1)
  
  # create the full frame with the new segments for large and small
  segOut <- rbind(largeSegs[1:(largePos-1),],
                  rbind(insertFrame[1,], smallSegs[smallPos,], insertFrame[2,]),
                  largeSegs[(largePos+1):nrow(largeSegs),])
  segOut$LOG2_RATIO_MEAN > 1.0
  segOut[,1] <- rep(sprintf("%s_%ikbp_insert%ikbp_%s", gene,
                            opt$largebins, opt$smallbins, opt$out),
                    dim(segOut)[1]) # change the name for igv
  
  write.table(segOut, 
              file = sprintf("%s_%ikbp_insert%ikbp_%s.seg", gene,
                             opt$largebins, opt$smallbins, opt$out), 
              quote = F, row.names = F, sep = "\t")
}

for (row in 1:length(geneDT$gene)) {
  combineSegs(geneDT$gene[row], geneDT$chrom[row], geneDT$start[row],
              geneDT$stop[row])
}

#### sub igv bins ####

smallOverlaps <- function(binTup, posTup) {
  # determines the percentage of overlap for a bin with position
  # "tuple"  means just c(1, 1), i know, i know
  # remove bins that don't have 10%? covered by the position of interest
  if ((binTup[1] < posTup[1] & binTup[2] <= posTup[1]) | 
      (binTup[1] >= posTup[2] & binTup[2] > posTup[2])) {
    return(0)
  } else if (binTup[1] < posTup[1] & binTup[2] >= posTup[1]) {
    return((binTup[2] - posTup[1])/(binTup[2] - binTup[1]))
  } else if (binTup[1] < posTup[2] & binTup[2] >= posTup[2]) {
    return((posTup[2] - binTup[2])/(binTup[2] - binTup[1]))
  } else if (binTup[1] >= posTup[1] & binTup[2] <= posTup[2]){
    return((binTup[2] - binTup[1])/(binTup[2] - binTup[1]))
  }
}

# statistics
runStats <- function(smallCN, smallPos, gene, chrom, start, stop){
  labels <- c(rep("nontarget", length(which(smallCN$end[smallPos] <= start))),
              rep("target", length(which(smallCN$end[smallPos] > start & 
                                           smallCN$start[smallPos] < stop))),
              rep("nontarget", length(which(smallCN$start[smallPos] >= stop))))
  
  stats <- data.frame(copynum = smallCN[smallPos,5], 
                      position = factor(labels, levels = c("nontarget", "target")))
  
  # perform and record some distribution comparisons
  cat(paste0("\nComparing copy number of bins in targeted region and the ",
             "two adjacent upstream and downstream areas of equal length."))
  
  # compare the distributions to one another
  cat(sprintf(paste0("\nNum target bins: %i non-target bins: %i\n"),
              length(which(stats$position == "target")),
              length(which(stats$position == "nontarget"))))
  eqVar <- var.test(stats$copynum ~ stats$position)$p.value >= 0.05
  ttest <- t.test(stats$copynum ~ stats$position, var.equal = eqVar)
  cat(sprintf(paste0("\nComparing the distribtions (target vs non-target).\n",
                     "Using %s"), ttest$method))
  # print lots of info
  if (ttest$p.value < 0.05) {
    cat(sprintf(paste0("\nA significant difference between the target and non-target",
                       " distibutions was detected (p-value: %f). Determining ",
                       "which region is greater or lesser..."),
                ttest$p.value))
    # target GREATER than nontarget
    greaterT <- t.test(stats$copynum[stats$position == "target"], 
                       mu=mean(stats$copynum[stats$position == "nontarget"]), 
                       alternative="greater", var.equal = eqVar)
    # target LESS than nontarget
    lesserT <- t.test(stats$copynum[stats$position == "target"], 
                      mu=mean(stats$copynum[stats$position == "nontarget"]), 
                      alternative="less", var.equal = eqVar)
    if (greaterT$p.value < 0.05) {
      wordChoice <- "greater than"
      pval <- greaterT$p.value
    } else {
      wordChoice <- "less than"
      pval <- lesserT$p.value
    }
    cat(sprintf(paste0("\nThe distribtion of copy numbers across the bins ",
                       "in the target area were found to be %s the non-target ",
                       "area. \np-value: %f, target bins: %i non-target bins: %i"),
                wordChoice, pval, length(which(stats$position == "target")),
                length(which(stats$position == "nontarget"))))
    cat(sprintf("\nTarget copy num mean: %f Nontarget mean: %f\n",
                mean(stats$copynum[stats$position == "target"]),
                mean(stats$copynum[stats$position == "nontarget"])))
  } else {
    cat(sprintf(paste0("\nNo signficant difference found in the disitribtion ",
                       "of copy numbers in bins between target and non-target ",
                       "regions.\n p-value: %f, target bins: %i non-target bins: %i"),
                ttest$p.value, length(which(stats$position == "target")),
                length(which(stats$position == "nontarget"))))
    cat(sprintf("\nTarget copy num mean: %f Nontarget mean: %f\n",
                mean(stats$copynum[stats$position == "target"]),
                mean(stats$copynum[stats$position == "nontarget"])))
  }
  
  if (opt$saveplots) {
    # plot the histograms
    title <- sprintf("%s_%ikbp_insert%ikbp_%s", gene, opt$largebins, 
                     opt$smallbins, opt$out)
    title <- sprintf("hist_bin_distribution_%s.png", title)
    if (length(title) < 1) {
      title <- "position_distribution.png"
    }
    png(filename = title, width = 960, height = 960, units = "px", type="cairo") 
    p <- ggplot(stats, aes(x = copynum, fill = position)) + 
      geom_histogram(bins = 60) +
      scale_fill_brewer(palette="Dark2") +
      theme_light() +
      theme(legend.position = "none") +
      xlab("Copy number of bin") +
      ylab("Number of bins") +
      facet_wrap(~ position, ncol = 1)
    print(p)
    dev.off()
  } 
}

combineBins <- function(gene, chrom, start, stop){
  largeCN <- read.table(largeNames$copynum, header = T)
  smallCN <- read.table(smallNames$copynum, header = T)
  
  largePos <- findOverlap(largeCN, chrom, start, stop)
  
  for (i in largePos) {
    percentCovered <- smallOverlaps(c(largeCN$start[i], largeCN$end[i]),
                                    c(start, stop))
    if (abs(percentCovered) < 0.03) {
      largePos <- largePos[-which(largePos == i)]
    }
  }
  
  if (length(largePos) == 0) {
    stop(sprintf(paste0("ERROR: No bins were found in the copy number ",
                        "file %s at the area specified, are you sure ",
                        "you have the right chromosome and positions?"), 
                 largeNames$copynum))
  }
  
  smallPos <- findOverlap(smallCN, opt$chromosome, 
                          min(largeCN$start[largePos]),
                          max(largeCN$end[largePos]))
  
  # warn if there are empty bins in the target area
  # if (opt$warnings) {
  #   for (i in c(2:length(smallPos))) {
  #     print(smallCN[smallPos,])
  #     if (smallCN[smallPos,][i,2] - smallCN[smallPos,][i-1,3] > 1){
  #       response <- readline(prompt=paste0("A bin with no reads was found ",
  #                                          "in the data. Try using larger ",
  #                                          "bin sizes. Continue? (Y/n): "))
  #       if (grepl("Y", response, ignore.case = T) | response == "") {
  #         break
  #       } else {
  #         stop("Bin found with no reads present, halting program...")
  #       }
  #     }
  #   }
  # }
  
  # warn if any bins are outside one SD
  mean <- mean(smallCN[smallPos,5])
  sd <- sqrt(var(smallCN[smallPos,5]))
  numOutsideSD <- length(which(smallCN[smallPos,5] < (mean - 1.5*sd)))
  if (numOutsideSD > 0){
    if (opt$warnings){
      cat(sprintf("\n%i bin(s) were found to be outside 1.5 SD.\n", numOutsideSD))
      response <- readline(prompt=paste0("\nThis could indicate the bin size ",
                                         "is too small. Consider using larger ",
                                         "bins. Continue? (Y/n): "))
      if (grepl("Y", response, ignore.case = T) | response == "") {
        next
      } else {
        stop("Bin found with low reads, halting program...")
      }
    }
  }
  
  # reassign large/small bins to not the stats areas
  largePos <- findOverlap(largeCN, chrom, start, stop)
  smallPos <- findOverlap(smallCN, chrom, 
                          min(largeCN$start[largePos]),
                          max(largeCN$end[largePos]))
  
  outbins <- rbind(largeCN[1:(min(largePos)-1),],
                   smallCN[smallPos,],
                   largeCN[(max(largePos)+1):nrow(largeCN),])
  
  igvOutFile <- sprintf("%s_%ikbp_insert%ikbp_CN_%s.igv", gene,
                        opt$largebins, opt$smallbins, opt$out)
  write.table(outbins, file = igvOutFile, 
              quote = F, row.names = F, sep = "\t")
  
  igvTemp <- readLines(igvOutFile)
  write(c("#type=COPY_NUMBER", "#track coords=1", igvTemp),
        file = igvOutFile)
  runStats(smallCN, smallPos, gene, chrom, start, stop)
}

for (row in 1:length(geneDT$gene)) {
  combineBins(geneDT$gene[row], geneDT$chrom[row], geneDT$start[row],
              geneDT$stop[row])
}
