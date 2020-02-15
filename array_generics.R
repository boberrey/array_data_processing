
## Generic functions for use in scATAC data processing and analysis
## Ben Ober-Reynolds
## boberrey@stanford.edu

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(magrittr)
  library(ggplot2)
})


#-------------
# IO Functions
#-------------


bamToFragmentRange <- function(bamFile, bamParam = ScanBamParam(mapqFilter = 20)){
  # Read in a paired-end bam file as a fragment genomic range (does NOT shift insertions!)
  # Convert bamFile to Rsamtools 'BamFile' class if not already done
  if(class(bamFile)[1] != "BamFile"){
    bamFile <- Rsamtools::BamFile(bamFile)
  }
  # Read in bam file as GenomicRanges
  # 'strandMode = 1' indicates that the 'strand of the pair is strand of its first alignment'
  fragments <- as(readGAlignmentPairs(
    bamFile, param = bamParam, use.names = TRUE,
    strandMode = 1), "GRanges")
  return(fragments)
}


#---------------------------
# Overlap analysis functions
#---------------------------


getInsertionEnrichment <- function(insertionRange, featureRange, 
                                   windowSize=2000, smoothWindow=51, normWindow=100){
  # Calculate insertion enrichment around a set of features
  overlaps <- findOverlapPairs(insertionRange, featureRange, 
                               ignore.strand = TRUE, maxgap = windowSize)
  # Correct for TSS directionality
  flips <- ifelse(strand(second(overlaps)) == "-", -1, 1)
  # Calculate strand-adjusted distance from insertion to feature center
  dists <- (start(first(overlaps)) - start(second(overlaps)))*flips
  distDF <- as.data.frame(table(dists))
  # Convert position from factor
  distDF$dists <- as.numeric(as.character(distDF$dists))
  # Normalize to mean of counts from positions +/- windowSize - normWindow
  normFactor <- mean(c(head(distDF$Freq, normWindow), tail(distDF$Freq, normWindow)))
  distDF$normed <- distDF$Freq / normFactor
  # Get smoothing line
  distDF$smoothed <- zoo::rollmean(distDF$Normed, smoothWindow, fill = 0)
  return(distDF)
}


nearestFeature <- function(x, subject, ignore.strand = TRUE, fix = "center"){
  # Find and report the single nearest feature in subject for each
  # entry in query. 'upstream' distances will be negative, and 'downstream'
  # distances will be positive.
  # Both x and subject are resized to have width 1 prior to finding distances
  x <- GenomicRanges::resize(x, width = 1, fix = fix)
  subject <- GenomicRanges::resize(subject, width = 1, fix = fix)
  
  nearest <- subject[GenomicRanges::nearest(x, subject, ignore.strand = ignore.strand, select = "arbitrary")]
  nearest$distance <- start(nearest) - start(x)
  return(nearest)
}


#---------------------------
# Miscellaneous functions
#---------------------------

loadBSgenome <- function(gnome){
  # Function to load appropriate BSgenome object
  if(gnome == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    g <- BSgenome.Hsapiens.UCSC.hg19
  }else if(gnome == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    g <- BSgenome.Hsapiens.UCSC.hg38
  }else if(gnome == "mm9"){
    library(BSgenome.Mmusculus.UCSC.mm9)
    g <- BSgenome.Mmusculus.UCSC.mm9
  }else if(gnome == "mm10"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    g <- BSgenome.Mmusculus.UCSC.mm10
  }else if(gnome == "dm6"){
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    g <- BSgenome.Dmelanogaster.UCSC.dm6
  }else{
    warning("Invalid genome selected. Please pick one of {hg19, hg38, mm9, mm10, dm6}")
    g <- NA
  }
  g
}


loadTxDb <- function(gnome){
  # Function to load appropriate TxDb object
  if(gnome == "hg19"){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    g <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(gnome == "hg38"){
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    g <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(gnome == "mm9"){
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    g <- TxDb.Mmusculus.UCSC.mm9.knownGene
  }else if(gnome == "mm10"){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    g <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(gnome == "dm6"){
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)  # No .knownGene for dm6?
    g <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  }else{
    warning("Invalid genome selected. Please pick one of {hg19, hg38, mm9, mm10, dm6}")
    g <- NA
  }
  g
}

loadOrgDb <- function(gnome){
  # Function to load appropriate TxDb object
  if(gnome == "hg19" || gnome == "hg38"){
    library(org.Hs.eg.db)
    g <- org.Hs.eg.db
  }else if(gnome == "mm9" || gnome == "mm10"){
    library(org.Mm.eg.db)
    g <- org.Mm.eg.db
  }else{
    warning("Invalid genome selected. Please pick one of {hg19, hg38, mm9, mm10}")
    g <- NA
  }
  g
}




baseContentFromGRange <- function(grange, bsgenome){
  ## Function to get the nucleotide frequencies from a GenomicRange object
  
  seqs <- Biostrings::getSeq(bsgenome, grange)
  baseCounts <- Biostrings::alphabetFrequency(seqs, as.prob = TRUE)
  return(baseCounts[,colSums(baseCounts != 0) > 0])
}


between <- function(x, low, high){return((x > low) & (x < high))}


takeSpread <- function(values, n){
  # Given an ordered vector of values, take n 'evenly' spaced values
  l <- length(values)
  selected <- lapply(1:n, function(x) values[ceiling(x*l/n)]) %>% unlist()
  return(selected)
}


# Function to get point density
getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

