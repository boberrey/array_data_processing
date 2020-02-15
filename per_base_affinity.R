#!/usr/bin/env Rscript

## Merge aligned fragment bed file with single cluster
## fit information.

# Usage:
# make_fragment_bed.R <bam_file> <bed_file> <fasta_file>

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(readr)
  library(magrittr)
  library(Rsamtools)
  library(rtracklayer)
})

## Read in command line arguments
# args <- commandArgs(trailingOnly = FALSE)
# scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
# args <- args[(grep("^--args", args, value = FALSE)[1] + 1): length(args)]
# bamFile <- args[1]
# outFile <- args[2]
# fastaFile <- args[3]

# Debugging
setwd("/oak/stanford/groups/wjg/boberrey/array_data/HB27_aligning/20200118_TtAgoDM_55C/")
bedFile <- "output/beds/HB27_S17_L001.cluster.bed.gz"
fitFile <- "curve_fitting/constrained_fits.txt"
fastaIdx <- "/share/PI/wjg/lab/genomes/HB27/HB27.fa.fai"
outFile <- "binned_median_energies.tsv"
posBed <- "binned_median_energies_pos.bed"
negBed <- "binned_median_energies_neg.bed"
posBW <- "binned_median_energies_pos.bw"
negBW <- "binned_median_energies_neg.bw"

binWidth = 1
slideBy = 1

# Functions
dGtoKd <- function(dG, temp=37, R=0.0019872){
  return(1/exp(-dG/(R*(temp + 273))))
}

KdtodG <- function(Kd, temp=37, R=0.0019872){
  return(R*(temp + 273)*log(Kd))
}

dGtoKt <- function(dG, temp=37, R=0.0019872){
  return(dG / (R*(273 + temp)))
}

# Read in bedfile to GR
fragDF <- suppressMessages(data.frame(readr::read_tsv(
  bedFile, col_names = c("chr","start","end","cluster_ID","score","strand"))))

# Read in fit dataframe
fitDF <- suppressMessages(data.frame(readr::read_tsv(fitFile)))
colnames(fitDF)[1] <- "cluster_ID"

# Merge dataframes on cluster_ID
mergedDF <- base::merge(fragDF, fitDF, by="cluster_ID", sort=FALSE)

# Convert to genomic range
fitRange <- GenomicRanges::makeGRangesFromDataFrame(
  df = mergedDF, keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE)

# Remove failed fits
fitRange <- fitRange[fitRange$ier >= 1.0 & fitRange$ier <= 4.0]
# Filter on rmse / rsq
fitRange <- fitRange[which(fitRange$rmse < 0.3)]

# Threshold dG
concentrations <- c(0.055, 0.131, 0.314, 0.754, 1.808, 4.34, 10.42)
lowLimBound <- 0.80
upperLimBound <- 0.20
ttagoMinKd = round(concentrations[1] / lowLimBound - concentrations[1], 3)
ttagoMaxKd = round(concentrations[length(concentrations)] / upperLimBound - concentrations[length(concentrations)], 3)

fitRange$Kd <- pmax(fitRange$Kd, ttagoMinKd)
fitRange$Kd <- pmin(fitRange$Kd, ttagoMaxKd)
fitRange$dG_kT <- round(dGtoKt(KdtodG(fitRange$Kd / 1e9, temp = 55), temp = 55),3)

# Now let's see if we can get median dG...
# Prepare a overlapping tiled GR, then find overlaps, then calculate on hits?
# (see: https://support.bioconductor.org/p/100222/)

# First, get chr lengths
faidx <- data.frame(readr::read_tsv(fastaIdx, col_names = FALSE))
faidx <- faidx[,c(1,2)]
colnames(faidx) <- c("chr", "len")
chromSizes <- GRanges(faidx$chr, IRanges(1, faidx$len))
tiledRange <- unlist(tile(chromSizes, width = slideBy)) %>% resize(., width = binWidth)

# First, split fitRange by strand
fitRangePos <- fitRange[strand(fitRange) == '+']
fitRangeNeg <- fitRange[strand(fitRange) == '-']

# Find overlaps
hitsPos <- findOverlaps(tiledRange, fitRangePos, minoverlap = binWidth)
hitsNeg <- findOverlaps(tiledRange, fitRangeNeg, minoverlap = binWidth)
# Calculate median energy for all overlapping fragments
aggPos <- aggregate(fitRangePos, hitsPos, med_dG=median(dG_kT))
aggNeg <- aggregate(fitRangeNeg, hitsNeg, med_dG=median(dG_kT))
# Add median energies back to tiled range
tiledRange$posMedEnergy[countQueryHits(hitsPos) > 0L] <- round(aggPos$med_dG, 3)
tiledRange$negMedEnergy[countQueryHits(hitsNeg) > 0L] <- round(aggNeg$med_dG, 3)
# Add coverage cols as well
tiledRange$posCov <- countQueryHits(hitsPos)
tiledRange$negCov <- countQueryHits(hitsNeg)

## Write to tab data
out <- data.frame(
  chr = seqnames(tiledRange),
  start = as.integer(start(tiledRange)) - 1,
  end = as.integer(end(tiledRange)),
  posMed = tiledRange$posMedEnergy,
  posCov = tiledRange$posCov,
  negMed = tiledRange$negMedEnergy,
  negCov = tiledRange$negCov
) %>% readr::write_tsv(
  x = .,
  append = FALSE,
  path = outFile
)

## Write positive strand to bed file
out <- data.frame(
  chr = seqnames(tiledRange),
  start = as.integer(start(tiledRange)) - 1,
  end = as.integer(end(tiledRange)),
  name = paste0("bin", 1:length(tiledRange)),
  posMed = tiledRange$posMedEnergy,
  strand = rep('+', length(tiledRange)),
  posCov = tiledRange$posCov
) %>% readr::write_tsv(
  x = ., col_names = FALSE,
  append = FALSE,
  path = posBed
)


## Write negative strand to bed file
out <- data.frame(
  chr = seqnames(tiledRange),
  start = as.integer(start(tiledRange)) - 1,
  end = as.integer(end(tiledRange)),
  name = paste0("bin", 1:length(tiledRange)),
  negMed = tiledRange$negMedEnergy,
  strand = rep('-', length(tiledRange)),
  negCov = tiledRange$negCov
) %>% readr::write_tsv(
  x = ., col_names = FALSE,
  append = FALSE,
  path = negBed
)


# Adding seqlengths is required to export bw
seqlengths(tiledRange) <- GenomicRanges::width(chromSizes[seqnames(chromSizes) %in% seqlevels(tiledRange)])
posRange <- granges(tiledRange)
strand(posRange) <- "+"
posRange$score <- tiledRange$posMedEnergy
posRange <- posRange[!is.na(posRange$score)]

# Save track in bigwig format
export.bw(posRange, posBW)

negRange <- granges(tiledRange)
strand(negRange) <- "-"
negRange$score <- tiledRange$negMedEnergy
negRange <- negRange[!is.na(negRange$score)]

# Save track in bigwig format
export.bw(negRange, negBW)


