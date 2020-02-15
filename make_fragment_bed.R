#!/usr/bin/env Rscript

## Convert a paired-end bam file into a fragment
## bed file while maintaining the read orientation

# Usage:
# make_fragment_bed.R <bam_file> <bed_file> <fasta_file>

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(readr)
})

## Read in command line arguments
args <- commandArgs(trailingOnly = FALSE)
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
args <- args[(grep("^--args", args, value = FALSE)[1] + 1): length(args)]
bamFile <- args[1]
outFile <- args[2]
fastaFile <- args[3]

# debugging
# setwd("/home/users/boberrey/oak/array_data/HB27_aligning/20200118_TtAgoDM_55C")
# scriptPath <- "/home/users/boberrey/git_clones/array_data_processing"
# bamFile <- "output/bams/filtered/HB27_S17_L001.filtered.bam"
# outFile <- "tst.bed"

# Functions
source(paste0(scriptPath, "/array_generics.R"))

# Get fragment GR
fragmentRange <- sort(bamToFragmentRange(bamFile), ignore.strand=TRUE)

# In order for bedtools getfasta to work, the strand needs to be in the 6th
# column. Make a fake score column to be column 5.
fragmentRange$score = "."

## Write to bed file
out <- data.frame(
  chr = seqnames(fragmentRange),
  start = as.integer(start(fragmentRange)) - 1,
  end = as.integer(end(fragmentRange)),
  names = names(fragmentRange),
  scores = fragmentRange$score,
  strands = strand(fragmentRange)
) %>% readr::write_tsv(
  x = .,
  append = TRUE,
  path = outFile
)

