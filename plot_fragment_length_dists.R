#!/usr/bin/env Rscript

## Plot the fragment length distributions for each sample individually
## and then all samples together

# Usage:
# plot_fragment_length_dists.R <bam_dir> <plot_dir> <n_cores>

# Ben Ober-Reynolds

suppressPackageStartupMessages({
  library(magrittr)
  library(Matrix)
  library(ggrastr)
  library(reshape2)
})

## Read in command line arguments
args <- commandArgs(trailingOnly = FALSE)
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
args <- args[(grep("^--args", args, value = FALSE)[1] + 1): length(args)]

bamDir <- args[1]
plotDir <- args[2]
nCores <- as.numeric(args[3])


## Debugging
# setwd("/raid/USRdirs/ben/phase_sep/hexanediol_ATAC/20190514_mESC_gradient/")
# scriptPath <- "/home/ben/git_clones/bulkATAC/"
# 
# bamDir <- "output/sample/bams/deduped/"
# plotDir <- "output/sample/plots/qc/insert_size/"
# nCores <- 10


# Plotting configuration file
source(paste0(scriptPath, "/plotting_config.R"))

# Functions
source(paste0(scriptPath, "/bulkATAC_generics.R"))

#----------------------
# Plotting functions:
#----------------------


plotSingleFragmentLenDist <- function(df, windowSize = 600, color = "red"){
  # Plot a single fragment length distribution
  # First column of df should be fragment lengths, second should be fraction of insertions
  # Window size indicates max x value (i.e. plot from 0 to windowSize)
  p <- (
    ggplot(data = df, aes(x=df[,1], y=df[,2]))
    + geom_line(color=color, size=0.8)
    + scale_x_continuous(limits = c(0, windowSize), 
                         breaks = seq(0, windowSize, windowSize/4), 
                         expand = c(0,0))
    + scale_y_continuous(expand = c(0.0001,0.0001, 0.0005, 0.0005))
    + xlab("Insertion size")
    + ylab("Fraction inserts")
    + theme_bw()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), aspect.ratio = 0.621) 
    
  )
  p
}


plotMultipleFragmentLenDist <- function(df, windowSize = 600, id.vars = "width", cmap = NULL){
  melted <- melt(df, id.vars = id.vars)
  p <- (
    ggplot(data = melted, aes(x=melted[,1], y=melted[,3], col=melted[,2]))
    + geom_line(size=1)
    + scale_x_continuous(limits = c(0, windowSize), 
                         breaks = seq(0, windowSize, windowSize/4), 
                         expand = c(0,0))
    + scale_y_continuous(expand = c(0.0001,0.0001, 0.0005, 0.0005))
    + xlab("Insertion size")
    + ylab("Fraction inserts")
    + theme_bw()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 0.621) # 1 / golden ratio
    + guides(colour = guide_legend(title="", override.aes = list(size=2)))
    
  )
  if(!is.null(cmap)){
    nsamp <- length(unique(melted[,2]))
    cmap <- getColorMap(cmap, n = nsamp)
    p <- p + scale_color_manual(values = cmap)
  }
  p
}


#---------
# Script
#---------

# Set colorscheme for all plots
plotCmap <- cmaps_BOR$circus

## Start by plotting fragment length distribution per sample
# Identify processed bam files:
bamFiles <- list.files(path = bamDir, pattern = "\\.bam$", full.names = TRUE)

## Get width table per fragment file
message("Getting fragment length distribution per sample...")
dfList <- mclapply(seq_along(bamFiles), function(i){
  bf <- bamFiles[i]
  name <- strsplit(basename(bf), '\\.')[[1]][1]
  lenDF <- bamToFragmentRange(bf) %>% width %>% base::table(.) %>% as.data.frame
  colnames(lenDF) <- c("width", name)
  lenDF[,1] <- lenDF[,1] %>% as.character %>% as.numeric
  lenDF[,2] <- lenDF[,2] / sum(lenDF[,2])
  lenDF
}, mc.cores = nCores)

lenDF <- Reduce(function(x,y) merge(x,y, all=TRUE), dfList)
lenDF[is.na(lenDF)] <- 0.0

## Plot all single fragment length distributions:
if(!endsWith("/", x = plotDir)){
  plotDir <- paste0(plotDir, "/")
}

for(i in seq_along(bamFiles)){
  df <- lenDF[,c(1, i+1)]
  outFile <- paste0(strsplit(basename(bamFiles[i]), '\\.')[[1]][1], "_fragment_dist.pdf")
  pdf(paste0(plotDir, outFile))
  print(plotSingleFragmentLenDist(df)) # Print call is necessary here for some reason?
  dev.off()
}


## Plot all fragment length distributions
pdf(paste0(plotDir, "all_sample_fragment_dist.pdf"))
plotMultipleFragmentLenDist(lenDF, cmap = plotCmap)
dev.off()

