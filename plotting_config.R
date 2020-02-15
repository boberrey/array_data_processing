# Plotting settings

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(magrittr)
  library(ggplot2)
  library(circlize)
})


theme_BOR <- function(base_size=14, base_family="Helvetica", border = TRUE) {
  library(grid)
  library(ggthemes)
  # Should plots have a bounding border?
  if(border){
    panel.border <- element_rect(fill = NA, color = "black", size = 0.7)
    axis.line <- element_blank()
  }else{
    panel.border <- element_blank()
    axis.line <- element_line(color = "black", size = 0.5)
  }
  
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = panel.border,
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = axis.line,
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))
  
}

scale_fill_BOR <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_BOR <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}




#---------------------------
# Colormaps
#---------------------------

cmaps_BOR <- list(
  ## Sequential colormaps:
  solarExtra = c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', "7"='#FDB31A',
                  "8"= '#E42A2A', "9"='#A31D1D'),
  
  # Jeff called this 'viridis', although its definitly different than the built in viridis
  sunrise = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  
  # Blue - Green - Gold/orange
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'),
  
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),
  
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  
  FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"), #wesanderson
  
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"), #wesanderson
  
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"), #wesanderson
  
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  
  ## Diverging colormaps:
  
  coolWarm = c("#00004c", "#0000ff", "#ffffff", "#ff0000", "#7f0000"),
  
  
  ## Qualitative colormaps:
  
  # see: https://carto.com/carto-colors/
  # I've removed the grey color from each palette
  cartoPrism = c("1"='#7F3C8D', "2"='#11A579', "3"='#3969AC', "4"='#F2B701', "5"='#E73F74', "6"='#80BA5A', "7"='#E68310', 
                  "8"='#008695', "9"='#CF1C90', "10"='#F97B72', "11"='#4B4B8F'),
  
  cartoSafe = c("1"='#88CCEE', "2"='#CC6677', "3"='#DDCC77', "4"='#117733', "5"='#332288', "6"='#AA4499', "7"='#44AA99',
                 "8"='#999933', "9"='#882255', "10"='#661100', "11"='#6699CC'),
  
  cartoBold = c("1"='#7F3C8D' ,"2"='#11A579', "3"='#3969AC', "4"='#F2B701', "5"='#E73F74', "6"='#80BA5A', "7"='#E68310',
                 "8"='#008695', "9"='#CF1C90', "10"='#f97b72', "11"='#4b4b8f'),
  
  cartoAntique = c("1"='#855C75', "2"='#D9AF6B', "3"='#AF6458', "4"='#736F4C', "5"='#526A83', "6"='#625377', "7"='#68855C',
                    "8"='#9C9C5E', "9"='#A06177', "10"='#8C785D', "11"='#467378'),
  
  cartoPastel = c("1"='#66C5CC', "2"='#F6CF71', "3"='#F89C74', "4"='#DCB0F2', "5"='#87C55F', "6"='#9EB9F3', "7"='#FE88B1',
                   "8"='#C9DB74', "9"='#8BE0A4', "10"='#B497E7', "11"='#D3B484'),
  
  cartoVivid = c("1"='#E58606', "2"='#5D69B1', "3"='#52BCA3', "4"='#99C945', "5"='#CC61B0', "6"='#24796C', "7"='#DAA51B',
                  "8"='#2F8AC4', "9"='#764E9F', "10"='#ED645A', "11"='#CC3A8E'),
  
  
  # 15 color
  circus = c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", 
              "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),
  
  iron_man = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
               "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),
  
  # The following 3 were designed by Ryan Corces. Jeff often uses this first one, 'stallion'
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),
  
  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
           "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),
  
  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
            "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13"),
  
  
  # Massive colorset for overly large things (30):
  tooManyClusters = c("#d52126", "#88ccee", "#fee52c", "#117733", "#cc61b0", "#99c945", "#2f8ac4", "#332288", "#e68316", "#661101", 
                      "#f97b72", "#ddcc77", "#11a579", "#6a00d8", "#635a63", "#d89896", "#004b00", "#b49ace", "#dcac6c", "#5e1645", 
                      "#544c19", "#3c807b", "#975064", "#9b1914", "#f3b420", "#8fb978", "#df3689", "#406582", "#8e8127", "#89288f")
)


#--------------------------
# Colormap helper functions
#--------------------------

mostDifferentColors <- function(cols, n = 20, firstpick = NULL, colorspace = "Lab", startingCols = NULL){
  stopifnot(length(cols) > n)
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  
  # First, convert sRGB to another colorspace (more 'perceptually uniform' colorspace)
  rgbCols <- t(col2rgb(cols))
  conv <- grDevices::convertColor(rgbCols, from = "sRGB", to = colorspace, scale.in = 255)
  
  # Now select n 'furthest neighbors' colors
  # (Note: I'm pretty sure this is an unsolved problem right now. 
  # This performs an iterative procedure for picking colors that maximize
  # 'distance' to already selected colors. The first color is picked randomly.)
  
  
  # If starting cols provided, add these to the list of picked cols (first pick and starting cols cannot be used together)
  if(!is.null(startingCols)){
    stConv <- grDevices::convertColor(t(col2rgb(startingCols)), from = "sRGB", to = colorspace, scale.in = 255)
    pickedColors <- list()
    for(i in seq_len(nrow(stConv))){
      pickedColors[[i]] <- stConv[i,]
    }
    remainingColors <- conv
  }else if(!is.null(firstpick)){
    message("using first pick...")
    stopifnot(firstpick %in% cols)
    idx <- which(cols == firstpick)
    pickedColors <- list(conv[idx,])
    remainingColors <- conv[-idx,]
  }else{
    message("using random start...")
    idx <- sample(1:nrow(conv), 1)
    pickedColors <- list(conv[idx,])
    remainingColors <- conv[-idx,]
  }
  pickedLen <- length(pickedColors)
  
  
  # Now, iteratively add the furthest color from the selected colors
  for(i in seq(pickedLen, n - 1)){
    distList <- list()
    for(j in seq_along(pickedColors)){
      colJ <- pickedColors[[j]]
      distMat <- dist(rbind(colJ, remainingColors), method = "euclidean") %>% as.matrix
      distList[[j]] <- distMat[2:nrow(distMat),1]
    }
    # What we actually want to maximize is the minimum distance between each color
    distMat <- do.call(cbind, distList)
    distMins <- apply(distMat, 1, FUN = min)
    idx <- which(max(distMins) == distMins)
    pickedColors[[i + 1]] <- remainingColors[idx,]
    remainingColors <- remainingColors[-idx,]
  }
  pickedLab <- do.call(rbind, pickedColors)
  pickedRgb <- round(grDevices::convertColor(pickedLab, from = colorspace, to = "sRGB", scale.out = 255),0)
  hex <- apply(pickedRgb, 1, rgb2hex)
  hex
}


pairwiseColorInterpolations <- function(cols, colorspace = "Lab"){
  # Get all pairwise interpolations between a vector of colors
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  interpolate <- function(c1, c2, colorspace){
    rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
  }
  paired <- sapply(cols, function(x) sapply(cols, function(y) interpolate(x, y, colorspace)))
  unique(as.vector(paired))
}


getColorMap <- function(cmap, n){
  stopifnot(n >= 1)
  # Return a character vector of n colors based on
  # the provided colormap. If n > length(cmap), do
  # some smart interpolation to get enough colors
  names(cmap) <- NULL # Having names on colors causes problems for some plotting routines
  if(length(cmap) < n){
    cmap <- mostDifferentColors(
      pairwiseColorInterpolations(cmap), 
      colorspace = "Apple RGB", n = n, startingCols = cmap
    )
  }
  cmap[1:n]
}

plotColorMap <- function(cols){
  # Plot each of the colors in a colormap
  cols <- base::unname(cols)
  n <- length(cols)
  df <- data.frame(
    x = seq_len(n),
    y = rep(1, n),
    z = factor(seq_len(n))
  )
  p <- (
    ggplot(df, aes(x=x,y=y,color=z)) 
    + geom_tile(aes(fill=z))
    + theme_BOR()
    + scale_color_manual(values = cols)
    + scale_fill_manual(values = cols)
  )
  p
}


# This is used primarily for making colormaps for ComplexHeatmap
makeColFun <- function(start, end, cmap, midpoint = NULL){
  # Make a color ramp function from provided start and end breaks,
  # and optionally a midpoint
  cmapLen <- length(cmap)
  if(!is.null(midpoint)){
    interpolate <- function(c1, c2, colorspace = "Lab"){
      rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
    }
    if(length(cmap) %% 2 == 0){
      # Interpolate middle colors if necessary to get midpoint
      preMidIdx <- floor(cmapLen / 2)
      midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
      cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
      cmapLen <- length(cmap)
    }
    midIdx <- ceiling(cmapLen / 2)
    breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
  } else {
    breaks <- seq(start, end, length.out = cmapLen)
  }
  colorRamp2(breaks, cmap)
}



#-------------------
# Plotting functions
#-------------------

plotUMAP <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = ""){
  # Given a df containing the UMAP x and y coords and a third column, 
  # plot the UMAP
  p <- (
    ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
    + geom_point_rast(size = 0.5)
    + theme_BOR()
    + theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1
    )
    + xlab("UMAP1")
    + ylab("UMAP2")
    + guides(fill = guide_legend(title=covarLabel), 
             colour = guide_legend(override.aes = list(size=5)))
  )
  # If colormap provided, update colors
  if(!is.null(cmap)){
    # Remove names (causes problems?)
    names(cmap) <- NULL
    if(dataType == "qualitative"){
      # check to make sure you have enough colors for qualitative mapping
      nvals <- length(unique(df[,3]))
      cmap <- getColorMap(cmap, n = nvals)
      p <- p + scale_color_manual(values = cmap,
                                  name = covarLabel)
    }else{
      p <- p + scale_color_gradientn(colors = cmap,
                                     name = covarLabel)
    }
  }
  p
}


qcBarPlot <- function(df, cmap = NULL){
  # Plot a bar plot
  nsamp <- nrow(df)
  # Fix colormap if provided
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = nsamp)
  }else{
    cmap <- "blue"
  }
  
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", color = cmap, fill = cmap)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_bw()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
}


qcViolinPlot <- function(df, cmap = NULL, makeLog = FALSE){
  # Plot a violin plot
  nsamp <- length(unique(df[,1]))
  aspectRatio <- 6/nsamp
  # Assume that the first column is the sample and the second column is the variable of interest
  if(makeLog){
    df[,2] <- log10(df[,2])
    colnames(df)[2] <- paste0("log10 ", colnames(df)[2]) 
  }
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], color = df[,1]))
    + geom_violin(aes(fill = df[,1]))
    + geom_boxplot(width = 0.8, alpha = 0)
    + scale_color_manual(values = cmap)
    + scale_fill_manual(values = alpha(cmap, 0.2))
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_bw()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
  
  # Adjust colors if necessary:
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = nsamp)
  }else{
    cmap <- rep("blue", times = nsamp)
  }
  p <- suppressMessages(p + scale_color_manual(values = cmap))
  p <- suppressMessages(p + scale_fill_manual(values = alpha(cmap, 0.3)))
  p
}


# Time to write a nice heatmap wrapper
BORHeatmap <- function(
  mat, # Data to plot (matrix or dataframe)
  limits = NULL, # Enforced limits for colormap (2 dimensional array)
  clusterCols = TRUE, # Should columns be clustered
  clusterRows = TRUE, # Should rows be clustered
  labelCols = FALSE, # Should columns be labeled
  labelRows = FALSE, # Should rows be labeled
  dataColors = NULL, # Colormap for plotting data
  dataColorMidPoint = NULL, # The data value to be the middle of the color map
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.15,
  useRaster = TRUE, # Should heatmap be rasterized
  rasterDevice = "CairoPNG", # For some reason 'png' doesn't work on cluster lustre
  rasterQuality = 5, # Raster quality. Higher is {better?}
  fontSize = 6, # Font size for labels
  showColDendrogram = FALSE, # Should the column dendrogram be shown
  showRowDendrogram = FALSE, # Should the row dendrogram be shown
  borderColor = NA, # Color for lines between cells
  mapname = " " # 'Name' to give heatmap
){
  
  #Packages
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  
  # Make sure mat is actually a matrix
  if(!is.matrix(mat)){
    message("'mat' needs to be a matrix. Converting...")
    mat <- as.matrix(mat)
  }
  
  # Prepare color function
  if(!is.null(limits)){
    ll <- limits[1]
    ul <- limits[2]
  }else{
    ll <- min(mat)
    ul <- max(mat)
  }
  # If no colormap provided, use solarExtra
  if(is.null(dataColors)){
    dataColors <- c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', 
                    "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', 
                    "7"='#FDB31A', "8"= '#E42A2A', "9"='#A31D1D')
  }
  dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)
  
  message("Preparing Heatmap...")
  hm <- Heatmap(
    # Main components:
    matrix = mat,
    name = mapname,
    col = dataColFun,
    
    # Legend options:
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = unit(5, "cm")
    ),
    rect_gp = gpar(col = borderColor), 
    
    # Column options:
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSize), 
    
    # Row options:
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    row_names_gp = gpar(fontsize = fontSize), 
    
    # Raster info:
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality
  )
  
  # Add row labels if provided:
  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    hm <- hm + rowAnnotation(
      link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
      width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
    )
  }
  
  return(hm)
}






