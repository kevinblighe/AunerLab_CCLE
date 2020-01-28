ccle_barplot <- function(
  ccledata,
  clinicaldata,
  keyword,
  gene,
  title,
  xlab,
  greyscale = TRUE,
  colour = NULL,
  titlesize = 24,
  axissize = 16) {

  # extract names of plasma cell myeloma / multiple myeloma lines
  lines <- sampleinfo[grep(keyword, clinicaldata$Hist_Subtype1),1]
  if (length(lines) < 2) {
    lines <- sampleinfo[grep(keyword, clinicaldata$Histology),1]

    if (length(lines) < 2) {
      stop('Error - too few cell-lines or nothing found')
    }
  }

  # filter expression data to only include multiple myeloma lines
  data <- ccledata[,which(colnames(ccledata) %in% lines)]

  # extract expression levels and convert to 'long' format
  ggdata <- reshape2::melt(data[which(rownames(data) %in% gene),])

  # change colnames and order by expression level
  colnames(ggdata) <- c('Line','Expression')
  ggdata <- ggdata[order(ggdata$Expression, decreasing = TRUE),]

  # tidy up cell line names
  #'$' matches end of field
  ggdata$Line <- gsub('_[A-Za-z0-9_]*$', '', ggdata$Line)
  ggdata$Line <- factor(ggdata$Line, levels = ggdata$Line)

  # Basic barplot
  p <- ggplot(data = ggdata, aes(x = Line, y = Expression, fill = Line)) +

    geom_bar(stat = 'identity') +

    # axis and main title(s) 
    xlab(xlab) +
    ylab(bquote(italic(.(gene))~RPKM~expression)) +
    ggtitle(title) +

    # set the size of the plotting window
    theme_bw(base_size=24) +

    # modify various aspects of the plot text and legend
    # NB - most of these, you will not have to edit
    theme(
      plot.title = element_text(angle=0, size = titlesize, face="bold", vjust = 1),
      axis.text.x = element_text(angle = 90, size = axissize, face = "bold", hjust = 1.0, vjust = 0.5),
      axis.text.y = element_text(angle = 0, size = axissize, face = "bold", vjust = 0.5),
      axis.title = element_text(size = 24, face = "bold"),

      # Legend
      legend.position = "none", # 'none' turns off the legend
      legend.background = element_rect(),
      legend.key = element_blank(),     #removes the border
      legend.key.size = unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text = element_text(size = 16),  #Text size
      title = element_text(size = 16)) +

    # change the size of the icons/symbols in the legend
    guides(colour = guide_legend(override.aes = list(size = 2.5)))

  # return the plot
  if (greyscale == TRUE) {
    p <- p + scale_fill_grey()
  }

  if (!is.null(colour)) {
    p <- p + scale_fill_manual(values = rep(colour, nrow(ggdata)))
  }

  return(p)
}
