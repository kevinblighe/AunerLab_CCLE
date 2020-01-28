ccle_correlation <- function(
  ccledata,
  clinicaldata,
  keyword,
  gene1,
  gene2,
  title,
  titlesize = 24,
  axissize = 16,
  corDigits = 3,
  pvalDigits = 4) {

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
  ggdata <- data.frame(
    t(data[which(rownames(data) %in% gene1),]),
    t(data[which(rownames(data) %in% gene2),]))
  colnames(ggdata) <- c('gene1', 'gene2')

  # calculate correlation coefficient and p-value
  r <- round(cor(ggdata[,1], ggdata[,2], method = 'spearman'), digits = corDigits)
  pval <- cor.test(ggdata[,1], ggdata[,2], method = 'spearman')$p.value
  pval <- ifelse(pval<0.0001, "p < 0.0001", paste0("p, ", round(pval, digits = pvalDigits)))

  # remove the dollar sign from the input variables (for visualisation)
  xlab <- gsub('\\$', '', gene1)
  ylab <- gsub('\\$', '', gene2)

  p <- ggplot(ggdata, aes(x = gene1, y = gene2)) +

    geom_point() +

    geom_smooth(method="lm", formula=y~x, level=0.95) +
    stat_smooth(method="lm", fullrange=TRUE, level=0.95, colour="red2") +

    # axis and main title(s) 
    xlab(bquote(italic(.(xlab))~expression)) +
    ylab(bquote(italic(.(ylab))~expression)) +
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
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +

    # add lines through the origin (0,0)
    #geom_hline(yintercept=0, linetype="dashed", color="black") +
    #geom_vline(xintercept=0, linetype="dashed", color="black") +

    # add the correlation coefficient and p-value as captions
    labs(caption = paste('Spearman rho, ', r, ' | ', pval, sep=""))

  return(p)
}
