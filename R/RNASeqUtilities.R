#` User defined functions for use in RNASeq analysis

#` Standardise commonly used commands

#` Create commonly used functions to simplify RNASeq analysis markdown reports.

#` Produce a set of descriptive graphs for DEG limma results (topTable)

#` Produce a histogram of p-values, volcano plot of logFC vs -log10 p-values and
#' a histogram of FDRs. Input the return value from topTable (limma) with column names
#' unchanged.
#' 
#' @param dat data.frame Return data.frame from topTable function in limma
#' @param fdrCut Significance cut-off value for FDR. Genes who pass cut-off will be 
#'                shown in red in volcano plot.
#'                
#' @return list with three ggplot2 elements - p-value histogram, volcano plot, and FDR
#'        histogram.                
createDEGGraphs <- function(dat, fdrCut=0.2) {
  
  topTableColNames <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  ## check that dat is non-null and a data.frame
  if(is.null(dat) | !is.data.frame(dat)) {
     stop("ERROR: dat must be a data.frame")
  }
  
  if(sum(colnames(dat) !=topTableColNames) != 0) {
     stop("ERROR: dat has invalid column names. They must be the same as topTable output")
  }
  
  ## label significant genes based on FDRCut
  dat$FDR <- ifelse(dat$adj.P.Val > fdrCut, "NonSig", "Sig")   
  
  ## return object
  graphSet <- list()
  
  ## p-value histogram
  graphSet[[1]] <- ggplot(dat, aes(x=P.Value)) + geom_histogram(binwidth=.05) +
    xlab("p-value") + ggtitle(label=labb) + theme_bw() + 
    theme(aspect.ratio=1, axis.text.x=element_text(size=10), 
          plot.title = element_text(size = 15)) + 
    scale_x_continuous(breaks=c(0, 0.5, 1))
  
  ## volcano plot
  graphSet[[2]] <- ggplot(dat, aes(x=logFC, y=-log10(P.Value), color=FDR)) +
    geom_point(shape=1) + ggtitle("logFC vs -log10 p-value") + theme_bw() + 
    theme(axis.title.y=element_text(size=10), axis.text.x=element_text(size=6), 
          plot.title = element_text(size = 13)) + 
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="blue") + 
    geom_hline(yintercept=-log10(c(0.01)), linetype="dashed", color="blue") +
    theme(legend.position="none") +  scale_colour_manual(values=c("black", "red")) +
    theme(aspect.ratio=1)
  
  ## FDR histogram
  graphSet[[3]] <- ggplot(ress1, aes(x=adj.P.Val)) + geom_histogram(binwidth=.05) +
    xlab("p-value") + ggtitle("FDR") + theme_bw() + 
    theme(aspect.ratio=1, plot.title = element_text(size = 15)) + 
    scale_x_continuous(breaks=c(0, 0.5, 1))
  
   return(graphSet)
} ## end createDEGGraphs


