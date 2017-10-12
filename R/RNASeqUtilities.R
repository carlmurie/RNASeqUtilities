#' User defined functions for use in RNASeq analysis
#'
#' Standardise commonly used commands
#'
#' Create commonly used functions to simplify RNASeq analysis markdown reports.
#' 
#' @docType package
#' @name RNASeqUtilities
#' @author Carl Murie
#' @import ggplot2
NULL


#' Produce a table of counts of all combinations of phenotype file.
#' 
#' Count number of instances for each combination of factor levels in relation to reference
#' factor. First column of return data.frame is total number of elements for each level 
#' of 'ref'. Used to calculate sample sizes for experimental designs.
#' 
#' @param dat data.frame of phenotype variables. Each row is one sample and each column
#'            is a factor of interest
#' @param ref character string or numeric column index referencing column for which 
#'            'factors' is tabulated on.
#' @param factors vector of column labels or column indices defining which factors are to
#'                be tabulated.
#'                
#' @return data.frame showing frequencies of 'factors' levels for 'ref'
#' @export
#' 
countSampleSizes <- function(dat, ref, factors) {
  
  ## check that 'dat' is non-null and a data.frame
  if(is.null(dat) | !is.data.frame(dat)) {
    stop("ERROR: 'dat' must be a data.frame")
  }
  
  if(is.null(ref) | is.null(factors)) {
    stop("ERROR: 'ref' or 'factors' is NULL")
  }
  
  ## check that if 'ref' is character it is also a column name in 'dat'
  if(is.character(ref) & !is.element(ref, colnames(dat))) {
    stop(paste("ERROR:", ref, "not a column name of 'dat'"))
  }
  
  ## check that if factors is character they are all column names in 'dat'
  if(is.character(factors) & sum(is.element(factors, colnames(dat))) != length(factors)) {
     stop("ERROR: There are elements of 'factors' that are not a column name of 'dat'")
  }
  
  ## check that if 'ref' is integer it is in range of 'dat' NOTE: I'm not checking for reals
  if(is.numeric(ref) & !is.element(ref, 1:ncol(dat))) {
    stop("ERROR: 'ref' is not a column index of 'dat'")
  }
  
  ## check that if factors are numeric they are in range of dat.
  if(is.numeric(factors) & sum(is.element(factors, 1:ncol(dat))) != length(factors)) {
    stop("ERROR: There are elements of 'factors' that are not a column indices of 'dat'")
  }
  
  ## check that 'ref' is not in 'factors'
  if(is.element(ref, factors)) {
     warning(paste(ref, "is an element of 'factors'. Is this what you want?"))
  }
  
  total <- table(dat[,ref])
  for(fac in factors) {
    total <- cbind(total, table(dat[,c(ref, fac)]))
  }
  return(data.frame(total))
} ## end countSampleSizes


#' Produce a set of descriptive graphs for DEG limma results (topTable)
#' 
#' Produce a histogram of p-values, volcano plot of logFC vs -log10 p-values and
#' a histogram of FDRs. Input the return value from topTable (limma) with column names
#' unchanged.
#' 
#' @param dat data.frame Return data.frame from topTable function in limma
#' @param title character string for title of graph.
#' @param fdrCut Significance cut-off value for FDR. Genes who pass cut-off will be 
#'                shown in red in volcano plot.
#'                
#' @return list with three ggplot2 elements - p-value histogram, volcano plot, and FDR
#'        histogram.        
#' @export        
#' 
createDEGGraphs <- function(dat, title=NULL, fdrCut=0.2) {
  
  topTableColNames <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  ## check that dat is non-null and a data.frame
  if(is.null(dat) | !is.data.frame(dat)) {
     stop("ERROR: dat must be a data.frame")
  }
  
  ## check that dat has the proper 6 columns
  if(ncol(dat) != 6 | sum(colnames(dat) !=topTableColNames) != 0) {
     stop("ERROR: dat is not constructed properly")
  }

  ## label significant genes based on FDRCut
  ## note that if all genes are significant then all genes are colored black not red.
  dat$FDR <- ifelse(dat$adj.P.Val > fdrCut, "NonSig", "Sig")   
  
  ## return object
  graphSet <- list()
  
  ## p-value histogram
  graphSet[[1]] <- ggplot(dat, aes_string(x="P.Value")) + geom_histogram(binwidth=.05) +
    xlab("p-value") + ggtitle(label=title) + theme_bw() + 
    theme(aspect.ratio=1, axis.text.x=element_text(size=10), 
          plot.title = element_text(size = 15)) + 
    scale_x_continuous(breaks=c(0, 0.5, 1))
  
  ## volcano plot
  graphSet[[2]] <- ggplot(dat, aes_string(x="logFC", y="-log10(P.Value)", color="FDR")) +
    geom_point(shape=1) + ggtitle("logFC vs -log10 p-value") + theme_bw() + 
    theme(axis.title.y=element_text(size=10), axis.text.x=element_text(size=6), 
          plot.title = element_text(size = 13)) + 
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="blue") + 
    geom_hline(yintercept=-log10(c(0.01)), linetype="dashed", color="blue") +
    theme(legend.position="none") +  scale_colour_manual(values=c("black", "red")) +
    theme(aspect.ratio=1)
  
  ## FDR histogram
  graphSet[[3]] <- ggplot(dat, aes_string(x="adj.P.Val")) + geom_histogram(binwidth=.05) +
    xlab("p-value") + ggtitle("FDR") + theme_bw() + 
    theme(aspect.ratio=1, plot.title = element_text(size = 15)) + 
    scale_x_continuous(breaks=c(0, 0.5, 1))
  
   return(graphSet)
} ## end createDEGGraphs



