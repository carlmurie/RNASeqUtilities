#' User defined functions for use in RNASeq analysis
#'
#' Standardise commonly used functionality in RNASeq analysis.
#'
#' Create commonly used functions to simplify RNASeq analysis markdown reports.
#'
#' @docType package
#' @name RNASeqUtilities
#' @author Carl Murie
#' @import ggplot2
#' @import limma
#' @import htmltools
#' @import DT
#' @import gridExtra
NULL



#'
#' Create a nicer kable table.
#'
#' Creates an html kable table with additional parameters, including a
#' multi-level column label if desired. This is designed to be used on a table
#' created by \link{countSampleSizes}. If 'levels' or 'labels' are NULL then
#' no header will be produced.
#'
#' @param tabley data.frame of sample sizes. Each column is a level and each row
#' is a level of the reference factor. The first column of the tabley will be total
#' and the remaining columns will be the levels.
#' @param labels vector of factor labels which will be the top level header.
#' @param levels vector of number of levels for each factor in 'labels'
#' @param title string for caption of table
#' @param align parameter for 'position' argument for kable_styling. default: left
#'
#' @return kable table
#'
#' @export
#'
headerKable <- function(tabley, labels=NULL, levels=NULL, title=NULL, align="left") {

  ## check that 'labels' and 'levels' are of same length
  if(length(labels) != length(levels)) {
    stop("parameters 'labels' and 'levels' must be of same length")
  }

  ## check that number of levels matchs dimenstion of 'tabley'
  if(!is.null(levels) && sum(levels) != ncol(tabley)-1) {
       stop("number of levels in 'levels' does not match number of levels in 'tabley'")
  }

  names(levels) <- labels
  tab <- kable(tabley, caption=title, align=rep("c", ncol(tabley))) %>%
         kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                       full_width=FALSE, position=align)

  ## add header if necessary
  if(!is.null(labels) && !is.null(levels)) {
     tab <- add_header_above(tab, c(" ", " ", levels))
  }

  tab
}


#'
#' Create markdown string showing deg results.
#'
#' Input values for parameters 'deg' and 'graphs' must be named deg and graphs.
#' TODO: explore htmltools:knitr_methods to output datatables instead.
#'
#' @param deg List of deg results from runDEG. Name of deg element will be used as title
#' @param graphs List of graphs from createDEGGraphs. Name of deg element will be used as title
#' @param doSlides Boolean: If TRUE then output markdown slides else produce standard markdown.
#' @param onlySigs Boolean: if TRUE output only genesets that pass the 'FDRCut' threshold
#'                 (output message if no significant gene sets) else output all genesets.
#'                 Default is FALSE.
#' @param cutOnPval Boolean: if TRUE use pvalCut for threshold else use FDRCut
#' @param FDRCut FDR threshold to use to determine significance. Default is 0.2.
#' @param pvalCut P value threshold to use to determine significance. Default is
#'                0.01
#' @param header Character string defining markdown header to ouput for each set
#'               of results. This is only used if doSlides=FALSE.
#' @return Vector of markdown strings
#'
#' @export
outputDEG <- function(deg, graphs, doSlides=TRUE, onlySigs=FALSE, cutOnPval=TRUE,
                      FDRCut=0.2, pvalCut=0.01, header="##") {

  ## check that 'deg' input is not null and a list
  if(is.null(deg) || !is.list(deg)) {
    stop("ERROR: 'deg' is either null or not a list.")
  }

  ## check that 'graphs' input is not null and a list
  if(is.null(graphs) || !is.list(graphs)) {
    stop("ERROR: 'graphs' is either null or not a list.")
  }

  ## check that doSlids is a boolean
  if(!is.logical(doSlides)) {
    stop("ERROR: 'doSlides' must be a boolean")
  }

  ## check that 'deg' and 'graphs' have the same length
  if(length(deg) != length(graphs)) {
     stop("ERROR: 'deg' and 'graphs' are not of same length")
  }

  for(i in 1:length(deg)) {

      graphy <- graphs[[i]]
      tab <- deg[[i]]

      ## reformat and rename deg results
      tab <- tab[,c("logFC", "P.Value", "adj.P.Val")]
      colnames(tab) <- c("logFC", "PValue", "FDR")
      namey <- names(deg)[i]

      ## output only significant results based on threshold
      if(onlySigs) {
          if(cutOnPval) {
            indo <- tab$PValue <= pvalCut
          } else {
             indo <- tab$FDR <= FDRCut
          }

          if(sum(indo) > 0) {
            tab <- signif(tab[indo,], 3)
          } else {
            tab <- -1
          }
      } ## end if onlySigs

      ## output slides
      if(doSlides) {

         ## output graphs
         knit_print.html(cat(paste("##", namey, " {.smaller}\n")))
         grid.arrange(graphy[[1]], graphy[[2]], graphy[[3]], nrow=1, ncol=3)
         knit_print.html(cat("  \n"))

         ## output deg table
         knit_print.html(cat(paste("##", namey, " {.smaller}\n")))

         if(is.null(dim(tab))) {
            knit_print.html(cat("No significant gene sets found  \n  \n"))
         }else {
            print(tagList(datatable(tab)))
         }
     } else { ## output normal markdown

        ## output graphs
        knit_print.html(cat(paste(header, namey, " \n  \n")))
        grid.arrange(graphy[[1]], graphy[[2]], graphy[[3]], nrow=1, ncol=3)
        knit_print.html(cat("  \n"))

        ## output deg table
        if(is.null(dim(tab))) {
          knit_print.html(cat("  \nNo significant gene sets found  \n  \n"))
        } else {
          print(tagList(datatable(tab)))
        }

     } ## end if doSlides

  } ## end for i

} ## end outputDEG


#'
#' Create markdown string showing gsea results.
#'
#' @param gsea List of gsea results from runGSEA. Name of gsea element will be used as title
#' @param doSlides Boolean: If TRUE then output markdown slides else produce standard markdown.
#' @param onlySigs Boolean: if TRUE output only genesets that pass the 'FDRCut' threshold
#'                 (output message if no significant gene sets) else output all genesets.
#'                 Default is FALSE.
#' @param FDRCut FDR threshold to use to determine significance. Default is 0.2.
#' @param header Character string defining markdown header to ouput for each set
#'               of results. This is only used if doSlides=FALSE.
#' @return Vector of markdown strings
#'
#' @export
outputGSEA <- function(gsea, doSlides=TRUE, onlySigs=TRUE, FDRCut=0.2,
                       header="##") {

  ## check that 'gsea' input is not null and a list
  if(is.null(gsea) || !is.list(gsea)) {
     stop("ERROR: 'gsea' is either null or not a list.")
  }

  ## check that 'doSlides' is a boolean
  if(!is.logical(doSlides)) {
     stop("ERROR: 'doSlides' must be a boolean")
  }

  ## check that 'onlySigs' is boolean
  if(!is.logical(onlySigs)) {
    stop("ERROR: 'onlySigs' parameter must be boolean")
  }

  ## check that 'FDRCut' is numeric between 0 and 1
  if(!is.numeric(FDRCut) || (FDRCut < 0 || FDRCut > 1)) {
    stop("ERROR: 'FDRCut' must be numeric and be between 0 and 1")
  }

  ## iterate across all results and output
  for(i in 1:length(gsea)) {

      tab <- gsea[[i]]
      namey <- names(gsea)[i]

      if(onlySigs) {
        indo <- tab$FDR <= FDRCut
        if(sum(indo) > 0) {
          tab <- cbind(tab[indo, c("NGenes", "Direction")],
                       signif(tab[indo, c("PValue","FDR")], 3))
        } else {
          tab <- -1
        }
      } ## end if onlySigs

      ## output slide markdown
      if(doSlides) {

           knit_print.html(cat(paste("##", namey, " {.smaller}\n")))

          if(is.null(dim(tab))) {
             knit_print.html(cat("No significant gene sets found  \n  \n"))
          }else {
            print(tagList(datatable(tab)))
          }
      } else { ## output normal markdown

          knit_print.html(cat(paste0("  \n", header, namey, "  \n  \n") ))
          if(is.null(dim(tab))) {
              knit_print.html(cat("  \nNo significant gene sets found  \n  \n"))
          } else {
             print(tagList(datatable(tab)))
          }
      } ## end of doSlides
  } ## end for i

} ## end outputGSEA



#'
#' Apply DEG analysis to voom normalised data for a specific model and coefficients.
#'
#' @param voomDat Eset of voom normalised data which must include the design. This is the
#'                output of voom function in limma.
#' @param coefs Vector of coefficients to test. They must map to a coefficient in the
#'              design matrix in 'voomDat'.
#' @param doRandomEffect Boolean indicating whether there is a random effect in the model. If
#'                       so then 'duplicateCorrelation' function is run based on 'blocker'.
#'                       The blocking factor must be identified by 'blocker' parameter.
#' @param blocker Vector showing factor of random effect. It is only used if doRandomEffect
#'                is TRUE.
#' @param labels Vector of character strings that labels each coefficient being tested. This
#'               allows for more human readable labels than the lmFit coefficients. Labels
#'               must map one to one to 'coefs'.
#' @param doWrite Boolean indicating whether to write the random correlation and fitted model
#'                to file defined in 'rdaPath'. Useful because computation time can be
#'                prohibitive when fitting model.
#' @param doRead Boolean indicating whether to read random corrleation and fitted model from
#'               file defined in 'rdaPath'. Data must have been previously saved with 'doWrite'
#'               parameter.
#' @param rdaPath Character string defining full path where rda object will be stored or read.
#'                Must be defined if 'doWrite' or 'doRead' are TRUE.
#'
#' @return List of DEG results. Each element in list corresponds to the results for the
#'         matching coefficients in 'coefs'.
#'
#' @export
runDEG <- function(voomDat, coefs, doRandomEffect=FALSE, blocker=NULL, labels=NULL,
                   doWrite=FALSE, doRead=FALSE, rdaPath=NULL) {

  ##check that all parameters are non-null
  if(is.null(voomDat) | is.null(coefs)) {
    stop("ERROR: passing null parameter to runDEG")
  }

  ## check that all elements of 'coefs' are in design matrix of 'voomDat'
  if(sum(is.element(coefs, colnames(voomDat$design))) != length(coefs)) {
    stop("ERROR: at least one element of 'coefs' is not found 'voomDat' design")
  }

  ## check that blocker is defined if 'doRandomEffect' is TRUE
  if(doRandomEffect && is.null(blocker)) {
     stop("ERROR: 'blocker' must be defined if 'doRandomEffect' is TRUE.")
  }

  ## check that blocker is same dimension as voomDat
  if(doRandomEffect && (length(blocker) != ncol(voomDat))) {
     stop("ERROR: 'blocker is not of same length as 'voomDat' samples")
  }

  ## check that labels is same length as coefs.
  if(!is.null(labels) && length(labels) != length(coefs)) {
    stop("ERROR: 'labels' is not of same length as 'coefs'")
  }

  if((doWrite | doRead) & is.null(rdaPath)) {
     stop("ERROR: 'rdaPath' must be defined if 'doWrite' or 'doRead' are TRUE.")
  }

  ## set list labels to either labels or coefs
  labbs <- labels
  if(is.null(labels)) {
      labbs <- coefs
  }

  ## read in previously saved data
  if(doRead) {
      if(!file.exists(rdaPath)) {
         stop("ERROR: 'rdaPath' is not found.")
      }
      load(rdaPath)

  } else {

      ## calculate random effects correlation if necessary
      ranCor <- NULL
      if(doRandomEffect) {
         ranCor <- duplicateCorrelation(voomDat, block=blocker)$consensus.correlation
      }

      ## fit coefficients
      fit1 <- lmFit(voomDat, block=blocker, correlation=ranCor)
      fit2 <- eBayes(fit1, trend=FALSE)

      ## write data to file for future use
      if(doWrite) {
         save(ranCor, fit2, file=rdaPath)
      }

  } ## end if doRead

  allTopTable <- list()
  labInd <- 0
  for(c in coefs) {
       labInd <- labInd + 1
       allTopTable[[labbs[labInd]]] <- topTable(fit2, number=Inf, coef=c, sort.by="P")
  } ##end for c

  return(allTopTable)

} ## end runDEG



#'
#' Apply GSEA to a voom normalized data set for a specific model and coefficients.
#'
#' @param voomDat Eset of voom normalised data which must include the design. This is the
#'                output of voom function in limma.
#' @param geneIndices List of gene indices corresponding to gene sets. It is created by
#'                     ids2indices function in limma.
#' @param coefs Vector of coefficients to test. They must map to a coefficient in the
#'              design matrix in 'voomDat'.
#' @param labels Vector of character strings that labels each coefficient being tested. This
#'               allows for more human readable labels than the lmFit coefficients. Labels
#'               must map one to one to 'coefs'.
#' @param onlySigs Boolean: if TRUE return only genesets that pass the 'FDRCut' threshold
#'                 (return -1 if no significant gene sets) else return all genesets. Default
#'                is FALSE.
#' @param FDRCut FDR threshold to use to determine significance. Default is 0.2.
#'
#' @return List of camera results. Each element is the result for a single coefficient from
#'         'coefs'
#'
#' @export
runGSEA <- function(voomDat, geneIndices, coefs, labels=NULL, onlySigs=FALSE, FDRCut=0.2) {

  ##check that all parameters are non-null
  if(is.null(voomDat) | is.null(geneIndices) | is.null(coefs)) {
     stop("ERROR: passing null parameter to runGSEACoef")
  }

  ## check that all elements of 'coefs' are in design matrix of 'voomDat'
  if(sum(is.element(coefs, colnames(voomDat$design))) != length(coefs)) {
    stop("ERROR: at least one element of 'coefs' is not found 'voomDat' design")
  }

  ## check that labels is same length as coefs.
  if(!is.null(labels) && length(labels) != length(coefs)) {
     stop("ERROR: 'labels' is not of same length as 'coefs'")
  }

  ## check that 'onlySigs' is boolean
  if(!is.logical(onlySigs)) {
    stop("ERROR: 'onlySigs' parameter must be boolean")
  }

  ## check that 'FDRCut' is numeric between 0 and 1
  if(!is.numeric(FDRCut) || (FDRCut < 0 || FDRCut > 1)) {
     stop("ERROR: 'FDRCut' must be numeric and be between 0 and 1")
  }

  ## set list labels to either labels or coefs
  labbs <- labels
  if(is.null(labels)) {
       labbs <- coefs
  }

  ## apply camera to each coefficient in coefs
  allOut <- list()
  labIndex <- 0
  for(coef in coefs) {
    labIndex <- labIndex+1
    res <-  camera(voomDat, index=geneIndices, contrast=coef, sort=TRUE)

    ## return only significant genes, if none return -1
    if(onlySigs) {
      indo <- res$FDR <= FDRCut
      if(sum(indo) > 0) {
        res <- cbind(res[indo, c("NGenes", "Direction")],
                                signif(res[indo, c("PValue","FDR")], 3))
      } else {
        res <- -1
      }
    } ## end if onlySigs

    allOut[[labbs[labIndex]]] <- res
  }## end for coef

  return(allOut)

} ## end runGSEA



#'
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


#'
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
  graphSet[[3]] <- ggplot(dat, aes_string(x="P.Value", y="adj.P.Val")) + geom_step() +
    xlab("FDR") + ylab("P value") +  ggtitle("FDR v P value") + theme_bw() +
    theme(aspect.ratio=1, plot.title = element_text(size = 15)) + ylim(0,1) + xlim(0,1)

   return(graphSet)
} ## end createDEGGraphs



