library(ggplot2)

test_that("createDEGGraphs handles various input", {
  
  ## test basic functionality
  set.seed(1)
  dat <- data.frame(matrix(runif(600, 0, 1), ncol=6)) 
  colnames(dat) <-  c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  expect_is(createDEGGraphs(dat), "list")
  
  ## test for bad input
  colnames(dat)[1] <- "fred"
  expect_error(createDEGGraphs(dat))
  
  ## test for all significant FDR results (this incorrectly results in all black points)
  colnames(dat)[1] <- "logFC"
  dat$adj.P.Val <- rep(.05, 100)
  expect_is(createDEGGraphs(dat), "list")
  
})