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

test_that("countSampleSizes handles various input", {
  
  dat <- data.frame(treat=rep(c("A", "B"), each=20), type=rep(c("C", "D"), 20),
                    site=rep(c("W", "X", "Y", "Z"), 10),
                    gender=rep(c("M", "F", "F", "M"), 10)) 
  
   truth <- data.frame(rbind(c(20,10,10,5,5,5,5,10,10), c(20,10,10,5,5,5,5,10,10)))
   dimnames(truth) <- list(c("A", "B"), c("total", "C","D", "W", "X", "Y", "Z", "F", "M"))
  
   ## test basic functionality
   expect_equal(countSampleSizes(dat, ref="treat", factors=c("type", "site", "gender")), truth)
   expect_equal(countSampleSizes(dat, ref=1, factors=2:4), truth)
   
   ## test bad inputs
   expect_error(countSampleSizes(dat, ref="fred", factors=c("type", "site", "gender")), 
                "ERROR: fred not a column name of 'dat'")
   
   expect_error(countSampleSizes(dat, ref="treat", factors=c("fred", "site", "gender")), 
                "ERROR: There are elements of 'factors' that are not a column name of 'dat'")
   expect_error(countSampleSizes(dat, ref=10, factors=2:4), 
                "ERROR: 'ref' is not a column index of 'dat'")
   expect_error(countSampleSizes(dat, ref=1, factors=2:10), 
             "ERROR: There are elements of 'factors' that are not a column indices of 'dat'")
   expect_warning(countSampleSizes(dat, ref="treat", factors=c("treat","type","site","gender")), 
                "treat is an element of 'factors'. Is this what you want?")
                
})
   
   
   
   