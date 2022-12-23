## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#install library
#devtools::install_github("ipatop/shRNAtools")
library(shRNAtools)

## -----------------------------------------------------------------------------
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "circs_totest.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/circs_totest.txt",header = T)

## -----------------------------------------------------------------------------
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "circs_totest.txt", package = "shRNAtools")

circRNASpliceOligoDesigner(input_coordinates  = "../test/circs_totest.txt",output = "../test/New_out.tsv")

## -----------------------------------------------------------------------------
circOligos<-circRNASpliceOligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)

## -----------------------------------------------------------------------------
circOligos

## -----------------------------------------------------------------------------
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "lin_totest.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/lin_totest.txt",header = T)

## -----------------------------------------------------------------------------
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "lin_totest.txt", package = "shRNAtools")

linRNASpliceOligoDesigner(input_coordinates  = "../test/circs_totest.txt",output = "../test/New_out_lin.tsv")

## -----------------------------------------------------------------------------
linOligos<-linRNASpliceOligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)

## -----------------------------------------------------------------------------
linOligos

