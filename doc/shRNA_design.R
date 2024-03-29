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
as.data.frame(circOligos)

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
as.data.frame(linOligos)

## -----------------------------------------------------------------------------
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "introns.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/introns.txt",header = T)

## -----------------------------------------------------------------------------
RCMs <- find_rcm(sp = "mm39",intron_coordinates  = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/introns.txt",outBlast = "outBlast.tsv",outBlast.text = "outBlast.txt",blastn = "/Users/inespatop/Documents/Scripts/blast-2.14.0-h23b05c9_0/bin/blastn",ret = T)

## -----------------------------------------------------------------------------
as.data.frame(RCMs)

