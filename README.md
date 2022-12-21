
# shRNAtools

<!-- badges: start -->
<!-- badges: end -->


The goal of shRNAtools is to create the oligos to knockdown specific circRNAs  

Contents
========

 * [Why?](#why)
 * [Installation](#installation)
 * [Usage](#usage)
 
### Why?

I wanted a tool that allows you to:

+ Get the specific _circRNA junction_.
+ Design shRNAs against circRNAs.
+ Design shRNAs shifts.

+ Get the specific _mRNA junction_.
+ Design shRNAs against linear RNAs.
+ Design shRNAs shifts.

### Installation
---

You can install the development version of shRNAtools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ipatop/shRNAtools")
```

### Usage

This is a basic example which shows you how to use it

# Run circRNA example

Input should look like this
```{r}
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "circs_totest.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/circs_totest.txt",header = T)
```

Run to create an output table
```{r}
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "circs_totest.txt", package = "shRNAtools")

circRNASpliceOligoDesigner(input_coordinates  = "../test/circs_totest.txt",output = "../test/New_out.tsv")
```

Run to create a DataFrame, if writetab = F
```{r}
circOligos<-circRNASpliceOligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)
```

Output
```{r}
circOligos
```

# Run linear RNA example

Input should look like this
```{r}
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "lin_totest.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/lin_totest.txt",header = T)
```

Run to create an output table
```{r}
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "lin_totest.txt", package = "shRNAtools")

linRNASpliceOligoDesigner(input_coordinates  = "../test/circs_totest.txt",output = "../test/New_out_lin.tsv")
```

Run to create a DataFrame, if writetab = F
```{r}
linOligos<-linRNASpliceOligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)
```

Output
```{r}
linOligos
```
