
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

+ Find _reverse complementary matches_ (RCM).
+ Output a comprehensive table with the results.
+ Generate grafical outputs with the RCM matches.

### Installation
---

You can install the development version of shRNAtools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ipatop/shRNAtools")
```

### Usage

This package contains two functions that will generate a Data Frame or tab separated file with the original table plus the shRNA Oligo design appended. 

The input should be a tab separated table with circRNA or linear RNA splicing junction coordinates and gene names in the format: 

<table>

Name    | Chr   | Start    | End      | Strand
------  | ----- | -------- | -------- | -------
circMbl | chr2R	| 17275410 | 17276063 |    +

</table>

This table can have an optional column with the Strand

So far the following species are available: 

+ Fly: 
  + dm3
  + dm6 
+ Human: 
  + hs19
  + hg38
+ Mice: 
  + mm10 
  + mm39
+ Rat: 
  + rn4
  
# Run circRNA example

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
# Find RCMs (reverse complementary repeats) in introns flanking circRNAs

Input should look like this
```{r}
#Copy and paste the full path of the test data into the readtable argument
system.file("extdata", "introns.txt", package = "shRNAtools")

read.table(file = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/introns.txt",header = T)
```

Run to create an output table
```{r}
RCMs <- find_rcm(sp = "mm39",intron_coordinates  = "/Users/inespatop/Documents/shRNA_design_circRNA/shRNAtools/inst/extdata/introns.txt",outBlast = "outBlast.tsv",outBlast.text = "outBlast.txt",blastn = "/Users/inespatop/Documents/Scripts/blast-2.14.0-h23b05c9_0/bin/blastn",ret = T)
```

**Output**
Table
```{r}
as.data.frame(RCMs)
```
