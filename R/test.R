setwd("~/Documents/shRNA_design_circRNA/shRNAtools/")

require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
tx<-TxDb.Dmelanogaster.UCSC.dm6.ensGene
require(BSgenome.Dmelanogaster.UCSC.dm6)
gen<-BSgenome.Dmelanogaster.UCSC.dm6



#read the input dat
coordinates <- read.delim("./test/lin_totest.tsv")

coordinates$junction <- paste0(as.data.frame(getSeq(gen, coordinates$Chr,coordinates$Start-20,coordinates$Start))[,1],as.data.frame(getSeq(gen, coordinates$Chr,coordinates$End,coordinates$End+20))[,1])

coordinates$junction


