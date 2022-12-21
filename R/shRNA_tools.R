#Ines Patop, 2022
#The things you can use for design of shRNAs

#function to get "not in" from %in% in dplyr
'%!in%' <- function(x,y)!('%in%'(x,y))

#' Generate shRNAs against circRNA back-splice junctions
#'
#'This script contains ONE funciton that will read a table with circRNA/exon-exon splicing junction coordinates and genenames in the format: Name Chr Start End.
#'Optional:a column named Strand
#'So far the following species are available: fly: dm3, dm6 and human: hs19, mice: mm10 and rat: rn4
#'
#'
#' @param sp character with Name of the genome to use, available are: dm6, dm3, hg19, hg38, mm10, rn4. Default is dm6
#' @param input_coordinates tab separated document with the coordinate of the circRNAs to generate shRNAs in TSV format tab seppatated. It requires a column with Chromosome, Start and End.
#' @param writetab a boolean, if writetab=T, then the function will not write any table, instead it will return the table to R
#' @param output if writetab = T, it will create the output in this path, default is "shRNAs.tsv"
#'
#' @param shift boolean, if shift=T, return also 5´and 3´shift oligos sequences
#'
#' @return if writetab=F, It generates a DataFrame with the original table with the shRNA deisng appended
#'
#'
#' @export
#'
#' @examples
#'
#' #circRNAoligoDesigner(input_coordinates = "../test/circs_totest.txt", writetab=F)
#'
#'
circRNASpliceOligoDesigner<-function(sp="dm6",input_coordinates="circ_totest.txt", writetab=T,output="shRNAs.tsv",shift=T){

  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(input_coordinates))
  )


  if (sp =="dm6") {
    require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm6.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm6)
    gen<-BSgenome.Dmelanogaster.UCSC.dm6

  }  else if (sp =="dm3") {
    require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm3.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm3)
    gen<-BSgenome.Dmelanogaster.UCSC.dm3

  } else if (sp =="hg19") {
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg19.knownGene
    require(BSgenome.Hsapiens.UCSC.hg19)
    gen<-BSgenome.Hsapiens.UCSC.hg19

  } else if (sp =="hg38") {
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg38.knownGene
    require(BSgenome.Hsapiens.UCSC.hg38)
    gen<-BSgenome.Hsapiens.UCSC.hg38

  } else if (sp =="mm10") {
    require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    tx<-TxDb.Mmusculus.UCSC.mm10.knownGene
    require(BSgenome.Mmusculus.UCSC.mm10)
    gen<-BSgenome.Mmusculus.UCSC.mm10

  }  else if (sp =="rn4") {
    require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    tx<-TxDb.Rnorvegicus.UCSC.rn4.ensGene
    require(BSgenome.Rnorvegicus.UCSC.rn4)
    gen<-BSgenome.Rnorvegicus.UCSC.rn4

  }

  #read the input data
  coordinates<-read.delim(input_coordinates)


  #create primers

  coordinates$junction_long <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End-20,coordinates$End))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start,coordinates$Start+20))[,1])

  if(any(coordinates$Strand=="-")){
    coordinates$junction_long[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction_long[coordinates$Strand=="-"])))
  }


  #create junction
  coordinates$junction <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End-9,coordinates$End))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start,coordinates$Start+10))[,1])

  if(any(coordinates$Strand=="-")){
    coordinates$junction[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction[coordinates$Strand=="-"])))
  }

  #create rev comp junction
  coordinates$junction_revcomp <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction)))

  #top and bottom strand of oligo
  coordinates$TopStrand <- paste0("ctagcagt",coordinates$junction,"tagttatattcaagcata",coordinates$junction_revcomp,"gcg")
  coordinates$BotStrand <- paste0("aattcgc",coordinates$junction,"tatgcttgaatataacta",coordinates$junction_revcomp,"actg")

  if(shift){

    #create junction
    coordinates$junction5p <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End-6,coordinates$End))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start,coordinates$Start+13))[,1])

    #reverse anything that is negative strand
    if(any(coordinates$Strand=="-")){
      coordinates$junction5p[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction5p[coordinates$Strand=="-"])))
    }


    #create rev comp junction
    coordinates$junction_revcomp5p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction5p)))


    #top and bottom strand of oligo
    coordinates$TopStrand5p <- paste0("ctagcagt",coordinates$junction5p,"tagttatattcaagcata",coordinates$junction_revcomp5p,"gcg")
    coordinates$BotStrand5p <- paste0("aattcgc",coordinates$junction5p,"tatgcttgaatataacta",coordinates$junction_revcomp5p,"actg")


    #create junction
    coordinates$junction3p <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End-11,coordinates$End))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start,coordinates$Start+8))[,1])

    if(any(coordinates$Strand=="-")){
      coordinates$junction3p[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction3p[coordinates$Strand=="-"])))
    }


    #create rev comp junction
    coordinates$junction_revcomp3p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction3p)))

    #top and bottom strand of oligo
    coordinates$TopStrand3p <- paste0("ctagcagt",coordinates$junction3p,"tagttatattcaagcata",coordinates$junction_revcomp3p,"gcg")
    coordinates$BotStrand3p <- paste0("aattcgc",coordinates$junction3p,"tatgcttgaatataacta",coordinates$junction_revcomp3p,"actg")

  }

  if(writetab==T){
    write.table(coordinates,quote = F,sep = "\t",file = output,row.names = F)
  } else {
    return(coordinates) }

}
#circRNASpliceOligoDesigner(input_coordinates = "./test/circs_totest.txt")

#' Generate shRNAs against linear RNA splice junctions
#'
#'This script contains ONE funciton that will read a table with circRNA/exon-exon splicing junction coordinates and genenames in the format: Name Chr Start End.
#'Optional:a column named Strand
#'So far the following species are available: fly: dm3, dm6 and human: hs19, mice: mm10 and rat: rn4
#'
#'
#' @param sp character with Name of the genome to use, available are: dm6, dm3, hg19, hg38, mm10, rn4. Default is dm6
#' @param input_coordinates document with the coordinate of the circRNAs to generate shRNAs in TSV format tab seppatated. It requires a column with Chromosome, Start and End.
#' @param writetab a boolean, if writetab=T, then the function will not write any table, instead it will return the table to R
#' @param output if writetab = T, it will create the output in this path, default is "shRNAs.tsv"
#'
#'
#' @param shift boolean, if shift=T, return also 5´and 3´shift sequences
#'
#' @return if writetab=F, It generates a DataFrame with the original table with the shRNA deisng appended
#'
#'
#' @export
#'
#' @examples
#'
#' #circRNAoligoDesigner(input_coordinates = "../test/circs_totest.txt", writetab=F)
#'
#'
linRNASpliceOligoDesigner<-function(sp="dm6",input_coordinates="lin_totest.txt", writetab=T,output="shRNAs.tsv",shift=T){

  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(input_coordinates))
  )


  if (sp =="dm6") {
    require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm6.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm6)
    gen<-BSgenome.Dmelanogaster.UCSC.dm6

  }  else if (sp =="dm3") {
    require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm3.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm3)
    gen<-BSgenome.Dmelanogaster.UCSC.dm3

  } else if (sp =="hg19") {
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg19.knownGene
    require(BSgenome.Hsapiens.UCSC.hg19)
    gen<-BSgenome.Hsapiens.UCSC.hg19

  } else if (sp =="hg38") {
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg38.knownGene
    require(BSgenome.Hsapiens.UCSC.hg38)
    gen<-BSgenome.Hsapiens.UCSC.hg38

  } else if (sp =="mm10") {
    require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    tx<-TxDb.Mmusculus.UCSC.mm10.knownGene
    require(BSgenome.Mmusculus.UCSC.mm10)
    gen<-BSgenome.Mmusculus.UCSC.mm10

  }  else if (sp =="rn4") {
    require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    tx<-TxDb.Rnorvegicus.UCSC.rn4.ensGene
    require(BSgenome.Rnorvegicus.UCSC.rn4)
    gen<-BSgenome.Rnorvegicus.UCSC.rn4

  }

  #read the input data
  coordinates<-read.delim(input_coordinates)


  #create junction long
  coordinates$junction_long <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start-20,coordinates$Start))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End,coordinates$End+20))[,1])

  if(any(coordinates$Strand=="-")){
    coordinates$junction_long[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction_long[coordinates$Strand=="-"])))
  }

  #create junction
  coordinates$junction <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start-9,coordinates$Start))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End,coordinates$End+10))[,1])

  #reverse negative
  if(any(coordinates$Strand=="-")){
    coordinates$junction[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction[coordinates$Strand=="-"])))
  }

  #create rev comp junction
  coordinates$junction_revcomp <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction)))

  #top and bottom strand of oligo
  coordinates$TopStrand <- paste0("ctagcagt",coordinates$junction,"tagttatattcaagcata",coordinates$junction_revcomp,"gcg")
  coordinates$BotStrand <- paste0("aattcgc",coordinates$junction,"tatgcttgaatataacta",coordinates$junction_revcomp,"actg")

  if(shift){

    #create junction
    coordinates$junction5p <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start-6,coordinates$Start))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End,coordinates$End+13))[,1])

    #reverse anything that is negative strand
    if(any(coordinates$Strand=="-")){
      coordinates$junction5p[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction5p[coordinates$Strand=="-"])))
    }


    #create rev comp junction
    coordinates$junction_revcomp5p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction5p)))

    #top and bottom strand of oligo
    coordinates$TopStrand5p <- paste0("ctagcagt",coordinates$junction5p,"tagttatattcaagcata",coordinates$junction_revcomp5p,"gcg")
    coordinates$BotStrand5p <- paste0("aattcgc",coordinates$junction5p,"tatgcttgaatataacta",coordinates$junction_revcomp5p,"actg")


    #create junction
    coordinates$junction3p <- paste0(as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$Start-11,coordinates$Start))[,1],as.data.frame(Biostrings::getSeq(gen, coordinates$Chr,coordinates$End,coordinates$End+8))[,1])

    if(any(coordinates$Strand=="-")){
      coordinates$junction3p[coordinates$Strand=="-"] <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction3p[coordinates$Strand=="-"])))
    }

    #create rev comp junction
    coordinates$junction_revcomp3p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction3p)))

    #top and bottom strand of oligo
    coordinates$TopStrand3p <- paste0("ctagcagt",coordinates$junction3p,"tagttatattcaagcata",coordinates$junction_revcomp3p,"gcg")
    coordinates$BotStrand3p <- paste0("aattcgc",coordinates$junction3p,"tatgcttgaatataacta",coordinates$junction_revcomp3p,"actg")

  }

  if(writetab==T){
    write.table(coordinates,quote = F,sep = "\t",file = output,row.names = F)
  } else {
    return(coordinates) }

}
#linRNASpliceOligoDesigner(input_coordinates = "./test/lin_totest.txt")



