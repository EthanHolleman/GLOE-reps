library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

main <- function(){

    # read the peak (bed file) from peak calling
    args = commandArgs(trailingOnly=TRUE)

    input.bed <- args[1]
    output.file <- args[2]
    files <- getSampleFiles()
    #peaks <- readPeakFile(files[[4]])
    peaks <- readPeakFile(args[1])
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    # annotate peaks
    peak.anno <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb,
                              annoDb='org.Hs.eg.db',
                              genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic")
                              )

    
    # add nubmber of peaks as extra column 
    anno.stat <- attributes(peak.anno)$annoStat
    num.peaks <-  attributes(peak.anno)$peakNum
    anno.stat$num_peaks <- rep(num.peaks, nrow(anno.stat))
    write.table(anno.stat, output.file, sep='\t')

}

if (! interactive()) {
    main()
}
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# main <- function(){
  
#   # read the peak (bed file) from peak calling
#   args = commandArgs(trailingOnly=TRUE)
  
#   input.bed <- args[1]
#   output.file <- 'test.anno.tsv'
#   files <- getSampleFiles()
  
#   #peaks <- readPeakFile(args[1])
#   peaks <- readPeakFile(files[[4]])
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
#   # annotate peaks
#   peak.anno <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb,
#                             annoDb='org.Hs.eg.db')
#   anno.stat <- attributes(peak.anno)$annoStat
#   num.peaks <-  attributes(peak.anno)$peakNum
#   anno.stat$num_peaks <- rep(num.peaks, nrow(anno.stat))
#   write.table(anno.stat, output.file)
  
# }


# main()