library(tidyr)
library(ggplot2)
library(ggpubr)
library(cluster)

read_meta_bed <- function(bed.path){
  
  bed.df <- read.table(
    bed.path, header=F, sep='\t'
  )
  colnames(bed.df) <- c(
    'chromosome', 'start', 'end', 'name', 'score', 'strand')
  bed.df <- separate(bed.df, 'name', c('name', 'index', 'bin_length'))
  bed.df$index <- as.numeric(bed.df$index)
  bed.df
}


metaplot <- function(bed.df){
  
  plt <- ggplot(bed.df, aes(x=index, y=score)) + 
    geom_point(size=1, alpha=0.5, color='dodgerblue')  + labs(x='Percent gene body') +
    ylim(min(bed.df$score), max(bed.df$score)) + theme_pubr()
  plt
}

clustering <- function(bed.df){
  # cluster genes by there nicks
  



}


main <- function(){
  
  args <- commandArgs(trailingOnly=TRUE)
  bed.path <- args[1]
  output.path <- args[2]
  bed.df <- read_meta_bed(bed.path)
  p <- metaplot(bed.df)
  ggsave(output.path, p, device = 'png', dpi=500)
  
}

if (! interactive()){
    main()
}