
library(ggplot2)
library(ggpubr)

read_bed <- function(bed.path){
  
  bed.df <- as.data.frame(read.table(
    bed.path, header=F
  ))
  bed.df <- remove_scaffolds(bed.df)
  bed.df
}


read_sample_csv <- function(sample.csv){
  
  samples.df <- as.data.frame(read.csv(sample.csv))
  return(samples.df)
  
}


extract_replicates_from_filepath <- function(file.path){
  
  file.path.split <- unlist(strsplit(file.path, "/"))[3]
  samples <- unlist(strsplit(file.path.split, "-"))
  return(samples)
  
}


get_modification_from_replicates <- function(samples, sample.df){
  
  sample_1 <- samples[1]
  sample_2 <- samples[2]
  
  sample_1.mod <- sample.df[sample.df$Sample.Name == sample_1, ]$modification
  sample_2.mod <- sample.df[sample.df$Sample.Name == sample_2, ]$modification
  
  print(sample_1.mod)
  print(sample_2.mod)
  stopifnot(sample_1.mod == sample_2.mod)
  
  return(sample_1.mod)
  
}

remove_scaffolds <- function(read_depths.df){
  names <- c(1:22, 'X', 'Y')
  chrs <- paste(rep('chr', length(names)), names, sep='')
  subset(read_depths.df, V1 %in% chrs)

}


make_boxplot <- function(bed.df){
  
  plt <- ggplot(bed.df, aes(x=`V13`)) + geom_boxplot() + 
    theme_pubr() + labs(y='', x='Distance to closest replicate nick') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  x.coord <- layer_scales(plt)$x$range$range[2] / 2
  #plt <- annotate("text", x=x.coord, y=mean(bed.df$V13) + 10, label=mean(bed.df$V13))
  plt
}

make_scatter <- function(bed.df, mod){

  ggplot(bed.df, aes(x=`V2`, y=`V8`, color=`V1`), size=15) + geom_point() +
    geom_abline(slope=1, intercept=0) + 
    facet_wrap(~V1) + theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x='', y='', title=mod) + theme(legend.position = "none") 

}

make_table <- function(bed.df){

  stable <- desc_statby(bed.df, measure.var = "V13", grps="V1")

  print(colnames(stable))
  stable <- stable[, c("V1", "median", "mean", "sd")]
  stable.p <- ggtexttable(stable, rows=NULL)
  stable.p



}



main <- function(){
  args <- commandArgs(trailingOnly=TRUE)
  bed.path <- args[1]
  samples.path <- args[2]
  output.path <- args[3]

  samples.df <- read_sample_csv(samples.path)
  message('Reading bed file')
  bed.df <- read_bed(bed.path)
  print(unique(bed.df$V1))
  replicates <- extract_replicates_from_filepath(bed.path)
  print(replicates)
  mod <- get_modification_from_replicates(replicates, samples.df)
  scatter.plt <- make_scatter(bed.df, as.character(mod))
  boxplot.plt <- make_boxplot(bed.df)
  table.plt <- make_table(bed.df)
  main.plt <- ggarrange(
                        ggarrange(scatter.plt, table.plt, nrow=1, ncol=2),
                        boxplot.plt, nrow=2, ncol=1, 
                        heights=c(2, 0.5), labels=c('A', 'B')
                        )
        
  message("Made boxplot")
  print(output.path)
  ggsave(output.path, main.plt, device='png', dpi=500, units='in', height=12, width=12)

}

if (! interactive()){
  main()
}








