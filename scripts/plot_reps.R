library(ggplot2)
library(tidyr)
library(ggpubr)
library(RColorBrewer)


# now we are also having to deal with forward and reverse strands as well


extract_replicates_from_filepath <- function(file.path){
  
  file.path.split <- unlist(strsplit(file.path, "/"))
  samples <- unlist(strsplit(file.path.split, "-"))
  return(samples)
  
}


get_modification_from_samples <- function(samples, sample.df){
  
  sample_1 <- samples[1]
  sample_2 <- samples[2]
  
  sample_1.mod <- sample.df[sample.df$Sample.Name == sample_1, ]
  sample_2.mod <- sample.df[sample.df$Sample.Name == sample_2, ]
  
  stopifnot(sample_1.mod == sample_2.mod)
  
  return(sample_1.mod)
  
}


read_sample_csv <- function(sample.csv){
  
  samples.df <- as.data.frame(read.csv(sample.csv))
  return(samples.df)
  
}


read_bed <- function(filepath){
  
  closest_nick <- list()
  i <- 1
  con <- file(filepath, "r")
  while ( TRUE ) {
    line <- read_bed_line(con)
    if ( length(line) == 0 ) {
      break
    }
    closest_nick[[i]] <- line$V5  # score which should be number aligned reads
  }
  close(con)
  return(closest_nick)
}


add_mod_to_closest_nick <- function(closest.nick.list, mod){
  
  df <- data.frame(
    'closest_nick'=unlist(closest.nick.list),
    'modification'=rep(mod, length(closest.nick.list))
  )
  return(df)
  
}



read_close_bed <- function(bed.path, samples.df){
  
  samples <- extract_replicates_from_filepath(bed.path)
  mod <- get_modification_from_samples(samples)
  closest.nick.list <- read_bed(bed.path)
  closest.nick.df <- add_mod_to_closest_nick(closest.nick.list, mod)
  return(closest.nick.df)
  
}

read_all_beds <- function(bed.path.list, samples.path){
  
  samples.df <- read_sample_csv(samples.path)
  bed.df.list <- list()
  for (bed.path in bed.path.list){
    bed.df.list[[bed.path]] <- read_close_bed(
      bed.path, samples.df
    )
  }
  
  return(bed.df.list)
  
}


concatenate_bed_dfs <- function(bed.df.list){
  
  do.call('rbind', bed.df.list)
  
}


plot_nick_dist_by_mod <- function(all.samples.df){
  
  ggplot(all.samples.df, 
         aes(x=modification, y=closest_nick, fill=modification)) +
    scale_fill_brewer(palette="Dark2") + 
    labs(x='Treatment', y='Closest replicate nick') +
    theme_pubclean()
    
}


main <- function(){
  
  args <- commandArgs(trailingOnly=TRUE)
  # first arg is the samples.csv
  samples.path <- args[1]
  # next is the output path
  output.path <- args[2]
  # all further arguments are replicate comparison bed files
  bed.path.list <- args[3:length(args)]
  
  bed.df.list <- read_all_beds(bed.path.list, samples.path)
  all.samples.df <- concatenate_bed_dfs(bed.df.list)
  plt1 <- plot_nick_dist_by_mod(all.samples.df)
  ggsave(output.path, plt1, dpi=500)
  
  
}


if (! interactive()){
  main()
}






# boxplot of distances between breaks between replicates
# want to do this for all replicates so will need to read them
# all into a file
# read all files and concatenate into one large one. Need to add
# a row that specifies the replicate comparison
# want another script that looks at signal strength as well number
# of reads for givben break 
