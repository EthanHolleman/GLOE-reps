options(scipen=999)

# https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line
read_bed <- function(filepath){
  bed.df <- as.data.frame(read.delim(filepath, header = F))
  colnames(bed.df) <- c('chr', 'end_position', 'depth')
  bed.df
}

# Remove reads from bed file that have fewer than min_sds * sd of number
# of reads and write to a new file
trim_coverage <- function(read_depths.df){
  
  sd_depth <- sd(read_depths.df$depth)
  min_depth <- median(read_depths.df$depth) - median(read_depths.df$depth) * 0.75
  # if (min_depth < 100){
  #   min_depth <- 100
  # }
  subset(read_depths.df, depth >= min_depth)
  
}


convert_to_bed <- function(read_depths.df, strand){
  # the coverage file will only have the base position. In bed file this is
  # the end position so need to readd the start position
  # add a name as well
  message(0)
  print(dim(read_depths.df))
  read_depths.df$name <- paste(rep('depth', nrow(read_depths.df)), 1:nrow(read_depths.df), sep='_')
  message(1)
  read_depths.df$start_position <- read_depths.df$end_position - 1
  message(2)
  read_depths.df$strand <- rep(strand, nrow(read_depths.df))
  message(3)
  read_depths.df <- read_depths.df[, c('chr', 'start_position', 'end_position', 'name', 'depth', 'strand')]
  message(4)
  read_depths.df
  
  
}

remove_scaffolds <- function(read_depths.df){
  names <- c(1:22, 'X', 'Y')
  chrs <- paste(rep('chr', length(names)), names, sep='')
  read_depths.df <- subset(read_depths.df, read_depths.df$chr %in% chrs)
  read_depths.df
}


main <- function(){
  
  args <- commandArgs(trailingOnly=TRUE)
  input.bed <- args[1]
  output.bed <- args[2]
  strand <- args[3]
  message(input.bed)
  read_depths.df <- read_bed(input.bed)
  message('Read file')
  read_depths.df.trim <- trim_coverage(read_depths.df)
  message('Trimmed coverage')
  read_depths.df.trim <- convert_to_bed(read_depths.df.trim, strand)
  message('Converted to bed')
  read_depths.df.trim <- remove_scaffolds(read_depths.df.trim)
  message('Removed scaffolds')
  write.table(read_depths.df.trim, output.bed, col.names = F, 
              row.names = F, quote = F, sep='\t')
  message('Wrote bed file')
  message('Done')
  
}

if (! interactive()){
  main()
}



