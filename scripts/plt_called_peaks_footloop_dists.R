library(ggpubr)
library(ggplot2)

read_bed <- function(bed.path){
  
  bed.df <- as.data.frame(read.table(
    bed.path, header=F
  ))
  bed.df
}


plot_distance_distrabution <- function(bed.df.fwd, bed.df.rev, xlim_bottom=0, 
                                       xlim_top=0, x.lab=''){
  
  big.bed.df <- rbind(bed.df.fwd, bed.df.rev)
  plt <- ggplot(big.bed.df, aes(x=`V12`, fill=V6)) + geom_density(alpha=0.7) + 
    labs(x=x.lab, y='Density', fill='Strand') + 
    theme_pubr() + scale_fill_manual(values = c('#eb3434', '#3493eb'))
  if (xlim_bottom != 0 & xlim_top != 0){
    plt <- plt + xlim(xlim_bottom, xlim_top)
  }
  plt
    
}

main <- function(){
  
  args = commandArgs(trailingOnly=TRUE)
  fwd.bed <- args[1]
  rev.bed <- args[2]
  output.path <- args[3]
  
  fwd.bed.df <- read_bed(fwd.bed)
  rev.bed.df <- read_bed(rev.bed)
  plt.full <- plot_distance_distrabution(fwd.bed.df, rev.bed.df)
  plt.zoom <- plot_distance_distrabution(fwd.bed.df, rev.bed.df, 
                                         xlim_bottom=-100, xlim_top=100, )
  plt.all <- ggpubr::ggarrange(plt.full, plt.zoom, nrow = 2, ncol=1, 
                               labels=c('A', 'B'), common.legend = TRUE)
  ggsave(output.path, plt.all, dpi=500)
  
  
}

if (! interactive()){
    main()
}