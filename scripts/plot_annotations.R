library(ggplot2)
library(ggpubr)


main <- function(){

    args <- commandArgs(trailingOnly=TRUE)

    if ( length(args) == 3 ){
        # passing in multible samples
        # pattern should be normal, reverse (peak calling) output
        normal.tsv <- args[1]
        reverse.tsv <- args[2]
        normal.df <- read.delim(normal.tsv, sep='\t')
        reverse.df <- read.delim(reverse.tsv, sep='\t')
        normal.df$identity <- rep('normal', nrow(normal.df))
        reverse.df$identity <- rep('reverse', nrow(reverse.df))
        freq.df <- rbind(normal.df, reverse.df)
        
    }else{
        freq.df <- read.delim(args[1], sep='\t')
        freq.df$identity <- rep('sample', nrow(freq.df ))
        
    }
    output.path <- args[length(args)]
    print(colnames(freq.df))
    barplot_freq <- ggplot(freq.df, aes(fill=Feature, x=identity, y=Frequency)) +
            geom_bar(position='dodge', stat='identity') +
            labs(y='Proportion', x='Peak calling method', fill='Annotation')
    barplot_num_peaks <- ggplot(freq.df, aes(x=identity, y=num_peaks)) +
                         geom_bar(stat='identity') + 
                         labs(y='Number of Peaks', x='Peak calling method') + 
                         theme_pubr()
    
    
    
    full_plt <- ggarrange(barplot_freq, barplot_num_peaks, nrow=1, ncol=2,
                        widths = c(2, 1), labels=c('A', 'B'))
    
    ggsave(output.path, full_plt, dpi=500, unit='in', width=10, height=7)

}


if (!interactive()){
    main()
}