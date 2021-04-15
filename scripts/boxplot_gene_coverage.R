library(ggplot2)
library(ggpubr)

read_intersect_bed_with_read_depth <- function(bed.path){

    con <- file(bed.path, 'r')
    read_depth <- as.numeric(unlist(strsplit(readlines(con, n=1), ' '))[1])
    sample_name <- readlines(con, n=1)
    group <- readlines(con, n=1)
    bed.df <- as.data.frame(read.table(bed.path, header=F))
    colnames(bed.df) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'count')
    bed.df$read_depth <- rep(read_depth, now(bed.df))
    bed.df$sample_name <- rep(sample_name, now(bed.df))
    bed.df$group <- rep(group, now(bed.df))
    close(con)
    bed.df

}

standardize_intersections <- function(bed.df){

    bed.df$stand_intersects <- scale(bed.df$count)
    bed.df

}

normalize_intersections_to_read_count <- function(bed.df, read_depth){

    bed.df$intersect_over_read_count <- bed.df$count / read_depth
    bed.df

}

reads_col_plot <- function(bed.df){

    temp.df <- cbind(bed.df$read_depth, bed.df$group, bed.df$sample_name)
    temp.df <- unique(temp.df)
    plt <- ggplot(temp.df, aes(x=sample_name, y=read_depth, fill=group)) +
           geom_col(color='black') + scale_fill_brewer(palette="Dark2") +
           coord_flip() + theme_pubr()
    plt
}

boxplot_intersections <- function(bed.df){

    plt <- ggplot(bed.df, aes(x=intersect_over_read_count, y=sample_name)) +
            geom_boxplot() + facet_wrap(~group) + theme_pubr()
    plt

}

scatter_plot <- function(bed.df){

    # for each mean do something
    # split by groups, average intersections for each gene
    # seperate by group
    groups <- unique(bed.df$groups)
    # should only be two groups
    bed.df.group_1 <- aggregate(subset(bed.df, group==groups[1])[, 'intersect_over_read_count', 'stand_intersects'], list(bed.df.group$name), mean)
    bed.df.group_2 <- aggregate(subset(bed.df, group==groups[2])[, 'intersect_over_read_count', 'stand_intersects'], list(bed.df.group$name), mean)
    # genes should then be the sampe in both
    # plot the corrindates of each
    
    
    



}