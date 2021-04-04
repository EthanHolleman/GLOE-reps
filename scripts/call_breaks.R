
# https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line
read_bed <- function(filepath){

    read_depths <- list()
    i <- 1
    con <- file(filepath, "r")
    while ( TRUE ) {
        line <- read_bed_line(con)
        if ( length(line) == 0 ) {
        break
        }
        read_depth[[i]] <- line$V5  # score which should be number aligned reads
    }
    close(con)
    }
    return(read_depths)
}


read_bed_line <- function(connection){

    line <- read.delim(con, nrows=1, header=F)
    return(line)

}

# Remove reads from bed file that have fewer than min_sds * sd of number
# of reads and write to a new file
trim_bed <- function(input_filepath, output_filepath, read_depths, min_sds=1){


    min_depth <- sd_depth * min_sds
    output_con <- file(output_filepath, "w")
    input_con <- file(input_filepath, "r")
    while (TRUE){

        line <- read_bed_line(input_con)
        if (length(line) > 0){
            if (line$V5 >= min_depth){
            write.table(output_con, line, header=F, 
                        quote=F, col.names=F, row.names=F
                        )
            }
        }else{
            break
        }

    }

}

main <- function(){

    args <- commandArgs(trailingOnly=TRUE)
    input.bed <- args[1]
    output.bed <- args[2]
    if (length(args) == 3){
        min_sds <- as.numberic(args[3])
    }else{
        min_sds <- 1
    }
    read_depths <- read_bed(input.bed)
    trim_bed <- trim_bed(input.bed, output.bed, read_depths, min_sds)
}

if (! interactive()){
    main()
}



