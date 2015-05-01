library(BSgenome)
library(Repitools)

#' calulate GC for several window sizes
#'
#' @param reference_genome_sequence actual reference genome (from BSgenome)
#' to get reference_genome_sequence use get_reference_genome() function
#'
#' @note start up time to load human genome is about 10 second
getGCpercentage <- function(
    sites, column_prefix, window_size, reference_genome_sequence
) {
    stopifnot(length(window_size) == length(names(window_size)))
    metadata <- mcols(sites)
    chrom_position <- convert_GRange_to_chr_pos(sites)
    sapply(seq(window_size), function(i) {
        val <- window_size[i]
        name <- names(window_size)[i]
        column_name <- paste0(column_prefix, ".", name)
        metadata[[column_name]] <<- get_gc_percentage_for_single_window(
            chrom_position, val, reference_genome_sequence)
    })
    mcols(sites) <- metadata
    sites
}

#' return genome seq for human readable UCSC format
#' 
#' format is: hg18, ...
get_reference_genome <- function(reference_genome) {
    pattern <- paste0("\\.", reference_genome, "$")
    match_index <- which(grepl(pattern, installed.genomes()))
    stopifnot(length(match_index) == 1)
    BS_genome_full_name <- installed.genomes()[match_index]
    get(BS_genome_full_name)
}

get_gc_percentage_for_single_window <- function(chrom_position, window, ref_genome_seq) {
    stopifnot(length(window) == 1)
    gcContentCalc(chrom_position, ref_genome_seq, window)
}

convert_GRange_to_chr_pos <- function(sites) {
    data.frame(
        "chr"=as.character(seqnames(sites)), 
        "position"=as.numeric(start(sites)),
        stringsAsFactors=FALSE
    )
}
