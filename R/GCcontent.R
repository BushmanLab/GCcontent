#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' calulate GC for several window sizes
#'
#' @param sites GRange locations on genome: chromosome, strand, position
#' @param column_prefix what prefix to add to each column
#' @param window_size vector with names: elements are window, names used as col names
#' @param reference_genome_sequence actual reference genome (from BSgenome.* object)
#' to get reference_genome_sequence use get_reference_genome() function from intSiteRetiever
#'
#' @return Grange object with added columns: prefix.name_window_size 
#' and value equal to GC content for a window
#' @export
#'
#' @note start up time to load human genome is about 10 second
getGCpercentage <- function(
    sites, column_prefix, window_size, reference_genome_sequence
){
    stopifnot(length(window_size) == length(names(window_size)))
    metadata <- mcols(sites)

    rangesToCalc <- .expand_trim_GRanges(sites, window_size)

   #seqs will take a lot of memory
    #could split at severe cpu time penelty
    seqs <- getSeq(reference_genome_sequence, rangesToCalc, as.character=F)

    letterFreqs <- letterFrequency(seqs, c("G", "C", "A", "T"))
    rm(seqs)

    GC <- letterFreqs[, c("G", "C")]
    ATGC <- letterFreqs[, c("A", "T", "G", "C")]

    gcContent <- rowSums(GC)/rowSums(ATGC)

    gcContent[!is.finite(gcContent)] <- NA #handled gracefully by pipeUtils

    gcContent <- DataFrame(matrix(gcContent, nrow=length(sites)))

    names(gcContent) <- paste(column_prefix, names(window_size), sep=".")

    mcols(sites) <- cbind(metadata, gcContent)

    sites
}

.expand_trim_GRanges <- function(sites, window_size) {
    nsites <- length(sites)
    strand(sites) = "+" #unimportant for GC and speeds up later calculations

    sites.seqinfo.original <- seqinfo(sites)
    isCircular(seqinfo(sites)) <- rep(FALSE, length(seqinfo(sites)))

    sites <- rep(sites, length(window_size))
    sites <- trim(suppressWarnings(flank(sites,
                                       rep(window_size/2, each=nsites),
                                       both=T)))

    seqinfo(sites) = sites.seqinfo.original
    sites
}
