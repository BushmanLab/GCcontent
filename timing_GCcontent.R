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

library(intSiteRetriever)
library(hiAnnotator)


reference <- get_reference_genome("hg19")

time_GC <- function() {
    sites <- lapply(c(10, 100, 1000, 10000, 100000), function(num) { 
        site <- get_random_positions(seq(1, num), reference, 'f', number_of_positions=1)
        site <- makeGRanges(site, soloStart=TRUE,
                    chromCol='chr', strandCol='strand', startCol='position')
        site
    })
    sapply(sites, function(site) { 
        window_size <- c(large=100) 
        GC_time <- system.time(
            getGCpercentage(site, "GC", window_size, reference))['user.self']
        list(length(site), GC_time)
    })
}
