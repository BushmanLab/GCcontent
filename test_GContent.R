source("intSiteRetriever/random_site.R")

context("Testing GC perecntage calculation and edge cases")

reference <- get_reference_genome("hg18")
human_seqinfo <- seqinfo(reference)

sites <- GRanges(
    seqnames=Rle(c("chr1", "chr11", "chr11")),
    ranges=IRanges(start=c(100, 200, 50000), width=rep(1, 3)),
    strand=Rle(c("+", "-", "+")),
    seqinfo=human_seqinfo
)

window_size <- c(small=10, large=1000)
sites_GC <- getGCpercentage(sites, "GC", window_size, reference)
GC_small <- sites_GC$GC.small
GC_large <- sites_GC$GC.large

test_that("can calculate GC for case with all bases known and valid window", {
    expect_true(GC_small[1] > 0)
    expect_true(GC_small[1] < 1)
    expect_true(GC_small[3] > 0)
    expect_true(GC_small[3] < 1)
})

test_that("single window works", {
    window_size_singleton <- c(small=10)
    getGCpercentage(sites, "GC", window_size_singleton, reference)
})

test_that("if all Ns need NA", {
    expect_true(is.na(GC_small[2]))
    expect_true(is.na(GC_large[2]))
})

test_that("can calculate GC for case with window larger than genome start/end", {
    expect_true(GC_large[1] > 0)
    expect_true(GC_large[1] < 1)
})

test_that("Ns are ignored", {
    expect_true(GC_small[3] > 0.1)
    expect_true(GC_small[3] < 0.3)
})

