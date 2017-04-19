#' ---
#' title: "chromVAR analysis"
#' author: "Put your name here idk"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#+ echo=FALSE, eval=FALSE
# Run these commands only once ever!
install.packages("devtools")
install.packages("plotly")
install.packages("heatmaply")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("tidyverse")

devtools::install_github("GreenleafLab/chromVAR")
devtools::install_github("caleblareau/BuenColors")

source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")
biocLite("SummarizedExperiment")
bioclite("BSgenome.Hsapiens.UCSC.hg19")

#' ### The goal of this document is to produce several interactive plots for ATAC data
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
# Import libraries
library(plotly)
library(heatmaply)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)
library(BuenColors)

#' ## Initialize parallel processing
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
library(BiocParallel)
register(MulticoreParam(2)) # adjust according to your machine

#' ## Import/Filter data
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
peakfile <- "../data/peaks/heme_peaks.bed"
peaks <- get_peaks(peakfile)

bamfiles <- list.files("../data/bams", full.names = TRUE)
raw_counts <- get_counts(bamfiles, peaks, paired =  TRUE, by_rg = TRUE, format = "bam",
                              colData = DataFrame(source = bamfiles))
# discuss paired, by_rg, etc.
counts_filtered <- filter_samples(raw_counts, min_depth = 500,
                                  min_in_peaks = 0.15, shiny = FALSE)
counts <- filter_peaks(counts_filtered)

#' ## See effect of filtering
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE
dim(raw_counts)
dim(counts_filtered)
dim(counts)

#' ## Get GC content/peak; get motifs from JASPAR/kmer; find overlaps in peaks
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
counts <- add_gc_bias(counts, genome = BSgenome.Hsapiens.UCSC.hg19)
motifs <- get_jaspar_motifs()
motif_ix <- match_motifs(motifs, counts, genome = BSgenome.Hsapiens.UCSC.hg19)
kmer_ix <- match_kmers(7, counts)

#' ## Compute deviations; typically time consuming
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
dev <- compute_deviations(object = counts, annotations = motif_ix)

#' ## Find variable motifs
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, fig.align='center', fig.cap = "Figure 1; variability across motifs"
variability <- compute_variability(dev)
plot_variability(variability, use_plotly = TRUE) 

#' ## Examine variable motifs in sample
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, fig.align='center', fig.cap = "Figure 2; motifs x samples"
mostvariable <- tail(sort(variability$variability, index.return = TRUE)$ix,30)
m <- assay(dev)[mostvariable,]
heatmaply(m, colors = jdb_palette("flame_light"))
