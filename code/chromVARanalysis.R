#' ---
#' title: "chromVAR analysis"
#' author: "Caleb Lareau"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

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
library(chromVARxx)
library(chromVARmotifs)

#' ## Initialize parallel processing
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
library(BiocParallel)
register(MulticoreParam(2)) # adjust according to your machine

#' ## Import/Filter data
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
peakfile <- "../data/peaks/heme_peaks.bed"
peaks <- getPeaks(peakfile)

bamfiles <- list.files("../data/bams", full.names = TRUE)
raw_counts <- getCounts(bamfiles, peaks, paired =  TRUE, by_rg = TRUE, format = "bam",
                              colData = DataFrame(source = bamfiles))
# discuss paired, by_rg, etc.
counts_filtered <- filterSamples(raw_counts, min_depth = 500,
                                  min_in_peaks = 0.15, shiny = FALSE)
counts <- filterPeaks(counts_filtered)

#' ## See effect of filtering
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE
dim(raw_counts)
dim(counts_filtered)
dim(counts)

#' ## Get GC content/peak; get motifs from chromVARmotifs package; find kmers
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
counts <- addGCBias(counts, genome = BSgenome.Hsapiens.UCSC.hg19)
data("human_pwms_v1") # also mouse_pwms_v1
motif_ix <- matchMotifs(human_pwms_v1, counts, genome = BSgenome.Hsapiens.UCSC.hg19)
kmer_ix <- matchKmers(5, counts, genome = BSgenome.Hsapiens.UCSC.hg19)

dim(motif_ix)
dim(kmer_ix)

#' ## Compute deviations; typically time consuming
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
dev <- computeDeviations(object = counts, annotations = motif_ix)

#' ## Find variable motifs
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, fig.align='center', fig.cap = "Figure 1a; variability across motifs"
variabilityAll <- computeVariability(dev)
plotVariability(variabilityAll, use_plotly = TRUE)

#' ## Bag related deviation scores
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE,
bagged <- bagDeviations(dev, 0.7, "human")
dim(dev)
dim(bagged)

#' ## Find variable bagged motifs
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, fig.align='center', fig.cap = "Figure 1b; variability across reasonably unique motifs"
variabilityBagged <- computeVariability(bagged)
plotVariability(variabilityBagged, use_plotly = TRUE)

#' ## Examine variable motifs in sample
#+ cache = FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10, echo = TRUE, fig.align='center', fig.cap = "Figure 2a; motifs x samples"
mostvariable <- tail(sort(variabilityAll$variability, index.return = TRUE)$ix,30)
m <- assay(dev)[mostvariable,]
heatmaply(m, colors = jdb_palette("flame_light"))

#' ## Examine variable bagged motifs in sample
#+ cache = FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10,  echo = TRUE, fig.align='center', fig.cap = "Figure 2b; bagged motifs x samples"
mostvariable <- tail(sort(variabilityBagged$variability, index.return = TRUE)$ix,30)
m <- assay(bagged)[mostvariable,]
heatmaply(m, colors = jdb_palette("flame_light"))
