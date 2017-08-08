
# Run these commands only once to make installation!
install.packages("devtools")
install.packages("plotly")
install.packages("heatmaply")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("tidyverse")

source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")
biocLite("SummarizedExperiment")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

devtools::install_github("GreenleafLab/chromVAR")
devtools::install_github("GreenleafLab/motifmatchr")
devtools::install_github("GreenleafLab/chromVARmotifs")

devtools::install_github("caleblareau/BuenColors")
devtools::install_github("caleblareau/chromVARxx")
