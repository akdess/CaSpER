#source("https://bioconductor.org/biocLite.R")
#biocLite(c('HMMcopy', 'GenomeGraphs', 'biomaRt', 'limma', 'GO.db', 'org.Hs.eg.db', 'GOstats'))

install.packages('devtools');

require(devtools)
install_github("akdess/CaSpER")
