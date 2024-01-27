install.packages("renv")

library(renv)
renv::init()

install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("doParallel")
install.packages('tidyverse')
install.packages('gsubfn')
install.packages('zoo')
install.packages('snow')
install.packages('plyr')
install.packages('gtools')
install.packages('ggsci')
install.packages('igraph')
install.packages('tidygraph')
install.packages('RColorBrewer')
install.packages('parallel')
install.packages('foreach')
install.packages("stringdist")
#install.packages('scatterpie') <- Commented out because of an error.

# run the fdllowing at the Terminal: 
# sudo apt install build-essential
# sudo apt install liblapack-dev libopenblas-dev
# sudo apt-get install gfortran 


install.packages('./rELA.v0.44.tar.gz')

# sudo apt-get install libglpk-dev # for installing igraph in RStudio

#renv::snapshot()

install.packages('devtools')
install.packages('ggtext')
devtools::install_github("hiroakif93/R-functions/AnalysisHelper", force=TRUE)

renv::snapshot()

