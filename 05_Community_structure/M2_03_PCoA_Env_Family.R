############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### 2023.12.21 Toju
#### R 4.3.2
#### 
############################################################################

set.seed(123)
setwd("../")


## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/M2_03_output')

# -- Load data table
sml <- readRDS("Table/sml_1474.rds") 
env <- readRDS("Table/env_1474.rds") 

dimrm <- readRDS('Table/PCoA_Family_rrarefy.rds')
pca <- list(dimrm[[1]]$points, dimrm[[2]]$points)
pca <- lapply(pca, function(x){colnames(x) <- paste('PcoA', 1:ncol(x)); return(x)})
pca_prk <- data.frame(SampleID=rownames(pca[[1]]), pca[[1]])
pca_fng <- data.frame(SampleID=rownames(pca[[2]]), pca[[2]])

prk <- merge(sml, pca_prk, by="SampleID")
fng <- merge(sml, pca_fng, by="SampleID")


############################################################################
# environmental gradient vs. PCoA axes
# Prokaryotes

glist <- list()

glist[[1]] <- ggplot(prk, aes(x=latitude_north, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(prk, aes(x=latitude_north, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x= pH_dry_soil, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(prk, aes(x= pH_dry_soil, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[5]] <- ggplot(prk, aes(x= log10(EC_electric_conductivity), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[6]] <- ggplot(prk, aes(x= log10(EC_electric_conductivity), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[7]] <- ggplot(prk, aes(x= CNratio, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[8]] <- ggplot(prk, aes(x= CNratio, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[9]] <- ggplot(prk, aes(x= log10(available_P), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[10]] <- ggplot(prk, aes(x= log10(available_P), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 15,
           filename=sprintf('%s/EnvGradient_Family_prk.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes
# Fungi

glist <- list()

glist[[1]] <- ggplot(fng, aes(x=latitude_north, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x=latitude_north, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(fng, aes(x= pH_dry_soil, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x= pH_dry_soil, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[5]] <- ggplot(fng, aes(x= log10(EC_electric_conductivity), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[6]] <- ggplot(fng, aes(x= log10(EC_electric_conductivity), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[7]] <- ggplot(fng, aes(x= CNratio, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[8]] <- ggplot(fng, aes(x= CNratio, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[9]] <- ggplot(fng, aes(x= log10(available_P), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[10]] <- ggplot(fng, aes(x= log10(available_P), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 15,
           filename=sprintf('%s/EnvGradient_Family_fng.pdf', dir$figdir))
           
############################################################################
# environmental gradient vs. PCoA axes: without latitude
# Prokaryotes

glist <- list()

glist[[1]] <- ggplot(prk, aes(x= pH_dry_soil, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(prk, aes(x= pH_dry_soil, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x= log10(EC_electric_conductivity), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(prk, aes(x= log10(EC_electric_conductivity), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[5]] <- ggplot(prk, aes(x= CNratio, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[6]] <- ggplot(prk, aes(x= CNratio, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[7]] <- ggplot(prk, aes(x= log10(available_P), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[8]] <- ggplot(prk, aes(x= log10(available_P), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 14,
           filename=sprintf('%s/EnvGradient_Family_prk.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes: without latitude
# Fungi

glist <- list()

glist[[1]] <- ggplot(fng, aes(x= pH_dry_soil, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x= pH_dry_soil, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(fng, aes(x= log10(EC_electric_conductivity), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x= log10(EC_electric_conductivity), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[5]] <- ggplot(fng, aes(x= CNratio, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[6]] <- ggplot(fng, aes(x= CNratio, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[7]] <- ggplot(fng, aes(x= log10(available_P), y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[8]] <- ggplot(fng, aes(x= log10(available_P), y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 14,
           filename=sprintf('%s/EnvGradient_Family_fng.pdf', dir$figdir))


############################################################################
# environmental gradient vs. PCoA axes
# Prokaryotes

glist <- list()

glist[[1]] <- ggplot(prk, aes(x = pH_dry_soil, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[2]] <- ggplot(prk, aes(x = pH_dry_soil, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x = log10(EC_electric_conductivity), y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[4]] <- ggplot(prk, aes(x = log10(EC_electric_conductivity), y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[5]] <- ggplot(prk, aes(x = CNratio, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[6]] <- ggplot(prk, aes(x = CNratio, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[7]] <- ggplot(prk, aes(x = log10(available_P), y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[8]] <- ggplot(prk, aes(x = log10(available_P), y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 14,
           filename=sprintf('%s/Contour_EnvGradient_Family_prk.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes
# Fungi

glist <- list()

glist[[1]] <- ggplot(fng, aes(x = pH_dry_soil, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x = pH_dry_soil, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[3]] <- ggplot(fng, aes(x = log10(EC_electric_conductivity), y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x = log10(EC_electric_conductivity), y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[5]] <- ggplot(fng, aes(x = CNratio, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[6]] <- ggplot(fng, aes(x = CNratio, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[7]] <- ggplot(fng, aes(x = log10(available_P), y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[8]] <- ggplot(fng, aes(x = log10(available_P), y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 14,
           filename=sprintf('%s/Contour_EnvGradient_Family_fng.pdf', dir$figdir))
           