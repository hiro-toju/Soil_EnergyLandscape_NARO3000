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

dimrm <- readRDS('Table/PCoA_Family_rrarefy.rds')
pca <- list(dimrm[[1]]$points, dimrm[[2]]$points)
pca <- lapply(pca, function(x){colnames(x) <- paste('PcoA', 1:ncol(x)); return(x)})
pca_prk <- data.frame(SampleID=rownames(pca[[1]]), pca[[1]])
pca_fng <- data.frame(SampleID=rownames(pca[[2]]), pca[[2]])

############################################################################
env_p <- readRDS("Table/env_pro_1771.rds") 
env_f <- readRDS("Table/env_fun_1664.rds") 

envcomp_p <- prcomp(env_p, scale.=TRUE, rank.=2)
envcomp_f <- prcomp(env_f, scale.=TRUE, rank.=2)

sink(sprintf('%s/summary_envcomp_pro1771.txt', dir$tabledir))
summary(envcomp_p)
cor(envcomp_p$x, env_p)
envcomp_p$rotation
sink()

sink(sprintf('%s/summary_envcomp_fun1664.txt', dir$tabledir))
summary(envcomp_f)
cor(envcomp_f$x, env_f)
envcomp_f$rotation
sink()

envPC_p <- envcomp_p$x
colnames(envPC_p) <- c('environment_PC1', 'environment_PC2')
envPC_p <- data.frame(SampleID=rownames(envPC_p), envPC_p)

envPC_f <- envcomp_f$x
colnames(envPC_f) <- c('environment_PC1', 'environment_PC2')
envPC_f <- data.frame(SampleID=rownames(envPC_f), envPC_f)

prk <- merge(envPC_p, pca_prk, by="SampleID")
fng <- merge(envPC_f, pca_fng, by="SampleID")

############################################################################
color.ramp.length <- 16
negative.length <- 8
positive.length <- 8
cols <- c(colorRampPalette(c("blue", "white"))(negative.length),
          colorRampPalette(c("white", "red"))(positive.length))

lev_p <- cor(envcomp_p$x, env_p)
lev_f <- cor(envcomp_f$x, env_f)

glist <- list()

glist[[1]] <- levelplot(lev_p,
          col.regions = cols,
          colorkey = list(col = cols, 
                          at = do.breaks(range(lev_p), 
                                         color.ramp.length)),
          xlab = "", ylab = "",              # remove axis titles
          scales = list(x = list(rot = 45),  # change rotation for x-axis text
                        cex = 0.8), main='Bacteria' )     

glist[[2]] <- levelplot(lev_f,
                        col.regions = cols,
                        colorkey = list(col = cols, 
                                        at = do.breaks(range(lev_f), 
                                                       color.ramp.length)),
                        xlab = "", ylab = "",              # remove axis titles
                        scales = list(x = list(rot = 45),  # change rotation for x-axis text
                                      cex = 0.8), main='Fungi' )  

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 8, h= 3, filename=sprintf('%s/Levelplot_envPCA.pdf', dir$figdir))


############################################################################
# environmental gradient vs. PCoA axes: without latitude
# Prokaryotes

glist <- list()

glist[[1]] <- ggplot(prk, aes(x=environment_PC1, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(prk, aes(x=environment_PC2, y=PcoA.1)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x=environment_PC1, y=PcoA.2)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(prk, aes(x=environment_PC2, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
           filename=sprintf('%s/EnvGradientPCA_Family_prk.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes: without latitude
# Fungi

glist <- list()

glist[[1]] <- ggplot(fng, aes(x=environment_PC1, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x=environment_PC2, y=PcoA.1)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(fng, aes(x=environment_PC1, y=PcoA.2)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x=environment_PC2, y=PcoA.2)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
       filename=sprintf('%s/EnvGradientPCA_Family_fng.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes
# Prokaryotes

glist <- list()

glist[[1]] <- ggplot(prk, aes(x =environment_PC1, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[2]] <- ggplot(prk, aes(x =environment_PC2, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x =environment_PC1, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[4]] <- ggplot(prk, aes(x =environment_PC2, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
           filename=sprintf('%s/Contour_EnvGradientPCA_Family_prk.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes
# Fungi

glist <- list()

glist[[1]] <- ggplot(fng, aes(x =environment_PC1, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x =environment_PC2, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[3]] <- ggplot(fng, aes(x =environment_PC1, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x =environment_PC2, y = PcoA.2)) + geom_density_2d_filled()  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
       filename=sprintf('%s/Contour_EnvGradientPCA_Family_fng.pdf', dir$figdir))

############################################################################
# environmental gradient vs. PCoA axes: without latitude

glist <- list()

glist[[1]] <- ggplot(prk, aes(x=environment_PC1, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x=environment_PC1, y=PcoA.1)) + geom_point(color="dodgerblue3", size=0.4)  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x=environment_PC2, y=PcoA.1)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x=environment_PC2, y=PcoA.1)) + geom_point(color="deeppink3", size=0.4)  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
       filename=sprintf('%s/EnvGradientPCA_Family_summary1.pdf', dir$figdir))

############################################################################

glist <- list()

glist[[1]] <- ggplot(prk, aes(x =environment_PC1, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[2]] <- ggplot(fng, aes(x =environment_PC1, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[3]] <- ggplot(prk, aes(x =environment_PC2, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

glist[[4]] <- ggplot(fng, aes(x =environment_PC2, y = PcoA.1)) + geom_density_2d_filled()  +theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=glist, ncol=2), w= 7, h= 7,
       filename=sprintf('%s/EnvGradientPCA_Family_summary2.pdf', dir$figdir))

