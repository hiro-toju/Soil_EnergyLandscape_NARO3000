############################################################################
####
#### R script for Fujita (2022)
####
#### Energy landscape analysis
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")


## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(mgcv)
library(rgl)

# -- Create directory to save
dir <- make.dir('06_ELA/03_output')


g <- list()

############################################################################

th <- "th5"
lv <- "Family"

data <- list()

data[[1]] <- readRDS(sprintf('Table/prok_fungi_%s/ELAstability_%s.rds', th, lv))
data[[2]] <- readRDS(sprintf('Table/prok_%s/ELAstability_%s.rds', th, lv))
data[[3]] <- readRDS(sprintf('Table/fungi_%s/ELAstability_%s.rds', th, lv))

############################################################################


mat <- rbind(data[[1]], data[[2]], data[[3]])

Dataset <- c(rep("1_Prokaryotes_Fungi", times=nrow(data[[1]])), rep("2_Prokaryotes", times=nrow(data[[2]])), rep("3_Fungi", times=nrow(data[[3]])))

df1 <- cbind(Dataset, mat)

g[[1]] <- ggplot(df1) + geom_boxplot(aes(x=Dataset, y=log10(as.numeric(nSS)), fill= Dataset)) + labs(x='', y='log10 (number of stable states)') + ggtitle(sprintf("%s level (%s)", lv, th)) + theme(legend.position="none") + scale_y_continuous(limits=c(0.2, 4.3))

	
	
############################################################################



############################################################################

th <- "th10"
lv <- "Family"

data <- list()

data[[1]] <- readRDS(sprintf('Table/prok_fungi_%s/ELAstability_%s.rds', th, lv))
data[[2]] <- readRDS(sprintf('Table/prok_%s/ELAstability_%s.rds', th, lv))
data[[3]] <- readRDS(sprintf('Table/fungi_%s/ELAstability_%s.rds', th, lv))[,c(1,2,3,5,4)]

############################################################################

mat <- rbind(data[[1]], data[[2]], data[[3]])

Dataset <- c(rep("1_Prokaryotes_Fungi", times=nrow(data[[1]])), rep("2_Prokaryotes", times=nrow(data[[2]])), rep("3_Fungi", times=nrow(data[[3]])))

df1 <- cbind(Dataset, mat)

g[[2]] <- ggplot(df1) + geom_boxplot(aes(x=Dataset, y=log10(as.numeric(nSS)), fill= Dataset)) + labs(x='', y='log10 (number of stable states)') + ggtitle(sprintf("%s level (%s)", lv, th)) + theme(legend.position="none") + scale_y_continuous(limits=c(0.2, 4.3))

	
############################################################################



############################################################################

th <- "th5"
lv <- "Order"

data <- list()

data[[1]] <- readRDS(sprintf('Table/prok_fungi_%s/ELAstability_%s.rds', th, lv))
data[[2]] <- readRDS(sprintf('Table/prok_%s/ELAstability_%s.rds', th, lv))
data[[3]] <- readRDS(sprintf('Table/fungi_%s/ELAstability_%s.rds', th, lv))

############################################################################


mat <- rbind(data[[1]], data[[2]], data[[3]])

Dataset <- c(rep("1_Prokaryotes_Fungi", times=nrow(data[[1]])), rep("2_Prokaryotes", times=nrow(data[[2]])), rep("3_Fungi", times=nrow(data[[3]])))

df1 <- cbind(Dataset, mat)

g[[3]] <- ggplot(df1) + geom_boxplot(aes(x=Dataset, y=log10(as.numeric(nSS)), fill= Dataset)) + labs(x='', y='log10 (number of stable states)') + ggtitle(sprintf("%s level (%s)", lv, th)) + theme(legend.position="none") + scale_y_continuous(limits=c(0.2, 4.3))



############################################################################



############################################################################

th <- "th10"
lv <- "Order"

data <- list()

data[[1]] <- readRDS(sprintf('Table/prok_fungi_%s/ELAstability_%s.rds', th, lv))
data[[2]] <- readRDS(sprintf('Table/prok_%s/ELAstability_%s.rds', th, lv))
data[[3]] <- readRDS(sprintf('Table/fungi_%s/ELAstability_%s.rds', th, lv))[,c(1,2,3,5,4)]

############################################################################

mat <- rbind(data[[1]], data[[2]], data[[3]])

Dataset <- c(rep("1_Prokaryotes_Fungi", times=nrow(data[[1]])), rep("2_Prokaryotes", times=nrow(data[[2]])), rep("3_Fungi", times=nrow(data[[3]])))

df1 <- cbind(Dataset, mat)

g[[4]] <- ggplot(df1) + geom_boxplot(aes(x=Dataset, y=log10(as.numeric(nSS)), fill= Dataset)) + labs(x='', y='log10 (number of stable states)') + ggtitle(sprintf("%s level (%s)", lv, th)) + theme(legend.position="none") + scale_y_continuous(limits=c(0.2, 4.3))


	
############################################################################


ggsave(plot=plot_grid(plotlist=g, ncol=4), w= 8, h= 2.3,
           filename=sprintf('%s/ELA_StableStates_Number_Comparison.pdf', dir$figdir)) 
           
            