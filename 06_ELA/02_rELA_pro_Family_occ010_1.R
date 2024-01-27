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

library('AnalysisHelper')

## -- Create directory to save
dir.create('06_ELA/02_output_pro_Family')
#dir <- make.dir('06_ELA/02_output')

library(renv)
library(parallel)
library(foreach)
library(renv)
library("Rcpp")
library("RcppArmadillo")
library("doParallel")
library('tidyverse')
library('gsubfn')
library('zoo')
library('snow')
library('plyr')
library('gtools')
library('ggsci')
library('igraph')
library('tidygraph')
library('RColorBrewer')
library("stringdist")
#library('scatterpie')

library("rELA")

detectCores()
n.core <- 32

# loading data
baseabtable <- readRDS('Table/ELA_matrix_each_pro_Family.rds')
env <- readRDS('Table/env_pro_1664.rds')
basemetadata <- data.frame(apply(env, 2, scale))
rownames(basemetadata) <- rownames(env)

#########################################################################
#########################################################################
## Analysis with environmental data (metadata)
## Setting 1

rel <- 0.01
occ <- 0.10
prun <- 0
qth <- 10^-5

dir <- make.dir(sprintf('06_ELA/02_output_pro_Family/rel%s_occ%s', rel, occ))

list[ocmat, abmat, enmat, samplelabel, specieslabel, factorlabel] <- Formatting(baseabtable, basemetadata, 1, c(rel, occ, 1 - occ))

saveRDS(ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.rds", dir$rdsdir, rel, occ))
write.csv(x = ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) 

dim(ocmat)
sum(ocmat)

sa <- runSA(ocmat=as.matrix(ocmat), enmat=as.matrix(enmat), qth=qth, rep=256, threads=n.core)
saveRDS(sa, file=sprintf("Table/runSA_pro_Family_1664_%s_%s.rds", rel, occ))

#########################################################################
# Disconnectivity graph

gela <- GradELA(sa=sa, eid="pH_dry_soil",  enmat=enmat, env=NULL, range=NULL, steps=32, th=prun, threads=n.core)# Specify the label or position of an environmental factor #[[1]]: return value of ELA function for each step, [[2]]: value of environmental factor for each step, [[3]]: specified environmental factor

pdf(file=sprintf("%s/showDG1_gela_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(gela[[1]][[1]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showDG2_gela_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(gela[[1]][[2]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showSSD_gela_pro_Family_1664_%s_%s_%s.pdf", dir$figdir, rel, occ, prun))
showSSD(gela)
dev.off()

#########################################################################
## GStability Manual

#gStability returns a list of 4 elements: the first two are the dataframe for pruned/non-pruned energy landscape, respectively. In addition to the dataframe of Stability it includes e.tipping (energy of tipping point) and energy.barrier (height of energy from observed state to the tipping point).The third output is a list of parameters (h, g, j, h+g*env) and a summary table of stable states, and the fourth output is a list encapsulating the inputs required for the various plots.

#output of gStability:
#[[1]]: data.frame(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping, state.id, stable.state.id)
#[[2]]: data.frame(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np, state.id.np, stable.state.id.np)
#[[3]]: w/ enmat: list(list(list(he, je, ge, hge), data.frame(sstable)), ...); w/o enmat: list(list(he, je, ge, hge), data.frame(sstable))
#[[4]]: w/ enmat: list(list(ocmat, env, sa, ela, elanp), ...); w/o enmat: list(ocmat, env, sa, ela, elanp)

#########################################################################
# gStability

# assuming mean values of environmental variables
gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core) # assuming mean values of environmental variables
saveRDS(gstb, file=sprintf("Table/gStability_meanEnv_pro_Family_1664_%s_%s.rds", rel, occ))
write.csv(x = gstb[[1]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s_prun%s.csv", dir$tabledir, rel, occ, prun)) #pruned
write.csv(x = gstb[[2]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) #non-pruned

#########################################################################
# PCA 

# if enmat=NULL, remove "[[sample.id]]"
ocmat <- gstb[[4]][[1]]
env <- gstb[[4]][[2]]
sa <- gstb[[4]][[3]]
ela <- gstb[[4]][[4]]

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
PCplot(ocmat, sa, ssrep=ela[[2]], pruned=FALSE)
dev.off()

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
PCplot(ocmat, sa, ssrep=ela[[2]])
dev.off()

#########################################################################

ocmat <- gstb[[4]][[1]]
ela <- gstb[[4]][[4]]
elanp <- gstb[[4]][[5]]

pdf(file=sprintf("%s/showDG_gstab4_env_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(ela[[1]], ocmat, "pruned")
dev.off()

pdf(file=sprintf("%s/showDG_gstab5_env_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(elanp, ocmat, "not-pruned")
dev.off()

#########################################################################
#########################################################################


#########################################################################
#########################################################################
## Analysis with environmental data (metadata)
## Setting 2

rel <- 0.005
occ <- 0.10
prun <- 0
qth <- 10^-5

dir <- make.dir(sprintf('06_ELA/02_output_pro_Family/rel%s_occ%s', rel, occ))

list[ocmat, abmat, enmat, samplelabel, specieslabel, factorlabel] <- Formatting(baseabtable, basemetadata, 1, c(rel, occ, 1 - occ))

saveRDS(ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.rds", dir$rdsdir, rel, occ))
write.csv(x = ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) 

dim(ocmat)
sum(ocmat)

sa <- runSA(ocmat=as.matrix(ocmat), enmat=as.matrix(enmat), qth=qth, rep=256, threads=n.core)
saveRDS(sa, file=sprintf("Table/runSA_pro_Family_1664_%s_%s.rds", rel, occ))

#########################################################################
# Disconnectivity graph

gela <- GradELA(sa=sa, eid="pH_dry_soil",  enmat=enmat, env=NULL, range=NULL, steps=32, th=prun, threads=n.core)# Specify the label or position of an environmental factor #[[1]]: return value of ELA function for each step, [[2]]: value of environmental factor for each step, [[3]]: specified environmental factor

pdf(file=sprintf("%s/showDG1_gela_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(gela[[1]][[1]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showDG2_gela_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(gela[[1]][[2]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showSSD_gela_pro_Family_1664_%s_%s_%s.pdf", dir$figdir, rel, occ, prun))
showSSD(gela)
dev.off()

#########################################################################
## GStability Manual

#gStability returns a list of 4 elements: the first two are the dataframe for pruned/non-pruned energy landscape, respectively. In addition to the dataframe of Stability it includes e.tipping (energy of tipping point) and energy.barrier (height of energy from observed state to the tipping point).The third output is a list of parameters (h, g, j, h+g*env) and a summary table of stable states, and the fourth output is a list encapsulating the inputs required for the various plots.

#output of gStability:
#[[1]]: data.frame(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping, state.id, stable.state.id)
#[[2]]: data.frame(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np, state.id.np, stable.state.id.np)
#[[3]]: w/ enmat: list(list(list(he, je, ge, hge), data.frame(sstable)), ...); w/o enmat: list(list(he, je, ge, hge), data.frame(sstable))
#[[4]]: w/ enmat: list(list(ocmat, env, sa, ela, elanp), ...); w/o enmat: list(ocmat, env, sa, ela, elanp)

#########################################################################
# gStability

# assuming mean values of environmental variables
gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core) # assuming mean values of environmental variables
saveRDS(gstb, file=sprintf("Table/gStability_meanEnv_pro_Family_1664_%s_%s.rds", rel, occ))
write.csv(x = gstb[[1]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s_prun%s.csv", dir$tabledir, rel, occ, prun)) #pruned
write.csv(x = gstb[[2]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) #non-pruned

#########################################################################
# PCA 

# if enmat=NULL, remove "[[sample.id]]"
ocmat <- gstb[[4]][[1]]
env <- gstb[[4]][[2]]
sa <- gstb[[4]][[3]]
ela <- gstb[[4]][[4]]

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
PCplot(ocmat, sa, ssrep=ela[[2]], pruned=FALSE)
dev.off()

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
PCplot(ocmat, sa, ssrep=ela[[2]])
dev.off()

#########################################################################

ocmat <- gstb[[4]][[1]]
ela <- gstb[[4]][[4]]
elanp <- gstb[[4]][[5]]

pdf(file=sprintf("%s/showDG_gstab4_env_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(ela[[1]], ocmat, "pruned")
dev.off()

pdf(file=sprintf("%s/showDG_gstab5_env_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(elanp, ocmat, "not-pruned")
dev.off()

#########################################################################
#########################################################################



#########################################################################
#########################################################################
## Analysis with environmental data (metadata)
## Setting 3

rel <- 0.001
occ <- 0.10
prun <- 0
qth <- 10^-5

dir <- make.dir(sprintf('06_ELA/02_output_pro_Family/rel%s_occ%s', rel, occ))

list[ocmat, abmat, enmat, samplelabel, specieslabel, factorlabel] <- Formatting(baseabtable, basemetadata, 1, c(rel, occ, 1 - occ))

saveRDS(ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.rds", dir$rdsdir, rel, occ))
write.csv(x = ocmat, file=sprintf("%s/ocmat_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) 

dim(ocmat)
sum(ocmat)

sa <- runSA(ocmat=as.matrix(ocmat), enmat=as.matrix(enmat), qth=qth, rep=256, threads=n.core)
saveRDS(sa, file=sprintf("Table/runSA_pro_Family_1664_%s_%s.rds", rel, occ))

#########################################################################
# Disconnectivity graph

gela <- GradELA(sa=sa, eid="pH_dry_soil",  enmat=enmat, env=NULL, range=NULL, steps=32, th=prun, threads=n.core)# Specify the label or position of an environmental factor #[[1]]: return value of ELA function for each step, [[2]]: value of environmental factor for each step, [[3]]: specified environmental factor

pdf(file=sprintf("%s/showDG1_gela_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(gela[[1]][[1]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showDG2_gela_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(gela[[1]][[2]][[1]], ocmat)
dev.off()

pdf(file=sprintf("%s/showSSD_gela_pro_Family_1664_%s_%s_%s.pdf", dir$figdir, rel, occ, prun))
showSSD(gela)
dev.off()

#########################################################################
## GStability Manual

#gStability returns a list of 4 elements: the first two are the dataframe for pruned/non-pruned energy landscape, respectively. In addition to the dataframe of Stability it includes e.tipping (energy of tipping point) and energy.barrier (height of energy from observed state to the tipping point).The third output is a list of parameters (h, g, j, h+g*env) and a summary table of stable states, and the fourth output is a list encapsulating the inputs required for the various plots.

#output of gStability:
#[[1]]: data.frame(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping, state.id, stable.state.id)
#[[2]]: data.frame(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np, state.id.np, stable.state.id.np)
#[[3]]: w/ enmat: list(list(list(he, je, ge, hge), data.frame(sstable)), ...); w/o enmat: list(list(he, je, ge, hge), data.frame(sstable))
#[[4]]: w/ enmat: list(list(ocmat, env, sa, ela, elanp), ...); w/o enmat: list(ocmat, env, sa, ela, elanp)

#########################################################################
# gStability

# assuming mean values of environmental variables
gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core) # assuming mean values of environmental variables
saveRDS(gstb, file=sprintf("Table/gStability_meanEnv_pro_Family_1664_%s_%s.rds", rel, occ))
write.csv(x = gstb[[1]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s_prun%s.csv", dir$tabledir, rel, occ, prun)) #pruned
write.csv(x = gstb[[2]], file=sprintf("%s/gStability_meanEnv_pro_Family_1664_%s_%s.csv", dir$tabledir, rel, occ)) #non-pruned

#########################################################################
# PCA 

# if enmat=NULL, remove "[[sample.id]]"
ocmat <- gstb[[4]][[1]]
env <- gstb[[4]][[2]]
sa <- gstb[[4]][[3]]
ela <- gstb[[4]][[4]]

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
PCplot(ocmat, sa, ssrep=ela[[2]], pruned=FALSE)
dev.off()

pdf(file=sprintf("%s/PCplot_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
PCplot(ocmat, sa, ssrep=ela[[2]])
dev.off()

#########################################################################

ocmat <- gstb[[4]][[1]]
ela <- gstb[[4]][[4]]
elanp <- gstb[[4]][[5]]

pdf(file=sprintf("%s/showDG_gstab4_env_pro_Family_1664_%s_%s_prun%s.pdf", dir$figdir, rel, occ, prun))
showDG(ela[[1]], ocmat, "pruned")
dev.off()

pdf(file=sprintf("%s/showDG_gstab5_env_pro_Family_1664_%s_%s.pdf", dir$figdir, rel, occ))
showDG(elanp, ocmat, "not-pruned")
dev.off()

#########################################################################
#########################################################################



