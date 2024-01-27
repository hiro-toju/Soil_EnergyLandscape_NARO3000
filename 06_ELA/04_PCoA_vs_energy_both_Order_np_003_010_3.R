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
setwd("../")

par(oma = c(0.5, 0.5, 0.5, 0.5)) 
par(mar = c(0.5, 0.5, 0.5, 0.5)) 

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(mgcv)
library(rgl)
library(tidyverse)
library(plot3D)

# -- Create directory to save

rel <- 0.03
occ <- 0.10
prun <- 0

dir.create('06_ELA/04_output_both_Order')
dir <- make.dir(sprintf('06_ELA/04_output_both_Order/rel%s_occ%s', rel, occ))

############################################################################
# -- Load data table

#mat <- readRDS("Table/ELA_matrix_each_both_Order.rds")
m.pro <- readRDS('Table/ELA_matrix_pro_Order.rds')
m.fun <- readRDS('Table/ELA_matrix_fun_Order.rds')
colnames(m.pro) <- gsub(" ", "_", colnames(m.pro))
colnames(m.fun) <- gsub(" ", "_", colnames(m.fun))
#cbind(rownames(m.pro), rownames(m.fun))

m.pro  <- m.pro/rowSums(m.pro)
m.fun  <- m.fun/rowSums(m.fun)
mat <- cbind(m.pro, m.fun)

dist <- vegdist(mat, method='bray')
pca <- cmdscale(dist, eig=TRUE)$points
colnames(pca) <- c("PCoA.1", "PCoA.2")
pca_both <- data.frame(SampleID=rownames(pca), pca)

saveRDS(pca, 'Table/PCoA_Order_both_rrarefy.rds')

env <- readRDS('Table/env_1474.rds')
env <- data.frame(SampleID=rownames(env), env)

eng <- read.csv(sprintf("06_ELA/02_output_both_Order/rel%s_occ%s/Table/gStability_meanEnv_both_Order_1664_%s_%s.csv", rel, occ, rel, occ), header=T, row.names=1)
eng <- data.frame(SampleID=rownames(eng), eng)

#table(eng$stable.state.id.np)

pcoa.env <- left_join(pca_both, env, by="SampleID")
df <- left_join(pcoa.env, eng, by="SampleID")

write.csv(df, sprintf('%s/Engergy_both_Order_%s_%s.csv', dir$tabledir, rel, occ))

############################################################################
# Color

col1 <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Paired"), palettes(), brewer.pal(n = 8, name = "Dark2"))

############################################################################
# PCoA vs. energy indices

g <- list()

g[[1]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= e.realize))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle='Energy') + theme(legend.position = "none")

g[[2]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= ss.entropy.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Stable-state entropy') + theme(legend.position = "none")

g[[3]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= energy.gap.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Energy gap') + theme(legend.position = "none")

g[[4]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= energy.barrier.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Energy barrier') + theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=g, ncol=2), w= 9, h= 9,
           filename=sprintf('%s/EnergyIndices_PCoA_Order_%s_%s.pdf', dir$figdir, rel, occ))  


# caption

g <- list()

g[[1]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= e.realize))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Energy') + theme(legend.position = "right")

g[[2]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= ss.entropy.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Stable-state entropy') + theme(legend.position = "right")

g[[3]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= energy.gap.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Energy gap') + theme(legend.position = "right")

g[[4]] <- ggplot(df, aes(x=PCoA.1, y=PCoA.2, fill= energy.barrier.np))+
  geom_point(shape=21, alpha=0.8)+
  scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
  theme_bw(base_size=15)+
  labs(x='PCoA 1', y='PCoA 2', subtitle='Energy barrier') + theme(legend.position = "right")

ggsave(plot=plot_grid( plotlist=g, ncol=2), w= 11.5, h= 9,
       filename=sprintf('%s/EnergyIndices_PCoA_Order_%s_%s_caption.pdf', dir$figdir, rel, occ))  
         
############################################################################
# levelplot

sml <- readRDS("Table/sml_1474.rds") 
sml <- data.frame(SampleID=rownames(sml), sml, DL.bin=sml$DL)

sml$DL.bin[which(sml$DL.bin=='Level 1')] <- 1
sml$DL.bin[which(sml$DL.bin=='Level 2')] <- 0
sml$DL.bin[which(sml$DL.bin=='Level 3')] <- 0
sml$DL.bin[which(sml$DL.bin=='Level 4')] <- 0
sml$DL.bin[which(sml$DL.bin=='Level 5')] <- 0
mode(sml$DL.bin) <- "numeric"

df2 <- left_join(df, sml[, c('SampleID', 'DL.bin')], by="SampleID")
df2_2 <- na.omit(df2[, c('SampleID', 'stable.state.id.np', 'DL.bin')])


# Excluding minor stable states

df3 <- matrix(NA, nrow=0, ncol=3)

for(i in 1:length(unique(df2_2$stable.state.id.np))) {
	a <- subset(df2_2, df2_2$stable.state.id.np ==unique(df2_2$stable.state.id.np)[i])
	b <- nrow(a)
		if (b >= 10 ) {
		df3 <- rbind(df3, a)
	} 
}

df4 <- df3[order(df3$SampleID), ]

comm <- mat[rownames(mat) %in%  df3$SampleID , ]
comm2 <- comm[order(rownames(comm)),]

#cbind(rownames(comm2), df4$SampleID)

comm3 <- matrix(NA, nrow=length(unique(df4$stable.state.id.np)), ncol=ncol(comm2))

for (i in 1:ncol(comm2)) {
	 comm3[,i] <- tapply(as.numeric(comm2[,i]), df4$stable.state.id.np, mean)
}

stdl <- tapply(df4$DL.bin, df4$stable.state.id.np, mean)
stateDL <- data.frame(cbind(ID =names(stdl), MeanDL=as.numeric(stdl)))

stateDLreorder <- stateDL[order(stateDL$MeanDL, decreasing=T),]
write.table(stateDLreorder, sprintf('%s/StableStates_MealDL_Order_%s_%s.csv', dir$tabledir, rel, occ))

#SS.ID <- paste("SS_",formatC(1:nrow(stateDLreorder), width=2, flag="0"), sep='')
#stateDLreorder2 <- cbind(SS.ID, stateDLreorder )

############################################################################
# Stable states DL
#colpl <- c("blue",  "mistyrose",  "red")
#colpl <- c("blue", rev(brewer.pal(3, "Blues")), "lightpink", "red")
#colpl <- c("blue", rev(brewer.pal(6, "Blues")), 'mistyrose', "lightpink", "deeppink2", "red", "brown4")
colpl <- c("blue4", rev(brewer.pal(5, "Blues")), brewer.pal(4, "Reds"),  "brown4")
#colpl <- c("blue", "skyblue1", "lightpink", "deeppink1")
#colpl <- c("blue", "skyblue1",  "deeppink1")

g <- ggplot(stateDLreorder, aes(x= reorder(ID,  -as.numeric(MeanDL)), y= as.numeric(MeanDL), fill=reorder(ID,  -as.numeric(MeanDL)))) + geom_bar(stat = "identity") + theme(axis.text=element_text(angle=45)) + xlab('Stable state ID') + ylab('prop. disease severity level < 20') + scale_fill_manual(values= colpl)

ggsave(g, w= 7, h= 3, filename=sprintf('%s/Barplot_StableStates_DL_Order_%s_%s.pdf', dir$figdir, rel, occ))  
############################################################################
# Stable states 

rownames(comm3) <- stateDL$ID
colnames(comm3) <- colnames(comm2)

comm3reorder <- comm3[order(stateDL$MeanDL, decreasing=T), ]

comm4 <- t(comm3reorder)
comm5 <- subset(comm4, rowSums(comm4) > 0)
comm6 <- subset(comm5, rownames(comm5)!="Unidentified")
#comm6 <- cbind(Order=rownames(comm5),comm5)

comm6 <- comm6[c(sort(rownames(comm6[229:313,]), decreasing=T), sort(rownames(comm6[1:228,]), decreasing=T)), ]


############################################################################
# levelplot all taxa

pdf(file=sprintf('%s/Levelplot_StableStates_Order_continuous_%s_%s.pdf', dir$figdir, rel, occ), width=6, height=8)
levelplot(as.matrix(t(comm6[, stateDLreorder$ID])), xlab="Stable state ID", ylab="Order",
          scales = list(x = list(rot = 45), cex = 0.8), col.regions=colorRampPalette
(c("white", "deeppink3 ")))      
dev.off()

pdf(file=sprintf('%s/Levelplot_StableStates_Order_continuous_incUnidentified_%s_%s.pdf', dir$figdir, rel, occ), width=6, height=8)
levelplot(as.matrix(t(comm5[, stateDLreorder$ID])), xlab="Stable state ID", ylab="Order",
          scales = list(x = list(rot = 45), cex = 0.8), col.regions=colorRampPalette
          (c("white", "deeppink3 ")))      
dev.off()


# levelplot selected taxa: all stable states

ss <- readRDS(sprintf("Table/gStability_meanEnv_both_Order_1664_%s_%s.rds", rel, occ))[[3]][[2]]
ss2 <- ss[,3:ncol(ss)]
rownames(ss2) <- ss$ID

ss3 <- t(as.matrix(ss2))
mode(ss3) <- "numeric"

ss4 <- subset(ss3, rowSums(ss3) > 0)
ss5 <- subset(ss4, rownames(ss4)!="Unidentified")

pdf(file=sprintf('%s/Levelplot_StableStates_Order_binary_%s_%s.pdf', dir$figdir, rel, occ), width=6, height=8)
levelplot(t(ss5[order(rownames(ss5), decreasing=T),]), xlab="Stable state ID", ylab="Order",
          scales = list(x = list(rot = 45), cex = 0.8), col.regions=colorRampPalette
          (c("white", "deeppink3 ")))      
dev.off()

# levelplot selected taxa: major stable states

ss6 <- ss5[, colnames(ss5) %in% stateDLreorder$ID]
ss7 <- subset(ss6, rowSums(ss6) > 0)
ss8 <- subset(ss7, rowSums(ss7) < ncol(ss7))
ss9 <- t(ss8)

ss9 <- ss9[, c(sort(colnames(ss9[ ,16:20]), decreasing=T), sort(colnames(ss9[ ,1:15]), decreasing=T))]

pdf(file=sprintf('%s/Levelplot_StableStates_Order_binary_majorSS_%s_%s.pdf', dir$figdir, rel, occ), width=6, height=8)
levelplot(ss9[stateDLreorder$ID, -1], xlab="Stable state ID", ylab="Order",
          scales = list(x = list(rot = 45), cex = 0.8), col.regions=colorRampPalette
          (c("white", "deeppink3 ")))      # "Unidentified.1" excluded
dev.off()


#######################################################
# Landscape plotting


pdf(file=sprintf('%s/Landscapes_3D_%s_%s_1.pdf', dir$figdir, rel, occ), width=12, height=6)

par(mfrow=c(1,2))
#par(oma = c(1, 1, 1, 1)) 
#par(mar = c(1, 1, 1, 1)) 

mod2 <- gam(e.realize ~ s(PCoA.1, PCoA.2), data= df2)
vis.gam(mod2, ticktype='detailed', color="bw", theta=40, phi=60, n.grid=50, zlab="energy")
#persp(PCoA_1, PCoA_2, Energy, theta = 20, phi=60, ticktype="detailed", zlab="energy", col = "khaki1")

mod1 <- gam(DL.bin ~ s(PCoA.1, PCoA.2), data= df2)
vis.gam(mod1, ticktype='detailed', color="heat", theta=40, phi=60, n.grid=50, zlab="prop. disease severity level < 20")

dev.off()

#######################################################

pdf(file=sprintf('%s/Landscapes_3D_Order_%s_%s_2.pdf', dir$figdir, rel, occ), width=12, height=12)

par(mfrow=c(2,2))
#par(oma = c(1, 1, 1, 1)) 
#par(mar = c(1, 1, 1, 1)) 

mod2 <- gam(e.realize ~ s(PCoA.1, PCoA.2), data= df2)

steps <- 40
PCoA_1 <- with(df2, seq(min(PCoA.1), max(PCoA.1), length = steps) )
PCoA_2 <- with(df2, seq(min(PCoA.2), max(PCoA.2), length = steps) )
newdat <- expand.grid(PCoA.1 = PCoA_1, PCoA.2 = PCoA_2)
Energy <- matrix(predict(mod2, newdat), steps, steps)

colors <- rep("grey80", times=nrow(df2))
states <- unique(stateDLreorder$ID)

for (i in 1:length(states)) {
	colors[df2$stable.state.id.np == states[i]] <- colpl[i]
}

p <- persp(PCoA_1, PCoA_2, Energy, theta = 40, phi=50, col = "khaki1", ticktype="detailed", zlim=c(min(df2$e.realize), max(df2$e.realize)), zlab="energy")
obs <- with(df2, trans3d(PCoA.1, PCoA.2, e.realize, p))
pred <- with(df2, trans3d(PCoA.1, PCoA.2, fitted(mod2), p))
points(obs, col = colors, pch = 16, cex=0.8)

p <- persp(PCoA_1, PCoA_2, Energy, theta = 40, phi=40, col = "khaki1", ticktype="detailed", zlim=c(min(df2$e.realize), max(df2$e.realize)), zlab="energy")
obs <- with(df2, trans3d(PCoA.1, PCoA.2, e.realize, p))
pred <- with(df2, trans3d(PCoA.1, PCoA.2, fitted(mod2), p))
points(obs, col = colors, pch = 16, cex=0.8)

p <- persp(PCoA_1, PCoA_2, Energy, theta = 40, phi=30, col = "khaki1", ticktype="detailed", zlim=c(min(df2$e.realize), max(df2$e.realize)), zlab="energy")
obs <- with(df2, trans3d(PCoA.1, PCoA.2, e.realize, p))
pred <- with(df2, trans3d(PCoA.1, PCoA.2, fitted(mod2), p))
points(obs, col = colors, pch = 16, cex=0.8)

p <- persp(PCoA_1, PCoA_2, Energy, theta = 40, phi=20, col = "khaki1", ticktype="detailed", zlim=c(min(df2$e.realize), max(df2$e.realize)), zlab="energy")
obs <- with(df2, trans3d(PCoA.1, PCoA.2, e.realize, p))
pred <- with(df2, trans3d(PCoA.1, PCoA.2, fitted(mod2), p))
points(obs, col = colors, pch = 16, cex=0.8)

dev.off()
#######################################################

