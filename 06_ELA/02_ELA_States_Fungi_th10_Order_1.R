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

target <- "fungi"
th <- "th10"
lv <- "Order"

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(mgcv)
library(rgl)

# -- Create directory to save
dir <- make.dir('06_ELA/02_output')


############################################################################


# -- Load data table
#dflist <- readRDS('Table/0303_rarefiedMatrix.rds') # should be replaced with 'Table/0303_rarefiedMatrix_rrarefy.rds'

#sml <- readRDS("Table/0301_sml.rds") 

#comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')

#commonrow <- intersect(rownames(comm16s), rownames(commits))

#commonmat <- cbind(comm16s[commonrow, colSums(comm16s)>0],
#commits[commonrow, colSums(commits)>0])

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa <- rbind(taxa16S, taxaITS)

score <- readRDS('05_Community_structure/05_output/RDS/PCoA.Table.Fng.rds')



############################################################################

taxalevel <- 'Family'
n.core <- parallel::detectCores()

#taxamat <-  Taxa.mat(commonmat,taxa=taxa, taxaLabel=taxalevel)
taxamat <-  Taxa.mat(commits,taxa=taxa, taxaLabel=taxalevel)



############################################################################

stb <- readRDS(sprintf('Table/%s_%s/ELAstability_%s.rds', target, th, lv))
SScomm <- readRDS(sprintf('Table/%s_%s/SScomm_%s.rds', target, th, lv))

############################################################################


stbt <- cbind(Sample.ID=rownames(stb), stb)

write.table(stbt, file=sprintf('%s/StableStates_Energy_%s_%s_%s.txt', dir$tabledir, target, th, lv), sep='\t', quote=F, row.names=F)

SScommt <- cbind(Sample.ID=rownames(SScomm), SScomm[, -ncol(SScomm)])

write.table(SScommt, file=sprintf('%s/StableStates_Communities_%s_%s_%s.txt', dir$tabledir, target, th, lv), sep='\t', quote=F, row.names=F)

SS <- merge(stbt, SScommt, by="Sample.ID")

SS3 <- merge(SS, score, by="Sample.ID")


write.table(SS3, file=sprintf('%s/StableStates_%s_%s_%s.txt', dir$tabledir, target, th, lv), sep='\t', quote=F, row.names=F)

saveRDS(SS3, file=sprintf('%s/StableStates_%s_%s_%s.rds', dir$rdsdir, target, th, lv))

############################################################################
# PCoA stable staets 1

col1 <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Paired"), palettes(), brewer.pal(n = 8, name = "Dark2"))

#mod1 <- gam(as.numeric(energy) ~ s(as.numeric(PCoA1), as.numeric(PCoA2)), data= SS3)

#vis.gam(mod1, color="cm", plot.type="contour")

          

############################################################################
# PCoA stable staets 2

g <- list()

g[[1]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= SSenergy))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle='Stable-state energy') + theme(legend.position = "none")

g[[2]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= energy))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle="Energy") + theme(legend.position = "none")     
        
g[[3]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= energy.gap))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle='Energy gap') + theme(legend.position = "none")     
         
        
ggsave(plot=plot_grid( plotlist=g, ncol=2), w= 9, h= 9,
           filename=sprintf('%s/StableStates_PCoA_%s_%s_%s.pdf', dir$figdir, target, th, lv))  
        

g <- list()

g[[1]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= SSenergy))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle='Stable-state energy') + theme(legend.position = "right")

g[[2]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= energy))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle="Energy") + theme(legend.position = "right")     
        
g[[3]] <- ggplot(SS3, aes(x=as.numeric(PCoA1), y=as.numeric(PCoA2), fill= energy.gap))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(9, 'YlOrRd'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(x='PCoA 1', y='PCoA 2', subtitle='Energy gap') + theme(legend.position = "right")     
        

ggsave(plot=plot_grid( plotlist=g, ncol=2), w= 11, h= 9,
           filename=sprintf('%s/StableStates_PCoA_%s_%s_%s_caption.pdf', dir$figdir, target, th, lv))  
         


############################################################################
# levelplot



df2 <- na.omit(SS3[, c('Sample.ID', 'StableState.id', 'DL.bin')])

df3 <- matrix(NA, nrow=0, ncol=3)

for(i in 1:length(unique(df2$StableState.id))) {
	a <- subset(df2, df2$StableState.id ==unique(df2$StableState.id)[i])
	b <- nrow(a)
		if (b >= 10 ) {
		df3 <- rbind(df3, a)
	} 
}


write.table(df3, file=sprintf('%s/StableStates_states_%s_%s_%s.txt', dir$tabledir, target, th, lv), sep='\t', quote=F, row.names=F)

df4 <- df3[order(df3$Sample.ID), ]

comm <- SS3[SS3$Sample.ID %in%  df3$Sample.ID , 8:(ncol(SS3)-35)]
rownames(comm) <- SS3[SS3$Sample.ID %in%  df3$Sample.ID , 1]
comm2 <- comm[order(rownames(comm)),]

write.table(cbind(rownames(comm), comm), file=sprintf('%s/StableStates_comm_%s_%s_%s.txt', dir$tabledir, target, th, lv), sep='\t', quote=F, row.names=F)

cbind(rownames(comm2), df4$Sample.ID)

comm3 <- matrix(NA, nrow=length(unique(df4$StableState.id)), ncol=ncol(comm2))

for (i in 1:ncol(comm2)) {
	 comm3[,i] <- tapply(as.numeric(comm2[,i]), df4$StableState.id, mean)
}

stdl <- tapply(df4$DL.bin, df4$StableState.id, mean)
stateDL <- data.frame(cbind(ID =names(stdl), MeanDL=as.numeric(stdl)))

stateDLreorder <- stateDL[order(stateDL$MeanDL, decreasing=T),]

SS.ID <- paste("SS_",formatC(1:nrow(stateDLreorder), width=2, flag="0"), sep='')
stateDLreorder2 <- cbind(SS.ID, stateDLreorder )

write.table(stateDLreorder2, file=sprintf('%s/Major_StableStates_%s_%s_%s.txt', dir$tabledir, target, th, lv), quote=F, sep='\t', row.names=F)

colpl <- c(rev(brewer.pal(5, "Blues")), "oldlace  ", "mistyrose")

g <- ggplot(stateDLreorder2, aes(x= reorder(SS.ID,  -as.numeric(MeanDL)), y= as.numeric(MeanDL), fill=SS.ID)) + geom_bar(stat = "identity") + theme(axis.text=element_text(angle=45)) + xlab('Stable state ID') + ylab('Proportion of disease level 1') + scale_fill_manual(values= colpl)

ggsave(g, w= 7, h= 3, filename=sprintf('%s/Barplot_StableStates_DL_%s_%s_%s.pdf', dir$figdir, target, th, lv))  

rownames(comm3) <- stateDL$ID
colnames(comm3) <- colnames(comm2)

comm3reorder <- comm3[order(stateDL$MeanDL, decreasing=T), ]

comm4 <- t(comm3reorder)
comm5 <- subset(comm4, rowSums(comm4) > 0)

comm6 <- cbind(Family=rownames(comm5),comm5)

famkin <- paste(taxa$Family, taxa$Kingdom, sep='_')
damkinunique <- unique(famkin)
vec <- unlist(strsplit(damkinunique, "_"))
famkinlist <- data.frame(Family=vec[seq(1,length(vec),2)], Kingdom=vec[seq(2,length(vec),2)])
colnames(famkinlist) <- c("Family", "Kingdom")

comm7 <- merge(comm6, famkinlist, by='Family', all=F, all.x=TRUE, all.y=FALSE)

comm7 
comm5

pdf(file=sprintf('%s/Levelplot_StableStates_%s_%s_%s.pdf', dir$figdir, target, th, lv), width=6, height=5)

levelplot(as.matrix(t(comm5[order(rownames(comm5), decreasing=T),])), xlab="Stable state ID", ylab="Family",
          scales = list(x = list(rot = 45), cex = 0.8), col.regions=colorRampPalette
(c("white", "deeppink3 ")))      
dev.off()



#######################################################
# Landscape plotting


pdf(file=sprintf('%s/Landscapes_3D_%s_%s_%s_1.pdf', dir$figdir, target, th, lv), width=4, height=8)

par(mfrow=c(2,1))
par(oma = c(1, 1, 1, 1)) 
par(mar = c(1, 1, 1, 1)) 

mod1 <- gam(DL.bin ~ s(PCoA1, PCoA2), data= SS3)
vis.gam(mod1, ticktype='detailed', color="cm", theta=20, phi=60, main="Disease level")

mod2 <- gam(energy ~ s(PCoA1, PCoA2), data= SS3)
vis.gam(mod2, ticktype='detailed', color="cm", theta=20, phi=60, main="Energy")

dev.off()


#######################################################

pdf(file=sprintf('%s/Landscapes_3D_%s_%s_%s_2.pdf', dir$figdir, target, th, lv), width=4, height=4)

par(mfrow=c(1,1))
par(oma = c(1, 1, 1, 1)) 
par(mar = c(1, 1, 1, 1)) 

mod2 <- gam(energy ~ s(PCoA1, PCoA2), data= SS3)

steps <- 40
PCoA_1 <- with(SS3, seq(min(PCoA1), max(PCoA1), length = steps) )
PCoA_2 <- with(SS3, seq(min(PCoA2), max(PCoA2), length = steps) )
newdat <- expand.grid(PCoA1 = PCoA_1, PCoA2 = PCoA_2)
Energy <- matrix(predict(mod2, newdat), steps, steps)

p <- persp(PCoA_1, PCoA_2, Energy, theta = 20, phi=60, col = "khaki1  ", ticktype="detailed")

colors <- rep("grey80", times=nrow(SS3))
states <- unique(stateDLreorder$ID)

for (i in 1:length(states)) {
	colors[SS3$StableState.id == states[i]] <- colpl[i]
}

obs <- with(SS3, trans3d(PCoA1, PCoA2, energy, p))
pred <- with(SS3, trans3d(PCoA1, PCoA2, fitted(mod2), p))
points(obs, col = colors, pch = 16, cex=0.5)

dev.off()
#######################################################





pdf(file=sprintf('%s/Landscapes_3D_suppl_%s_%s_%s.pdf', dir$figdir, target, th, lv), width=12, height=4)

par(mfrow=c(1,3))
par(oma = c(1, 1, 1, 1)) 
par(mar = c(1, 1, 1, 1)) 

mod1 <- gam(energy ~ s(PCoA1, PCoA2), data= SS3)
vis.gam(mod1, ticktype='detailed', color="cm", theta=20, phi=60, main="Energy")

mod2 <- gam(SSenergy ~ s(PCoA1, PCoA2), data= SS3)
vis.gam(mod2, ticktype='detailed', color="cm", theta=20, phi=60, main="Stable-state energy")

mod3 <- gam(energy.gap ~ s(PCoA1, PCoA2), data= SS3)
vis.gam(mod3, ticktype='detailed', color="cm", theta=20, phi=60, main="Energy gap")


dev.off()
