############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

devtools::install_github("hiroakif93/R-functions/AnalysisHelper", force=TRUE)

set.seed(123)
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr'))

# -- Create directory to save
dir <- make.dir('03_shaping_data/01_output')

# -- Load data table
sml <- readRDS("Table/sample_info_01.rds") ; dim(sml)
head(sml)

write.table(cbind(Sample.ID=rownames(sml), sml), file=sprintf('%s/SampleInfo_2903samples.txt', dir$tabledir), quote=F, sep='\t', row.names=F)


## ========================================================== ##
## -- Make color palette

crop <- sort(unique(na.omit(sml$crop)))
cropcol <- c('firebrick', 'darkgreen', 'darkolivegreen1', 'olivedrab3', 'palegreen',
			  'slateblue1', 'navajowhite1', 'palegreen4', 'lightgreen', 'bisque', 'bisque3',
			  'palevioletred1', 'seashell', 'plum3', 'rosybrown4', 'seagreen4', 'red', 'khaki','tomato')
names(cropcol) <- crop

## ========================================================== ##
## -- Append new column
sml$CNratio <- as.numeric(as.vector(sml $Total_C))/as.numeric(as.vector(sml$Total_N))

# -- Compile longtitude/latitude to decimal sysitem
lat <- do.call(rbind, sapply(sml$latitude_north, strsplit, "\\." )) ; lat2 <- apply(lat, 2, as.numeric)
lon <- do.call(rbind, sapply(sml$longitude_east, strsplit, "\\." )) ;lon2 <- apply(lon, 2, as.numeric)

sml$latitude <- lat2[,1] + lat2[,2]/60 + lat2[,3]/3600
sml$longitude <- lon2[,1] + lon2[,2]/60 + lon2[,3]/3600

sml$site2 <- paste( format(floor(sml$latitude * 10^1) /10^1, nsmall=1),
					format(floor(sml$longitude * 10^1) /10^1, nsmall=1) )

siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

## ========================================================== ##

# -- Remove outlier sample
num <- c()
for(i in 1:ncol(sml)){	if(is.numeric(sml[,i])){num <- c(num,i)}	}

lf <- gather(sml, key=env, value=value, num)

checkgg <- ggplot(lf)+
			 geom_histogram(aes(x = value, fill= crop))+
			 facet_wrap(~env, scales='free')+
			 scale_fill_manual(values= cropcol)+
			 labs(subtitle='Colored by soil.groups') +
			 theme_minimal()
ggsave(plot= checkgg, filename=sprintf('%s/Check_outlier_sample.pdf', dir$figdir), width=21,heigh=29)

rmOutlier <- function(x){ rm=which(x > (mean(x,na.rm=TRUE)+5*sd(x,na.rm=TRUE))); x[rm] <- NA; return(x) }

smltmp <- sml
smltmp[,num] <- apply(smltmp[,num], 2, rmOutlier)

lf <- gather(smltmp, key=env, value=value, num)

checkgg <- ggplot(lf)+
			 geom_histogram(aes(x =value, fill= crop))+
			 facet_wrap(~env, scales='free')+
			 scale_fill_manual(values= cropcol)+
			 labs(subtitle='Colored by soil.groups') +
			 theme_minimal()
ggsave(plot= checkgg, filename=sprintf('%s/Check_outlier_sample2.pdf', dir$figdir), width=21,heigh=29)

# -- Manually removing
ex1 <- which(smltmp$EC_electric_conductivity>20)
smltmp[ex1, 'EC_electric_conductivity'] <- NA

## ========================================================== ##
## -- Check Multi-colinearlinity

cl <- c()
for(j in 1:ncol(smltmp)) cl <- c(cl, class(smltmp[,j]))
nas <- apply(smltmp[,cl=="numeric"], 2, function(x){ sum(is.na(x))/length(x)})
    
expV <- c(names(which(nas < 0.3)))
cormat <- cor(na.omit(smltmp[,expV]), use='pairwise.complete.obs', method='pearson')
tree <- hclust(dist(cormat))

corlf <- gather(data.frame(x=colnames(cormat), cormat), y, value, -1)
corlf$x <- factor(corlf$x, levels=tree$label[tree$order])
corlf$y <- factor(corlf$y, levels=tree$label[tree$order])

corgg <- ggplot(corlf)+
		 geom_tile(aes(x=x, y=y, fill=value))+
		 scale_fill_gradient2()
ggsave(plot=corgg, filename=sprintf('%s/Correlation_between_env.pdf', dir$figdir))

pdf(sprintf('%s/scatter_pairs.pdf', dir$figdir),w=21,h=21)
pairs(smltmp[,expV]); dev.off()

smlomit <- na.omit(smltmp[,expV[c(1,3,4,5,10,11)]]); dim(smlomit)
smlscale <- apply(smlomit, 2, scale)
pca <- prcomp(smlscale)
saveRDS(pca, sprintf('%s/PCA.rds', dir$rdsdir))
sink(sprintf('%s/PCA.txt', dir$tabledir)); print(pca); sink()

## ========================================================== ##

pc <- data.frame(pca$x[,1:5])
rownames(pc) <- rownames(smlomit)

nas <- apply(smltmp, 2, function(x){ sum(is.na(x))/length(x)})
expV <- c(names(which(nas < 0.3)))

smlnew <- cbind(smltmp[,c(expV, 'DSI', 'DiseasePlant')], pc[rownames(smltmp), ])
rownames(smlnew) <- gsub('S_0','S_', rownames(smlnew))

############################################################################
## -- Categorize Disease indices
disease <- smlnew[,c('DSI', 'DiseasePlant')]

## -- disease information
diseaseRow <- apply(disease, 1, function(x) sum(is.na(x)) )

diseaseRmNA <- disease[diseaseRow<2,]
diseaseRmNA[is.na(diseaseRmNA)] <- 0

diseaseLev <- diseaseSum <- data.frame(row.names=names(which(diseaseRow<2)), 
                         disease=rep(0,nrow(diseaseRmNA)) )
anyD <- names(which(diseaseRow==1))
allD <- names(which(diseaseRow==0))

diseaseSum[anyD,1] <- diseaseRmNA[anyD, 1]
diseaseSum[anyD,1] <- diseaseRmNA[anyD, 2]
diseaseSum[allD,1] <- diseaseRmNA[allD, 2]

diseaseLev[which(100 >= diseaseSum[,1] & diseaseSum[,1]>=80), ] <- 'Level 5'
diseaseLev[which(80 > diseaseSum[,1] & diseaseSum[,1]>=60), ] <- 'Level 4'
diseaseLev[which(60 > diseaseSum[,1] & diseaseSum[,1]>=40), ] <- 'Level 3'
diseaseLev[which(40 > diseaseSum[,1] & diseaseSum[,1]>=20), ] <- 'Level 2'
diseaseLev[which(20 > diseaseSum[,1] & diseaseSum[,1]>=0), ] <- 'Level 1'

diseaseLev[diseaseLev==0] <- 'Level 5'

sml <- cbind(smlnew, DL=diseaseLev[rownames(smlnew),1])

############################################################################

write.table(cbind(Sample.ID=rownames(sml), sml), file=sprintf('%s/SampleInfo_RmOutliers.txt', dir$tabledir), quote=F, sep='\t', row.names=F)

saveRDS(sml, 'Table/0301_sml.rds')

############################################################################
