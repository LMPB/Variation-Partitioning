library(ggplot2)
library(spaa)
library(reshape2)
library(adespatial)
library(spdep)

######################################
###Beta-diversity Partitioning Step###
######################################
spe <- read.csv("C:/The/way/to/your/data/asv_table_cleared_rarefied.csv", header=TRUE, row.names = 1, sep =",")
spe <- t(spe)
spe_rarefy<-rrarefy(spe,min(rowSums(spe)))
##as we were using taxonimic inference as ASVs, the rarer organisms need to be cut. 
#For more details see the DADA2 documentation
spe.cut<-spe_rarefy[,colSums(spe_rarefy) > 10] 
dim(spe.freq); sum(spe.freq)
dim(spe.cut);sum(spe.cut)
spe.pa<- decostand(spe.cut,"pa")

#Beta-diversity partition
#The traditional betadiversity partition uses presence-absence data
#We also used the betapair abund using Bray-Curtis
beta.pair=beta.pair(spe.pa, index.family = "jaccard")

nestedness <- data.frame(t(combn(rownames(spe.pa),2)), as.numeric(beta.pair$beta.jne))
names(nestedness) <- c("S1", "S2", "nestedness")
turnover <- data.frame(t(combn(rownames(spe.pa),2)), as.numeric(beta.pair$beta.jtu))
names(turnover) <- c("S1", "S2", "turnover")
dissimilarity <- data.frame(t(combn(rownames(spe.pa),2)), as.numeric(beta.pair$beta.jac))
names(dissimilarity) <- c("S1", "S2", "dissimilarity")

beta.pair.abund=beta.pair.abund(spe.cut, index.family = "bray")
abund.variation <- data.frame(t(combn(rownames(spe.cut),2)), as.numeric(beta.pair.abund$beta.bray.bal))
names(abund.variation) <- c("S1", "S2", "variation")
abund.gradient <- data.frame(t(combn(rownames(spe.cut),2)), as.numeric(beta.pair.abund$beta.bray.gra))
names(abund.gradient) <- c("S1", "S2", "gradient")
abund.dissimilarity <- data.frame(t(combn(rownames(spe.cut),2)), as.numeric(beta.pair.abund$beta.bray))
names(abund.dissimilarity) <- c("S1", "S2", "dissimilarity")


#Boxplot
total <- merge (dissimilarity,turnover,by = c("S1","S2"), all = T)
complete <- merge (total,nestedness,by = c("S1","S2"), all = T)

total_abund <- merge (abund.variation,abund.gradient,by = c("S1","S2"), all = T)
complete_abund <- merge (total_abund,abund.dissimilarity,by = c("S1","S2"), all = T)

#Change "complete" to "complete_abund" will show the partition of abundance variation
complete %>% 
  gather(key = "betadiv", value = "value", -S1, -S2) %>% 
  transform(betadiv=factor(betadiv,levels=c("dissimilarity","turnover","nestedness"))) %>% 
  ggplot(aes(x=betadiv,y=value, fill=betadiv))+
  geom_boxplot(outlier.shape = NA, width = 0.5,size=0.5, fill = "white")+
  geom_jitter(width = 0.2,alpha=0.2)+
  #scale_color_manual(name="", labels=c("dissimilarity","variation","gradient"), 
  #                  values = c("black","black","black"))+
  scale_fill_manual(name="", labels=c("dissimilarity","turnover","nestedness"), 
                    values = c("Grey90","Grey70","Grey50"))+
  theme_bw()+
  #stat_summary(fun=mean,fill="black", colour="black", geom="point", 
  #             shape=18, size=3,show.legend = FALSE) +
  #geom_text(data = means, aes(label = mean, y = mean), color="black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+
  theme(legend.position = "none")

#Dotplot
par(mar=rep(2,4))
plot(beta.pair$beta.sim, beta.pair$beta.sor)
text(beta.pair$beta.sim,cex=0.6, pos=3, col="black")

#plot Dissimilarity by distance
library(fields)

#geolocation set as long-lat
#numerical form
geo <- as.matrix(read.csv("C:/The/way/to/your/data/geography.csv", header=TRUE, row.names = 1, sep =","))
#Obtain a geographic distance matrix
geo.mat <- rdist.earth(geo, geo, miles = FALSE, R = NULL)
geo.mat[upper.tri(geo.mat, diag=TRUE)]<- NA
#
rownames(spe.cut)==rownames(geo.mat)
identical(rownames(spe.cut),rownames(geo.mat))

geo.dist<-na.omit(subset(melt(geo.mat)))
colnames(geo.dist)<- c("S1","S2","geo_dist")
geo.dist<-arrange(geo.dist, S1, S2)

complete<-arrange(complete, S1, S2)
identical(rownames(complete),rownames(geo.dist))
final<-cbind(complete,geo.dist)
final<-final[,-6]

ggplot(final,aes(x=geo_dist,y=dissimilarity))+
  geom_jitter ()+
  #xlim(-20,15000)+
  #geom_hline(yintercept = -2, linetype = "dashed")+
  #geom_hline(yintercept = 2, linetype = "dashed")+
  geom_smooth(method = "lm")+
  theme_classic()+
  labs(x= paste(" ",sep="\n"),color = " ")