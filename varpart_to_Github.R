library(vegan)
library(betapart)
library(adespatial)

#Variation Partitioning Step
#import data
#ASVs table: rows are asvs, columns are sites
#Environmental table: rows are sites, environmental variables are columns
#the environmental table should be normalyzed
#As we were using taxonimic inference as ASVs, the rarer organisms need to be cut. 
#For more details see the DADA2 documentation
spe <- read.csv("C:/The/way/to/your/data/asv_table_cleared_rarefied.csv", header=TRUE, row.names = 1, sep =",")
spe <- t(spe)
spe_rarefy<-rrarefy(spe,min(rowSums(spe)))
spe.cut<-spe[,colSums(spe) > 10]
dim(spe.freq); sum(spe.freq)
dim(spe.cut);sum(spe.cut)
spe.pa<- decostand(spe.cut,"pa")
#Here, the total dissimilarity recovered in the betadiversity partition will be used as the biological data
beta.pair=beta.pair(spe.pa, index.family = "jaccard") #Presence-absence analysis
beta.pair.abund=beta.pair.abund(spe.cut, index.family = "bray") #Abundance analysis
#From this point, just change "beta.pair" to "beta.pair.abund" if you want to analyse abundances

env <- read.csv("C:/The/way/to/your/data/environmental_table.csv",header=TRUE, row.names = 1, sep =",")
denv=decostand(env, na.rm=T, method="standardize")
denv[,2]<-env[,2] #pH is already normalized

geo <- read.csv("C:/The/way/to/your/data/geographic_table.csv", header=TRUE, row.names = 1, sep =",")
rownames(spe.pa)==rownames(denv)
identical(rownames(spe.pa),rownames(denv))
identical(rownames(geo),rownames(denv))
identical(rownames(spe.pa),rownames(geo))

time <- read.csv("C:/The/way/to/your/data/collection_day_table.csv", header=TRUE, row.names = 1, sep =",")

#Forward selection of environmental data
rda2env <- capscale(beta.pair$beta.jac ~ ., data=denv,distance='jaccard') #distance = "bray" if using beta.pair.abund
rda1env <- capscale(beta.pair$beta.jac ~ 1, data=denv, distance='jaccard') #distance = "bray" if using beta.pair.abund
step.forward <- ordiR2step(rda1env,scope=formula(rda2env), direction="forward", perm.max=200,pstep=999)
step.forward$anova
#Parsimonious subsets of explanatory variables (based on #forward selection)
names(denv)
env.fraction <- denv[, c(1,2,6)] #Chose significant results from "step.forward" selection
head(env.fraction)

#Forward selection of geographic data
#http://www.hiercourse.com/docs/variationPartitioning.html

#SPACE FRACTION
geo.dbmem <- dbmem(geo, store.listw = TRUE,MEM.autocor = "positive")
if(require("adegraphics", quietly = TRUE)){
  s.value(geo,geo.dbmem[,1:9])
  plot(geo.dbmem[,1:6], geo, pSp.cex = 3)
}
str(scores(geo.dbmem))
ordisurf(geo,scores(geo.dbmem,choi=1),bubble=4, main='dbMEM 1')

cap.dbmem <- capscale(beta.pair$beta.jac ~ ., data=as.data.frame(scores(geo.dbmem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
mod0.dbmem <- capscale(beta.pair$beta.jac ~ 1, data=as.data.frame(scores(geo.dbmem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
step.dbmem <- ordiR2step(mod0.dbmem, scope=formula(cap.dbmem), direction="forward", perm.max=200,pstep=999)
step.dbmem$anova
space.fraction <- scores(geo.dbmem, choices=c(13,23)) #Chose significant results from "step.forward" selection


#NEIGHBORHOOD ANALYSIS
#Performing MEM (Moran I Eningenvectors Analysis)
nbtri <- tri2nb(geo) #make triangulation considering distances
dist.tri <- nbdists(nbtri, as.matrix(geo))
fdist.tri <- lapply(dist.tri, function(x) 1 - x/max(dist(as.matrix(geo))))
listw.tri <- nb2listw(nbtri, glist = fdist.tri, style = "W")

geo.mem <- scores.listw(nb2listw(nbtri, listw.tri$weights), MEM.autocor = "positive")
if(require("adegraphics", quietly = TRUE)){
  s.value(geo,geo.mem[,1:9])
  plot(geo.mem[,1:6], geo, pSp.cex = 3)
}

str(scores(geo.mem))
ordisurf(geo,scores(geo.mem,choi=1),bubble=4, main='dbMEM 1')

cap.mem <- capscale(beta.pair$beta.jac ~ ., data=as.data.frame(scores(geo.mem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
mod0.mem <- capscale(beta.pair$beta.jac ~ 1, data=as.data.frame(scores(geo.mem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
step.mem <- ordiR2step(mod0.mem, scope=formula(cap.mem), direction="forward", perm.max=200,pstep=999)
step.mem$anova
neighborhood.fraction <- scores(geo.mem, choices=c(11,15)) #Chose significant results from "step.forward" selection


#REGIONAL ANALYSIS
#Performing MEM (Moran I Eningenvectors Analysis)
nbtri <- tri2nb(geo) #make triangulation considering distances
l<- read.csv("C:/The/way/to/your/data/weight_matrix.csv", header=TRUE, row.names = 1, sep =",") #Weight matrix, site by site showing the connection level each pair of sites should have
l.dist<-as.dist(l)
listw.l<-mat2listw(as.matrix(l.dist))

reg.mem <- scores.listw(nb2listw(listw.l$neighbours, listw.l$weights, style = "W"),MEM.autocor = "positive")
if(require("adegraphics", quietly = TRUE)){
  s.value(geo,reg.mem[,c(1,2,3,4)])
  plot(reg.mem[,c(1,2,3,4)], geo, pSp.cex = 3)
}

str(scores(reg.mem))
ordisurf(geo,scores(reg.mem,choi=4),bubble=4,main='dbMEM 4')

cap.mem <- capscale(beta.pair$beta.jac ~ ., data=as.data.frame(scores(reg.mem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
mod0.mem <- capscale(beta.pair$beta.jac ~ 1, data=as.data.frame(scores(reg.mem)), distance='jaccard') #distance = "bray" if using beta.pair.abund
step.mem <- ordiR2step(mod0.mem, scope=formula(cap.mem), direction="forward", perm.max=200,pstep=999)
step.mem$anova
regional.fraction <- scores(reg.mem, choices=c(4)) #Chose significant results from "step.forward" selection


#CONNECTIVITY ANALYSIS
#ADD a citação aqui
connect<-read.csv("C:/The/way/to/your/data/connections_matrix.csv", header=TRUE, row.names = 1, sep =",") #Weight matrix, site by river branch, set values as "1" if the water flows from the specific branch in direction to the site
(connect.aem <- aem(binary.mat = as.matrix(connect)))
connect.aem.vec <- connect.aem$vectors

if(require("adegraphics", quietly = TRUE)){
  s.value(geo,connect.aem.vec[,1:3])
  plot(connect.aem.vec[,1:3], geo, pSp.cex = 3)
}

str(scores(connect.aem.vec))
ordisurf(geo,scores(connect.aem.vec,choi=3),bubble=4,main='AEM 3')

cap.mem <- capscale(beta.pair$beta.jac ~ ., data=as.data.frame(scores(connect.aem.vec)), distance='jaccard') #distance = "bray" if using beta.pair.abund
mod0.mem <- capscale(beta.pair$beta.jac ~ 1, data=as.data.frame(scores(connect.aem.vec)), distance='jaccard') #distance = "bray" if using beta.pair.abund
step.mem <- ordiR2step(mod0.mem, scope=formula(cap.mem), direction="forward", perm.max=200,pstep=999)
step.mem
step.mem$anova
connectivity.fraction <- scores(connect.aem.vec, choices=c(23,13)) #Chose significant results from "step.forward" selection

#Variation Partitioning
#par(mfrow=c(2,4))
var.part <- varpart(beta.pair$beta.jacy,env.fraction,space.fraction,time) #change to the fractions you want to compare
plot(var.part, digits=2, Xnames = c('Environment', 'Space', 'time'))
var.part
showvarparts(3)
#Test of all testable fractions
#test of fractions [a+b+c+d+e+f+g]
all.pars <-cbind(space.fraction,neighborhood.fraction, time)
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(all.pars)),step=1000)

#test of fractions [a+d+f+g] => Space fraction
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(space.fraction)),step=1000)
#test of fractions [b+d+b+g] => Neighborhood fraction
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(neighborhood.fraction)),step=1000)
#test of fractions [a+d+b+g] => Regional fraction
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(time)),step=1000)

#test of fractions [a+d+b+g] => Space vs Neighborhood
all.pars <-cbind(space.fraction,neighborhood.fraction)
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(all.pars)+ Condition(as.matrix(time)),sqrt.dist = T),step=1000)
#test of fractions [b+e+c+g] => Neighborhood vs Region
all.pars <-cbind(space.fraction,time)
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(all.pars)+ Condition(as.matrix(neighborhood.fraction)),sqrt.dist = T),step=1000)
#test of fraction [a+f+c+g] => Space vs Region
all.pars <-cbind(neighborhood.fraction,time)
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(all.pars)+ Condition(as.matrix(space.fraction)),sqrt.dist = T),step=1000)

#test of fraction [a]
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(space.fraction)+ Condition(as.matrix(neighborhood.fraction))+ Condition(as.matrix(time)),sqrt.dist = T),step=1000)
#test of fraction [b]
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(neighborhood.fraction)+ Condition(as.matrix(space.fraction))+ Condition(as.matrix(time)),sqrt.dist = T),step=1000)
#test of fraction [c]
anova.cca(dbrda(beta.pair$beta.jac~as.matrix(time)+ Condition(as.matrix(space.fraction))+ Condition(as.matrix(neighborhood.fraction)),sqrt.dist = T),step=1000)