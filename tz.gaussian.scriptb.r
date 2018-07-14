#This script is the backbone of a set of processes to generate BRT models from a dataframe of 1000s of zOTU abundance profiles.

#The script loads the required libraries, imports the input table, imports the datasets to predict against,
#set up a matrix to receive each set of predictions, a list to store the species names as 
#they are extracted from the model and a set of vectors to store the evaluation data
#and contributions.

#Recent additions include an error catching trycatch that allows the loop to continue when
#gbm.step throws an error

#To do
#Split this up into a batch file or incorporate a parrallel process. 
#The models cannot be stored in memory because they are too large

library(tidyverse)
library(dismo)
library(gbm)

#input.table.3568.nona.csv", h=T, sep=',')
zotus1000<-read.csv("input.table.sqrt.nona.csv", h=T)

#zotus.batch<-sqrt(zotus1000[21:3577])
#zotus1000<-cbind(zotus1000[,c(1:20)], zotus.batch)

#phb<-read.table("PH_Historical_Env_through_2017_mo.txt", sep='\t', h=T)

phb<-read.csv("ph.historical.nona.mo.csv", h=T)
gs<-read.csv("gs.line.12months.csv", h=T)
tz<-read.csv("tzline.csv", h=T)

#phb$n.p<-phb$nox/phb$po4
#colnames(phb)<-c("Sample", "Day" ,"Month" ,"Year","Bottom.depth","Lat","Lon" , "depth","DL","strat","temp","sal","sil", "nox", "po4")
#write.table(phb, "ph.historical.nona.mo.csv", sep=',', quote=F)



nzotu<-ncol(zotus1000)-21
my.zotus<-vector('list', nzotu)
for (i in 1:nzotu){
my.zotus[[i]]<-colnames(zotus1000)[i+21]
}

sp.names<-vector('double', nzotu)
preds.phb<-matrix(nrow=nrow(phb), ncol=nzotu)
preds.tz<-matrix(nrow=nrow(tz), ncol=nzotu)
preds.gs<-matrix(nrow=nrow(gs), ncol=nzotu)
#preds.cars<-matrix(nrow(cars), ncol=nzotu)

#my.mods<-vector('list', nzotu)

eval.data.out<-vector('double', nzotu)
cv<-vector('double', nzotu)
cvse<-vector('double', nzotu)
ntrees<-vector('double', nzotu)
dev<-vector('double', nzotu)
contributions<-matrix(nrow=8, ncol=nzotu)
rownames(contributions)<-c("temp","sal","depth","strat","nox", "po4", "n.p", "DL", "sil")




for (i in 2821:2899){

my.mod<-try(update(gbm.step(data=zotus1000,gbm.x = c(13:21), gbm.y = i+21, family = "gaussian", tree.complexity = 10, learning.rate = 0.001, bag.fraction = 0.5)), TRUE)

if(isTRUE(class(my.mod)=="try-error")) { next } else {


	sp.names[i]<-my.mod[[28]][5][[1]]

	#preds.mai[,i]<-predict.gbm(my.mod, mai[,c(7,10:12,8:9)], n.trees=my.mod$gbm.call$best.trees, type="response")

	preds.phb[,i]<-predict.gbm(my.mod, phb, n.trees=my.mod$gbm.call$best.trees, type="response")

	#preds.rot[,i]<-predict.gbm(my.mod, rot[,c(7,10:12,8:9)], n.trees=my.mod$gbm.call$best.trees, type="response")

	#preds.cars[,i]<-predict.gbm(my.mod, cars, n.trees=my.mod$gbm.call$best.trees, type="response")

	preds.tz[,i]<-predict.gbm(my.mod, tz, n.trees=my.mod$gbm.call$best.trees, type="response")

	#preds.gs[,i]<-predict.gbm(my.mod, gs, n.trees=my.mod$gbm.call$best.trees, type="response")

  eval.data.out<-unlist(sp.names)
  cv[i]<-my.mod$cv.statistics$correlation.mean
  cvse[i]<-my.mod$cv.statistics$correlation.se
  ntrees[i]<-my.mod$n.trees
  contributions[,i]<-my.mod[[32]][c("temp","sal","depth","strat","nox", "po4", "n.p", "DL"),2]
  
  tmp.result<-cbind(eval.data.out, cv, cvse, ntrees, t(contributions))
  write.table(tmp.result[(tmp.result[,"ntrees"]!=0),], sep='\t', quote=F)
  save.image(paste(my.mod[[28]][5][[1]], ".gaussian.RData", sep=""))
rm(my.mod)

colnames(preds.tz)<-eval.data.out
colnames(preds.phb)<-eval.data.out
write.table(preds.tz[,colSums(is.na(preds.tz))==0], "preds.tz.13.tmp", sep=',', quote=F)
write.table(preds.phb[,colSums(is.na(preds.phb))==0], "preds.phb.13.tmp", sep=',', quote=F)
 }
}

write.table(preds.tz[,colSums(is.na(preds.tz))==0], "preds.tz.5.tmp", sep=',', quote=F)
write.table(preds.phb[,colSums(is.na(preds.phb))==0], "preds.phb.5.tmp", sep=',', quote=F)
write.table(preds.phb, file=paste(batch.number, "phb.g.preds", sep=""), sep=',', quote=F)
write.table(preds.gs, file=paste(batch.number, "gs.g.preds", sep=""), sep=',', quote=F)
write.table(preds.tz, file=paste(batch.number, "ts.g.preds", sep=""), sep=',', quote=F)
write.table(preds.cars, file=paste(batch.number, "cars.g.preds", sep=""), sep=',', quote=F)
save.image(paste(batch.number, "gaussian.RData", sep=""))



write.table(preds.tz[,colSums(is.na(preds.tz))==0], "preds.tz.5.tmp", sep=',', quote=F)
write.table(preds.phb[,colSums(is.na(preds.phb))==0], "preds.phb.5.tmp", sep=',', quote=F)
write.table(preds.phb, file=paste(batch.number, "phb.g.preds", sep=""), sep=',', quote=F)
write.table(preds.gs, file=paste(batch.number, "gs.g.preds", sep=""), sep=',', quote=F)
write.table(preds.tz, file=paste(batch.number, "ts.g.preds", sep=""), sep=',', quote=F)
write.table(preds.cars, file=paste(batch.number, "cars.g.preds", sep=""), sep=',', quote=F)
save.image(paste(batch.number, "gaussian.RData", sep=""))