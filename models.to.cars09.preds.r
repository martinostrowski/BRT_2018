#####
# This script loads a set of BRT models, one by one, and runs gbm predictions against each set of 
# monthly CARS 2009 surface data frames. The predictions are read into a a list of monthly 
# matrices. The matrices can be converted into a cars matrix one column at a time. Code for 
# plotting the relabunds of individual zOTUs is included. 
#
# The monthly matrices will be combined and analysed by cluster analyses to identify 
# ecologically relevant clusters
#####

library(oce)
library(rgdal)
library(gbm)
library(dismo)
source("brt.functionsImageP2018.R")
library(maps)
library(ncdf4)

zotus1000<-read.csv("~/data/brt_201807/tz_models_210618_DRV/input.table.3568.nona.csv", h=T)

#sort models by total zOTU abundance in the input data

#write.csv(names(rev(sort(colSums(zotus1000[,21:3577])))), "~/data/brt_201807//model.names.by.order.csv", quote=F)

zotu.names<-names(rev(sort(colSums(zotus1000[,21:3577]))))

#make a list of potential models

model.filenames<-vector('double', length(zotu.names))

for (i in 1:length(model.names)){
	model.filenames[i]<-paste(zotu.names[i], ".gaussian.RData", sep="")
}

models<-list.files("~/data/brt_201807/tz_models_210618_DRV/", pattern="AMD*")

#check which models are actually available

model.set<-model.filenames[model.filenames %in% models]


#make a load list
model.filenames<-vector('double', length(model.set))
for (i in 1:length(model.set)){

	model.filenames[i]<-paste("~/data/brt_201807/tz_models_210618_DRV/",model.set[i], sep="")
}

nmods<-length(model.filenames)


########

d0cars.filenames<-list.files(pattern='cars.matrixd0')
d0cars.files<-vector('list', 12)

for (i in 1:12){

  d0cars.files[[i]]<-read.csv(paste("~/data/brt_201807/CARS/",d0cars.filenames[[i]], sep=""), h=F)
  colnames(d0cars.files[[i]])<-c('lat','lon', 'month', 'depth','temp', 'sal', 'nox', 'po4', 'sil', 'n.p','DL','strat') 
	
} 

cars09<-nc_open("cars2009.nc")

clat<-ncvar_get(cars09, 'LATITUDE')
clon<-ncvar_get(cars09, 'LONGITUDE')

#make a list of prediction matrices

preds.list<-vector('list',12)

for (m in 1:12){

  preds.list[[m]]<-matrix(-9999, nrow=238651, ncol=nmods)

}



######

pb <- txtProgressBar(min = 0, max = 500, style = 3)
k<-0

zotu.names<-vector('double', nmods)
sp.names<-vector('double', nmods)


cv<-vector('double', nmods)
cvse<-vector('double', nmods)
ntrees<-vector('double', nmods)
dev<-vector('double', nmods)
contributions<-matrix(nrow=9, ncol=nmods)

rownames(contributions)<-c("temp","sal","depth","strat","nox", "po4", "n.p", "DL", "sil")


for (l in 1002:1499){

	load(model.filenames[l])

	zotu.names[l]<-my.mod[[28]][[5]][[1]]
	
	cv[l]<-my.mod$cv.statistics$correlation.mean
  	cvse[l]<-my.mod$cv.statistics$correlation.se
  	dev[l]<-my.mod$cv.statistics$deviance
  	ntrees[i]<-my.mod$n.trees
  	contributions[,i]<-my.mod[[32]][c("temp","sal","depth","strat","nox", "po4", "n.p", "DL", "sil"),2]
  
  	tmp.result<-cbind(zotu.names, cv, cvse, dev, ntrees, t(contributions))
  	write.table(tmp.result[(tmp.result[,"ntrees"]!=0),], file="evaluation.data.set.csv", sep=',', quote=F)
	
	

		for (n in 1:12){
  		
  			preds.list[[n]][,l]<-gbm.predict.grids(my.mod, d0cars.files[[n]], want.grids=T, pred.vec=rep(-9999,238651), filepath="~/", num.row=331, num.col=721, xll=0, yll=0, cell.size=100, no.data=-9999, plot=F)
		}

k<-k+1
setTxtProgressBar(pb, k)
	
	}


month<-c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12")

for (n in 1:12){

write.csv(preds.list[[n]], file=paste("cars.", month[n], "zotus.gaussian.201807.csv", sep=','), quote=F) 

}





