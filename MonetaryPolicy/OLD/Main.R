##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/piotr.dybka/Desktop/Mrozu")

##### BIBLIOTEKA #####
library(MSBVAR)
library(vars)
library(forecast)
library(tseries)
library(urca)
library(Hmisc)
library(mFilter)
library(gtrendsR)
library(spikeslab)
library(bsts)
library(lubridate)
library(plyr)
library(parallel)
library(pryr)
library(magrittr)
library(doParallel)
library(iterators)
library(foreach)
library(tsDyn)
library(nnet)
library(reshape)


source("./holidays.R")
source("./PrivateLib.R")

##### PARAMETRY #####
sob<-0
pocz<-as.Date("01-01-2007","%d-%m-%Y")
p_lag<-3
assign("p_lag",p_lag,envir = .GlobalEnv)
psst<-0.4

u_oS<-matrix(0,6,6)
u_oS[1,c(2)]<-NA
u_oS[2,c(1)]<-NA
u_oS[3,c(1,2)]<-NA
u_oS[4,c(1,5)]<-NA
u_oS[5,c(2,4)]<-NA
u_oS[6,c(3,4,5)]<-NA
diag(u_oS)<-NA
##### DATA PREP #####
#szbsvar

source("./DataPrep.R")
lSENTY<-list(NULL)
lIS<-list(NULL)

qwe<-NULL
for (i in 1:32){
	qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7],u_oS,6,12,144,c(0.9,0.1,0.3,1000,1,10,10),12)$sent  #c(lambda..)
		)
}
qwe<-apply(qwe,1,median)
plot(cumsum(qwe)/sd(cumsum(qwe)),type="l",ylim=c(-1,4))
lines(cumsum(sent[-c(1:5)])/sd(cumsum(sent[-c(1:5)])),col="red")

####################################################################################
qwe1<-cumsum(qwe)/sd(cumsum(qwe))
plot(qwe1, type="l")

h<-24 #horyzont œredniej
z<-h-1
me<-seq(length.out=length(qwe1)-z)

for (i in 1:length(me)){
  j<-i+h-1
  me[i]<-mean(qwe1[i:j])
  
}

qwe2<-qwe1[-c(1:z)]

for (i in 1:length(me)){
  qwe2[i]<-qwe2[i]-me[i]
}

plot(qwe2,type="l")



qwe3<-qwe2
for (i in 1:length(qwe2)){
  if (qwe2[i] >0){
    qwe3[i]<-1}
    else {
      qwe3[i]<-0}
}

qwe3
plot(qwe3, type ="l")
plot(qwe4,type ="l")

qweX <- (qwe1 - mean(qwe1))/sd(qwe1)

qwe4 <- (1+exp(qweX))^(-1)
plot(qwe4,type ="l")
write.csv2(qwe1, file="qwa1.csv")
qwel <- read.csv2("mar.csv")
plot(qwel, type="l")
qweX <- (qwel[,2] - mean(qwel[,2]))/sd(qwel[,2])

qwe4 <- (1+exp(qweX))^(-1)
plot(qwe4,type ="l")

####################################################################################
if(FALSE){
#(nrow(dane_var)-59-p_lag)
for (i in 1:(nrow(dane_var)-59-p_lag)){
	print(i)
	ds<-NULL
	for(o in 1:256){
		ds<-cbind(
			ds,
			tryCatch(
				{getSentsB(ts(dane_var[i:(i+59+p_lag),-7],freq=12),u_oS,6,12,2,c(0.95,0.95,1,1000,0.0001,0.00001,10000),12)$sent},
				error=function(a){return(rep(-Inf,58))}
			)
		)
	}
	ds<-apply(ds[,apply(ds,2,min)>-Inf],2,function(x){x/sd(x)})
	lSENTY[[i]]<-apply(ds,1,median)
}
}




pds<-apply(ds,1,median)
rds<-apply(ds-pds,2,rmse)

(ds%*%rds)/sum(rds)

plot(c(ds[,1]),type="l",ylim=c(-8,8))
for (i in 2:ncol(ds)){
	lines(c(ds[,i]))
}
lines(apply(ds,1,median),col="red",lwd=2)


hist(ds[1,])

median(ds[1,])
qwe<-read.csv("./Temp/sentiments.csv")
for(i in 1:(ncol(qwe)-1)){
	lSENTY[[i]]<-qwe[,i+1]

}

dev.new()
par(bg="grey")
plot(c(lSENTY[[1]],rep(NA,nrow(dane_var)-59-p_lag)),type="l",ylim=c(-4,4))
for (i in 2:(nrow(dane_var)-59-p_lag)){
	lines(c(rep(NA,i-1),lSENTY[[i]],rep(NA,nrow(dane_var)-59-p_lag-i+1)),type="l")
}
lines(c(sent/sd(sent)),col="orange",lwd=3)





lIS[[i]]<-InformationSet(
		lSENTY[[i+125]],
		tail(dane_gg[i:(i+59+p_lag-1),],58),
		115,sqrt(3),10e-04,32,
		)
for (i in 2:(nrow(dane_var)-59-p_lag-125)){
	cat(paste("###",i,"/",(nrow(dane_var)-59-p_lag-125),"###"))
	lIS[[i]]<-InformationSet(
				lSENTY[[i+125]],
				tail(dane_gg[i:(i+59+p_lag-1),],58),
				115,sqrt(3),10e-04,4,
				substr(
					rownames(lIS[[i-1]]$model$summary),
					regexpr("]",rownames(lIS[[i-1]]$model$summary))+1,
					nchar(rownames(lIS[[i-1]]$model$summary))
				)
			)
}


for (i in 1:length(lIS)){
Save.InformationSet(
	cbind(
		substr(
			rownames(lIS[[i]]$model$summary),
			regexpr("]",rownames(lIS[[i]]$model$summary))+1,
			nchar(rownames(lIS[[i]]$model$summary))
			),
		lIS[[i]]$model$summary$bma
	),
	paste("./TempTXT/kwerendy_",i,".txt",sep="")
)
}












if(FALSE){
##### FORECAST #####
dzialaj<-c(5,5)
NEW_PROGNOZA<-mc_prognoza2(dane_var,dane_gg_sent,5,60,115,sqrt(3)*100000,TRUE)

NEW_PROGNOZA$prognoza<-cbind(
	NEW_PROGNOZA$prognoza,
	c(tail(dane_var[,2],55)),
	c(tail(dane_var[,2],53),rep(NA,2)),
	c(tail(dane_var[,2],50),rep(NA,5)),
	c(tail(dane_var[,2],44),rep(NA,11))
)
hor<-3:ncol(NEW_PROGNOZA$prognoza)

NEW_PROGNOZA_1M<-NEW_PROGNOZA$prognoza[,c(1,2,hor[(hor %% 4)==3])]
NEW_PROGNOZA_3M<-NEW_PROGNOZA$prognoza[-c((nrow(NEW_PROGNONEW_PROGNOZA_1M<-NEW_PROGNOZA$prognozZA_1M)-1):(nrow(NEW_PROGNOZA_1M))),c(1,2,hor[(hor %% 4)==0])]
NEW_PROGNOZA_6M<-NEW_PROGNOZA$prognoza[-c((nrow(NEW_PROGNOZA_1M)-4):(nrow(NEW_PROGNOZA_1M))),c(1,2,hor[(hor %% 4)==1])]
NEW_PROGNOZA_12M<-NEW_PROGNOZA$prognoza[-c((nrow(NEW_PROGNOZA_1M)-10):(nrow(NEW_PROGNOZA_1M))),c(1,2,hor[(hor %% 4)==2])]

NAZWY_PROGNOZA<-NEW_PROGNOZA$kwerendy
#rm(NEW_PROGNOZA)

NEW_PROGNOZA_1M<-apply(NEW_PROGNOZA_1M,2,as.numeric)
NEW_PROGNOZA_3M<-apply(NEW_PROGNOZA_3M,2,as.numeric)
NEW_PROGNOZA_6M<-apply(NEW_PROGNOZA_6M,2,as.numeric)
NEW_PROGNOZA_12M<-apply(NEW_PROGNOZA_12M,2,as.numeric)

blad_1M_new<-NEW_PROGNOZA_1M[,-ncol(NEW_PROGNOZA_1M)]-NEW_PROGNOZA_1M[,ncol(NEW_PROGNOZA_1M)]
blad_3M_new<-NEW_PROGNOZA_3M[,-ncol(NEW_PROGNOZA_3M)]-NEW_PROGNOZA_3M[,ncol(NEW_PROGNOZA_3M)]
blad_6M_new<-NEW_PROGNOZA_6M[,-ncol(NEW_PROGNOZA_6M)]-NEW_PROGNOZA_6M[,ncol(NEW_PROGNOZA_6M)]
blad_12M_new<-NEW_PROGNOZA_12M[,-ncol(NEW_PROGNOZA_12M)]-NEW_PROGNOZA_12M[,ncol(NEW_PROGNOZA_12M)]

new_blad<-cbind(
	apply(blad_1M_new[,],2,rmse),
	apply(blad_3M_new[,],2,rmse),
	apply(blad_6M_new[,],2,rmse),
	apply(blad_12M_new[,],2,rmse)
)


##### RESULTS #####
colnames(new_blad)<-c("1M","3M","6M","12M")
rownames(new_blad)<-c("naive","RW",
			"ARIMA_00","ARIMA_10","ARIMA_20","ARIMA_01","ARIMA_11","ARIMA_21","ARIMA_02","ARIMA_12","ARIMA_22",
			"ARIMAXS_00","ARIMAXS_10","ARIMAXS_20","ARIMAXS_01","ARIMAXS_11","ARIMAXS_21","ARIMAXS_02","ARIMAXS_12","ARIMAXS_22",
			"ARIMAXW_00","ARIMAXW_10","ARIMAXW_20","ARIMAXW_01","ARIMAXW_11","ARIMAXW_21","ARIMAXW_02","ARIMAXW_12","ARIMAXW_22",
			"LSTAR_11","LSTAR_21","LSTAR_21","LSTAR_22",
			"ARIMAXGT_00","ARIMAXGT_10","ARIMAXGT_20","ARIMAXGT_01","ARIMAXGT_11","ARIMAXGT_21","ARIMAXGT_02","ARIMAXGT_12","ARIMAXGT_22",
			"ARNNETGT_00","ARNNETGT_10","ARNNETGT_20","ARNNETGT_01","ARNNETGT_11","ARNNETGT_21","ARNNETGT_02","ARNNETGT_12","ARNNETGT_22",
			"NNET_00","NNET_10","NNET_20","NNET_01","NNET_11","NNET_21","NNET_02","NNET_12","NNET_22",
			"NNETGT_00","NNETGT_10","NNETGT_20","NNETGT_01","NNETGT_11","NNETGT_21","NNETGT_02","NNETGT_12","NNETT_22",
			"VAR_1","VAR_2","VAR_3"
	)

rownames(new_blad)<-paste(1:nrow(new_blad),rownames(new_blad),sep="_")

### In respect to RW
(zzz<-t(apply(new_blad,1,function(x){x*100/new_blad[2,]})))
rbind(
	c("ARMA",apply(zzz[3:11,],2,function(x){sum(x<100)}),"||",apply(zzz[3:11,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXS",apply(zzz[12:20,],2,function(x){sum(x<100)}),"||",apply(zzz[12:20,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXW",apply(zzz[21:29,],2,function(x){sum(x<100)}),"||",apply(zzz[21:29,],2,function(x){round(mean(x[x<100]),2)})),
	c("LSTAR",apply(zzz[30:33,],2,function(x){sum(x<100)}),"||",apply(zzz[30:33,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXGT",apply(zzz[34:42,],2,function(x){sum(x<100)}),"||",apply(zzz[34:42,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){sum(x<100)}),"||",apply(zzz[43:51,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNET",apply(zzz[52:60,],2,function(x){sum(x<100)}),"||",apply(zzz[52:60,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){sum(x<100)}),"||",apply(zzz[61:69,],2,function(x){round(mean(x[x<100]),2)})),
	c("VAR",apply(zzz[70:72,],2,function(x){sum(x<100)}),"||",apply(zzz[70:72,],2,function(x){round(mean(x[x<100]),2)}))
)


### In respect to best ARMA
(zzz<-t(apply(new_blad,1,function(x){x*100/apply(new_blad[3:11,],2,min)})))
rbind(
	c("ARMAXS",apply(zzz[12:20,],2,function(x){sum(x<100)}),"||",apply(zzz[12:20,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXW",apply(zzz[21:29,],2,function(x){sum(x<100)}),"||",apply(zzz[21:29,],2,function(x){round(mean(x[x<100]),2)})),
	c("LSTAR",apply(zzz[30:33,],2,function(x){sum(x<100)}),"||",apply(zzz[30:33,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXGT",apply(zzz[34:42,],2,function(x){sum(x<100)}),"||",apply(zzz[34:42,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){sum(x<100)}),"||",apply(zzz[43:51,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNET",apply(zzz[52:60,],2,function(x){sum(x<100)}),"||",apply(zzz[52:60,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){sum(x<100)}),"||",apply(zzz[61:69,],2,function(x){round(mean(x[x<100]),2)})),
	c("VAR",apply(zzz[70:72,],2,function(x){sum(x<100)}),"||",apply(zzz[70:72,],2,function(x){round(mean(x[x<100]),2)}))
)


### In respect to best WWUK
(zzz<-t(apply(new_blad,1,function(x){x*100/apply(new_blad[21:29,],2,min)})))
rbind(
	c("ARMAXS",apply(zzz[12:20,],2,function(x){sum(x<100)}),"||",apply(zzz[12:20,],2,function(x){round(mean(x[x<100]),2)})),
	c("LSTAR",apply(zzz[30:33,],2,function(x){sum(x<100)}),"||",apply(zzz[30:33,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARMAXGT",apply(zzz[34:42,],2,function(x){sum(x<100)}),"||",apply(zzz[34:42,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){sum(x<100)}),"||",apply(zzz[43:51,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNET",apply(zzz[52:60,],2,function(x){sum(x<100)}),"||",apply(zzz[52:60,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){sum(x<100)}),"||",apply(zzz[61:69,],2,function(x){round(mean(x[x<100]),2)})),
	c("VAR",apply(zzz[70:72,],2,function(x){sum(x<100)}),"||",apply(zzz[70:72,],2,function(x){round(mean(x[x<100]),2)}))
)

### In respect to best Linear Sent (ARMAXS & ARMAXGT)
(zzz<-t(apply(new_blad,1,function(x){x*100/apply(new_blad[c(12:20,34:42),],2,min)})))
rbind(
	c("LSTAR",apply(zzz[30:33,],2,function(x){sum(x<100)}),"||",apply(zzz[30:33,],2,function(x){round(mean(x[x<100]),2)})),
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){sum(x<100)}),"||",apply(zzz[43:51,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){sum(x<100)}),"||",apply(zzz[61:69,],2,function(x){round(mean(x[x<100]),2)}))
)

### In respect to best VAR
(zzz<-t(apply(new_blad,1,function(x){x*100/apply(new_blad[61:63,],2,min)})))
rbind(
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){sum(x<100)}),"||",apply(zzz[43:51,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNET",apply(zzz[52:60,],2,function(x){sum(x<100)}),"||",apply(zzz[52:60,],2,function(x){round(mean(x[x<100]),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){sum(x<100)}),"||",apply(zzz[61:69,],2,function(x){round(mean(x[x<100]),2)}))
)

### BEST among the BEST
(zzz<-t(apply(new_blad,1,function(x){x*100/new_blad[2,]})))
rbind(
	c("ARMA",apply(zzz[3:11,],2,function(x){round(min(x),2)})),
	c("ARMAXS",apply(zzz[12:20,],2,function(x){round(min(x),2)})),
	c("ARMAXW",apply(zzz[21:29,],2,function(x){round(min(x),2)})),
	c("LSTAR",apply(zzz[30:33,],2,function(x){round(min(x),2)})),
	c("ARMAXGT",apply(zzz[34:42,],2,function(x){round(min(x),2)})),
	c("ARNNETGT",apply(zzz[43:51,],2,function(x){round(min(x),2)})),
	c("NNET",apply(zzz[52:60,],2,function(x){round(min(x),2)})),
	c("NNETGT",apply(zzz[61:69,],2,function(x){round(min(x),2)})),
	c("VAR",apply(zzz[70:72,],2,function(x){round(min(x),2)}))
)
rbind(
	c("ARMA",round(apply(zzz[3:11,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("ARMAXS",round(apply(zzz[12:20,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("ARMAXW",round(apply(zzz[21:29,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("LSTAR",round(apply(zzz[30:33,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("ARMAXGT",round(apply(zzz[34:42,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("ARNNETGT",round(apply(zzz[43:51,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("NNET",round(apply(zzz[52:60,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("NNETGT",round(apply(zzz[61:69,],2,min)/apply(zzz[3:11,],2,min),2)),
	c("VAR",round(apply(zzz[70:72,],2,min)/apply(zzz[3:11,],2,min),2))
)

##### PLOTS #####
dev.new()
plot(NEW_PROGNOZA$pers,type="l",ylim=c(0,1))
###
dev.new(width=36,height=36)
par(mfrow=c(2,1))
qwe<-NEW_PROGNOZA$senty

ord<-NULL
plot(ts(qwe[1,][-c(1:3,57:59)],start=c(2008,3),freq=12),ylim=c(min(qwe[!is.na(qwe)]),max(qwe[!is.na(qwe)])),type="l",ylab="Sentiment",main="Sentiments for each rolling window draw")
for(q in 2:nrow(qwe)){
	lines(ts(qwe[q,][-c((q):(q+2),(56+q):(58+q))],start=c(2008,3),freq=12))
	ord<-rbind(ord,auto.arima(ts(qwe[q,][-c((q):(q+2),(56+q):(58+q))]a,start=c(2008,3),freq=12))$arma)
}

plot(qwe[1,],ylim=c(min(qwe[!is.na(qwe)]),max(qwe[!is.na(qwe)])),type="l")
for(q in 2:nrow(qwe)){
	lines(qwe[q,])
}
###
NAZWY<-list(NULL)
NAZWY_print<-matrix(0,nrow=ncol(dane_gg),ncol=length(3:57))
rownames(NAZWY_print)<-colnames(dane_gg)

persystencja<-NULL
for (p in 1:55){
	NAZWY[[p]]<-as.character(read.table(paste("./TempTXT/Kwerendy_",p+2,".txt",sep=""))[,1])
	NAZWY_print[,p]<-as.numeric(colnames(dane_gg) %in% NAZWY[[p]])
	
	if(p!=1){
		persystencja<-c(persystencja,sum(NAZWY_print[,p-1]+NAZWY_print[,p]==2)/sum(NAZWY_print[,p-1]+NAZWY_print[,p]>0))	
	}
}

popularne<-apply(NAZWY_print,1,sum)
popularne<-popularne[order(-popularne)]
popularne<-popularne[popularne>0]
head(popularne,15)

Np<-NAZWY_print[rownames(NAZWY_print) %in% names(popularne),]
Np<-Np[order(apply(Np,1,sum)),]

cols <- c(
'0' = "#FFFFFF",
'1' = "#FF0000"
)

dev.new(width=45,height=45)
image(1:ncol(Np),1:nrow(Np),t(Np),col=cols,main=paste("Persistency:",psst),xlab="Period",ylab="")
dev.new()
plot(ts(c(rep(NA,3),persystencja)),ylab="Pers",main="Persistency in given period",ylim=c(0,1))
abline(h=psst,lwd=2)
abline(h=mean(persystencja),lty=4)
dev.new()
plot(ts(c(NA,NA,apply(Np,2,sum))),ylab="# of queries",main="# of Google Queries")
}
