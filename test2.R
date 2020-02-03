 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/42874/OneDrive - Bain/PhD/Workspace/Monetary_policy")

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
library(plotly)
library(readr)

source("./holidays.R")
source("./PrivateLib.R")
source("./GridSearch.R")

##### PARAMETRY #####
p_lag<-2
assign("p_lag",p_lag,envir = .GlobalEnv)

u_oS<-matrix(0,6,6)
u_oS[1,c(2)]<-NA
u_oS[2,c(1)]<-NA
u_oS[3,c(1,2)]<-NA
u_oS[4,c(1,5)]<-NA
u_oS[5,c(2,4)]<-NA
u_oS[6,c(3,4,5)]<-NA
diag(u_oS)<-NA

n_loop<-128
##### DATA PREP #####
#szbsvar
#?szbsvar()

source("./DataPrep.R")
wwuk<-read.csv("./Dane/ESIPL.csv")
wwuk<-ts(wwuk[,2],start=c(1993,10),freq=12)

lSENTY<-list(NULL)
lIS<-list(NULL)

##### SENTIMENTS #####
#lm.szbsvar<-c(0.5,0.5,1,100,1,100,100) 
lm.szbsvar<-c(0.95,0.15,1,0.15,1,6.3,1600) 
lm.szbsvar<-c(0.92,0.15,1,0.175,1,6.7,1600)

#range_1<-seq(0.05,1,by=0.05)
#range_2<-seq(0.25,5,by=0.25)
#range_3<-10^tail(seq(-4,4,length.out=21),-1)

range_1<-seq(0.89,0.99,length.out=11)
range_2<-seq(0.05,1,length.out=39)
range_3<-seq(1000,2000,length.out=11)
range_4<-seq(0.05,0.3,length.out=11)
range_5<-seq(6,7,length.out=11)

lm.range<-cbind(
		 range_2
		,range_2
		)

lm.order<-2

for(lm.step in lm.order){
	stahp<-NULL
	for(jj in lm.range[,lm.step]){
		lm.szbsvar.temp<-lm.szbsvar
		lm.szbsvar.temp[lm.step]<-jj
		qwe<-NULL	
		for (i in 1:n_loop){
			qwe<-cbind(
				qwe,
				getSentsB(dane_var[,-7] 		#Dane
					,u_oS 				#Macierz A
					,6 					#Lagi w VARze
					,12 					#Lag w BVARze
					,1 					# szlif
					,lm.szbsvar.temp			#lambda
					,12 					#Okresowosc
					)$sent  				#c(lambda..)
			)

		}
	qwe<-window(qwe,end=c(2017,12))
	qwe<-apply(qwe,1,median)
	stahp<-c(stahp,cor(-cumsum(qwe),window(wwuk,start=c(1994,8),end=c(2017,12))))
	}
	lm.szbsvar[lm.step]<-lm.range[,lm.step][which.max(stahp)]
}

print(lm.szbsvar)

qwe<-NULL	
for (i in 1:256){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1 					# szlif
			,lm.szbsvar				#lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)















plot_qwe<-ts((cumsum(qwe)-mean(cumsum(qwe)))/sd(cumsum(qwe)),end=c(2017,12),freq=12)
#plot_qwe<-ts(read.csv("./sents.csv")[,2],end=c(year(enddate),month(enddate)),freq=12)
plot_wwuk<-ts((wwuk-mean(wwuk))/sd(wwuk),start=start(wwuk),freq=12)

dev.new()
plot(0,type="n",xlim=c(1990,2018),ylim=c(-5,3),lwd=2,ylab="Normalized values",xlab="Time")
rect(1990,-7,2020,5,col="#EEEEEE")
#abline(h=0,col="#555555")
lines(plot_wwuk,col="black",lwd=2,lty=4)
#lines(plot_sent,col="red",lwd=2)
lines(plot_qwe,col="red",lwd=2)
legend("bottomright",legend=c("ESI","Extracted sentiment"),lty=c(4,1),col=c("black","red"),bg="white",horiz=FALSE,xpd=TRUE, cex = 1.2)
dev.print(pdf,"./PDF/SENTY.pdf")
dev.off()