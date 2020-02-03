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
u_oS[1,c(2,3)]<-NA
u_oS[2,c(1,3)]<-NA
u_oS[3,c(1,2)]<-NA
u_oS[4,c(1,5,6)]<-NA
u_oS[5,c(2,4,6)]<-NA
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

#### lambda 1
stahp<-NULL
for(jj in seq(0.05,0.95,by=0.05)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,2 					# szlif
			,c(0.15,jj,0.8,0.1,0.1,0.1,0.001)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}
optL1<-seq(0.05,0.95,by=0.05)[which.max(stahp)]

#### lambda 2
stahp<-NULL
for(jj in seq(0.25,5,by=0.25)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1 					# szlif
			,c(0.15,optL1,jj,0.1,0.1,0.1,0.001)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}

optL2<-seq(0.25,5,by=0.25)[which.max(stahp)]

#### lambda 4
stahp<-NULL
for(jj in 10^seq(-3,3,by=1)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1					# szlif
			,c(0.15,optL1,optL2,jj,0.1,0.1,0.001)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}
optL4<-10^seq(-3,3,by=1)[which.max(stahp)]


#### lambda 5
stahp<-NULL
for(jj in 10^seq(-3,3,by=1)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1					# szlif
			,c(0.15,optL1,optL2,optL4,jj,0.1,0.001)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}
optL5<-10^seq(-3,3,by=1)[which.max(stahp)]


#### mu0
stahp<-NULL
for(jj in 10^seq(-3,3,by=1)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1					# szlif
			,c(0.15,optL1,optL2,optL4,optL5,jj,0.001)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}
optM0<-10^seq(-3,3,by=1)[which.max(stahp)]


#### mu1
stahp<-NULL
for(jj in 10^seq(-3,3,by=1)){
qwe<-NULL	
for (i in 1:n_loop){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1					# szlif
			,c(0.15,optL1,optL2,optL4,optL5,optM0,jj)   #lambda
			,12 					#Okresowosc
			)$sent  				#c(lambda..)
		)

}
qwe<-window(qwe,end=c(2017,12))
qwe<-apply(qwe,1,median)
stahp<-c(stahp,cor(qwe,window(diff(wwuk),start=c(1994,8),end=c(2017,12))))
}
optM1<-10^seq(-3,3,by=1)[which.max(stahp)]


optL1
optL2
optL4
optL5
optM0
optM1

qwe<-NULL	
for (i in 1:256){
qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 		#Dane
			,u_oS 				#Macierz A
			,6 					#Lagi w VARze
			,12 					#Lag w BVARze
			,1 					# szlif
			,c(0.15,optL1,optL2,optL4,optL5,optM0,optM1)   #lambda
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