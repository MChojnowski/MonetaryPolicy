 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/42874/OneDrive - Bain/PhD/Workspace/Monetary_policy")
library(devtools)
install_github("Mchojnowski/SentR")
set.seed(2019)

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
library(colorRamps)
library(matlib)
library(tidyverse)
library(SentR)
library(RColorBrewer)

lag <- stats::lag

source("./holidays.R")
source("./PrivateLib.R")


##### PARAMETRY #####
poczatek_proby <- c(1995,1)
poczatek_model <- c(1998,2)

u_oS<-matrix(0,6,6)
u_oS[1,c(2)]<-NA
u_oS[2,c(1)]<-NA
u_oS[3,c(1,2)]<-NA
u_oS[4,c(1,5)]<-NA
u_oS[5,c(2,4)]<-NA
u_oS[6,c(3,4,5)]<-NA
diag(u_oS)<-1

### FUNCTIONS ###
plot_closed_irf<-function(model,title=""){

 	plot.ts(model$model
		,ylim=c(min(unlist(model)),max(unlist(model)))
		,lwd=2
		,main=title)

 	lines(model$calc,col="red",lwd=2,lty=2)

  	lines(model$closed,col="blue",lwd=2)

  	lines(model$Psi,col="orange",lty=4)

	legend("topright"
		,legend=c("vars::irf","calc IRF","IRF closed channel","Psi")
		,lty=c(1,1,1,4)
		,pch=c(NA,NA,NA,NA)
		,col=c("black","red","blue","orange")
		,bg="white"
		,horiz=FALSE
		,xpd=TRUE
		,cex = 0.7)
}

##### DATA PREP #####
#szbsvar
#?szbsvar()
#c(0.7,0.2,2,1000,1,0.1,100)

source("./DataPrep.R")
lSENTY<-list(NULL)
lIS<-list(NULL)
#0.95,0.15,1,0.1584893,1,6.3095734,1584.8931925

#result <- NULL
#for(i in seq(0.2,2,by=0.2)){
sent_org <- getSents(
  dane_var[,-7] # data
  ,3 # lagS
  ,3 # lagB
  ,32 # n_rep
  ,u_oS # Amat
  , # Bmat
  , # exo
  ,c(0.7,0.2,1.6,1,1,10,10) # SZBSVAR
)

qwe <- ts(sent_org$sent,end=end(dane_var),freq=12)

qwe1<-ts(cumsum(qwe)/sd(cumsum(qwe)),freq=12,end=end(dane_var))
qwe1<-(qwe1-mean(qwe1))/sd(qwe1)
#result <- rbind(result,c(i,cor(
#  window(qwe1,start=c(1999,1),end=c(2014,12))
#  ,window(esi_norm,start=c(1999,1),end=c(2014,12))
#)))
#}

dane_monet<-window(dane_monet,start=poczatek_proby)
qwe1<-window(qwe1,start=poczatek_proby)
esi_sent <- window(esi_norm,start=poczatek_proby,end=end(dane_monet))

##### MODEL SEARCH #####
grid_sent <- gridSearch(
    dane_monet # Data
    ,1 # P Higher
    ,1 # P Lower
		,qwe1			# Sentymenty
		,mean(qwe1)			# srednia sentymentu
		,sd(qwe1)			# sd sentymentu
		,0.1			# Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
		,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
		,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
	)

grid_esi <- gridSearch(
  dane_monet # Data
  ,1 # P Higer
  ,1 # P Lower
  ,esi_sent			# Sentymenty
  ,0			# srednia sentymentu
  ,1			# sd sentymentu
  ,0.1			# Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
  ,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
  ,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
)

plotGrid(grid_sent)	
plotGrid(grid_esi)		

prob_sent <- pnorm(qwe1,grid_sent$optparam[2],grid_sent$optparam[1])
prob_esi <- pnorm(esi_sent,grid_esi$optparam[2],grid_esi$optparam[1])

var_sent <- STSVAR(dane_monet,1,1,prob_sent)
var_esi <- STSVAR(dane_monet,1,1,prob_esi)

AmatMP <- diag(4)
AmatMP[lower.tri(AmatMP)] <- NA
BmatMP <- diag(4)
diag(BmatMP) <- NA

model_D<-VAR(dane_monet,p=1)

##### PLOTS IRF #####
zmienne <- colnames(dane_monet)
# Extracted sentiments
for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_sent,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/",zmienne[i],"Shock_SENT.pdf",sep=""))
  dev.off()
}

# Extracted sentiments
for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_esi,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/",zmienne[i],"Shock_ESI.pdf",sep=""))
  dev.off()
}

## IRF 3D for sentiments
IRF3D <- irf3D(var_sent,c("WIBOR","CPI"),c("GDP","CPI","REER"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"WIBOR","CPI",48,100)
dev.print(pdf,paste("./PDF/IRF3D_SENT.pdf",sep=""))
dev.off()

## IRF 3D for ESI
IRF3D <- irf3D(var_esi,c("WIBOR","CPI"),c("GDP","CPI","REER"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"WIBOR","CPI",48,100)
dev.print(pdf,paste("./PDF/IRF3D_ESI.pdf",sep=""))
dev.off()

##### OTHER PLOTS #####
## Probability SENT
dev.new()
plot(prob_esi, ylab="Probability", main="Probability of higher regime - SENT")
dev.print(pdf,paste("./PDF/Probability_SENT.pdf",sep=""))
dev.off()

## Probability ESI
dev.new()
plot(prob_esi, ylab="Probability", main="Probability of higher regime - ESI")
dev.print(pdf,paste("./PDF/Probability_ESI.pdf",sep=""))
dev.off()

## Plot ESI, plot sent

## extracted sentiments
dev.new()
plot(qwe1, ylab="", main="Extracted sentiments")
dev.print(pdf,paste("./PDF/SENT.pdf",sep=""))
dev.off()

## Probability SENT
dev.new()
plot(esi_sent, ylab="Probability", main="Probability of higher regime - ESI")
dev.print(pdf,paste("./PDF/ESI.pdf",sep=""))
dev.off()