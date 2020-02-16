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
library(colorRamps)

source("./holidays.R")
source("./PrivateLib.R")
source("../../SentiR/GridSearch.R")

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

##### DATA PREP #####
#szbsvar
#?szbsvar()
#c(0.7,0.2,2,1000,1,0.1,100)

source("./DataPrep.R")
lSENTY<-list(NULL)
lIS<-list(NULL)
#0.95,0.15,1,0.1584893,1,6.3095734,1584.8931925

if(TRUE){
qwe<-NULL
for (i in 1:1024){
	qwe<-cbind(
		qwe,
		getSentsB(dane_var[,-7] 			#Dane
			,u_oS 					#Macierz A
			,6 						#Lagi w VARze
			,12 						#Lag w BVARze
			,1 						# ?
			,c(0.92,0.125,1,0.175,1,6.7,1600)	#lambda
			,12 						#Okresowosc
			)$sent  					#c(lambda..)
		)
}

qwe<-apply(qwe,1,median)
qwe1<-ts(cumsum(qwe)/sd(cumsum(qwe)),freq=12,end=end(dane_var))
plot(cumsum(qwe)/sd(cumsum(qwe)),type="l",ylim=c(-1,4))

lines(lm(cumsum(qwe)/sd(cumsum(qwe))~c(1:length(qwe)))$fitted,col="red")

if(FALSE){
	asd<-cumsum(qwe)/sd(cumsum(qwe))
	asd<-c(
		forecast(arima(rev(asd),order=c(3,0,0),method="CSS"),h=3)$mean
		,asd
		,forecast(arima(asd,order=c(3,0,0),method="CSS"),h=3)$mean
		)
	
	dev.new()
	plot(hpfilter(ts(asd,freq=12)))
	dev.print(pdf,paste("./PDF/SentOrg.pdf",sep=""))
	dev.off()

	qwe1<-hpfilter(ts(asd,freq=12))$cycle
	#qwe1<-cumsum(qwe)/sd(cumsum(qwe))-lm(cumsum(qwe)/sd(cumsum(qwe))~c(1:length(qwe)))$fitted
	qwe1<-ts(head(tail(qwe1,-3),-3),freq=12,end=end(dane_var))
}

}
######### for external value of sentiments
if(FALSE){
	qwe1<-read.csv("sentfull2.csv")
	qwe1<-qwe1[,2]
	qwe1<-ts(qwe1[-c(1:3)],freq=12,end=end(dane_var))
}

if(FALSE){
	qwe1<-read.csv("./Dane/ESIPL.csv")
	qwe1<-qwe1[,2]
	qwe1<-ts(qwe1,freq=12,start=c(1993,10))
}

qwe1<-(qwe1-mean(qwe1))/sd(qwe1)

###### Exponential smoothing
qwe_sm <- holt(qwe1, alpha=0.1,  beta=FALSE, initial="simple", h=12) 
fit<-fitted(qwe_sm)
#plot(fit)

###### liczymy ró¿nice
qwe2<-fit-qwe1
#plot(qwe2, type="l")

sdp<-sd(qwe2)

###### liczymy prawdopodobieñstwa
#prob1<-pnorm(qwe2,  mean = 0, sd = sdp)

#z<-length(wynagrodzenie)-length(prob)

if(TRUE){
	start.dane<-c(1995,1)
	pkb<-window(pkb,start=start.dane)
	cpi<-window(cpi,start=start.dane)
	wibor<-window(wibor,start=start.dane)
	reer<-window(reer,start=start.dane)
	qwe1<-window(qwe1,start=start.dane)
}

dane_monet<-na.omit(
	cbind(
  		diff(pkb)+1,
  		diff(cpi)+1,
  		diff(wibor)+1,
  		diff(reer)+1
	)
)


colnames(dane_monet)<-c("GDP","CPI","WIBOR","REER")
var_dflt<-VAR(dane_monet,p=p_lag)

#####################################################################################
###### Przekszta³cenie danych do modelu VAR     #####################################
################################################################################################################################

#################### PreGrid ######################
resDfltVar<-residuals(var_dflt)
	resDfltVar[,1]<-resDfltVar[,1]/apply(dane_monet,2,sd)[1]
	resDfltVar[,2]<-resDfltVar[,2]/apply(dane_monet,2,sd)[2]
	resDfltVar[,3]<-resDfltVar[,3]/apply(dane_monet,2,sd)[3]
	resDfltVar[,4]<-resDfltVar[,4]/apply(dane_monet,2,sd)[4]


if(FALSE){
dev.new()
par(mfrow=c(2,4))

#################### SD1 ########################
gs<-gridSearch(var_dflt,qwe1,"sd",0,1)

#bestSdFct<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(apply(gs$aic[gs$roots[,1]<1 & gs$roots[,2]<1,],1,sum))]
bestSdFct<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1])]

if(identical(bestSdFct,numeric(0))){bestSdFct<-sd(qwe1)}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,mean)/abs(AIC(var_dflt))
	,col=kolor
	,pch=19
	,main=paste("sd#1 || min:",bestSdFct)
	,ylab="Relative SSE"
	)

abline(h=1)

#################### MEAN1 ########################
gs<-gridSearch(var_dflt,qwe1,"mean",,bestSdFct)

bestMeanFct<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(gs$sse[gs$roots[,1]<1 & gs$roots[,2]<1])]
#bestMeanFct<-gs$grid[apply(gs$roots,1,max)<1][which.min(apply(gs$aic[apply(gs$roots,1,max)<1,],1,sum))]

if(identical(bestMeanFct,numeric(0))){bestMeanFct<-mean(qwe1)}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("mean#1 || min:",round(bestMeanFct,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### SD2 ########################
gs<-gridSearch(var_dflt,qwe1,"sd",bestMeanFct,)

#bestSdFct2<-gs$grid[apply(gs$roots,1,max)<1][which.min(apply(gs$aic,1,sum)[apply(gs$roots,1,max)<1])]
bestSdFct<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1])]

if(identical(bestSdFct2,numeric(0))){bestSdFct2<-bestSdFct}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,mean)
	,col=kolor
	,pch=19
	,main=paste("sd#2 || min:",round(bestSdFct2,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### MEAN2 ########################
gs<-gridSearch(var_dflt,qwe1,"mean",,bestSdFct2)

bestMeanFct2<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(gs$sse[gs$roots[,1]<1 & gs$roots[,2]<1])]
#bestMeanFct2<-gs$grid[apply(gs$roots,1,max)<1][which.min(apply(gs$aic,1,sum)[apply(gs$roots,1,max)<1])]

if(identical(bestMeanFct2,numeric(0))){bestMeanFct2<-bestMeanFct}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("mean#2 || min:",round(bestMeanFct2,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### SD3 ########################
gs<-gridSearch(var_dflt,qwe1,"sd",bestMeanFct2,)

#bestSdFct3<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(apply(gs$aic,1,sum)[gs$roots[,1]<1 & gs$roots[,2]<1])]
bestSdFct3<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1])]

if(identical(bestSdFct3,numeric(0))){bestSdFct3<-bestSdFct2}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("sd#3 || min:",round(bestSdFct3,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### MEAN3 ########################
gs<-gridSearch(var_dflt,qwe1,"mean",,bestSdFct3)

#bestMeanFct3<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(apply(gs$aic,1,sum)[gs$roots[,1]<1 & gs$roots[,2]<1])]
bestMeanFct3<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1])]

if(identical(bestSdFct,numeric(0))){bestMeanFct3<-bestMeanFct2}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("mean#3 || min:",round(bestMeanFct3,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### SD4 ########################
gs<-gridSearch(var_dflt,qwe1,"sd",bestMeanFct3,)

#bestSdFct4<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(apply(gs$aic,1,sum)[gs$roots[,1]<1 & gs$roots[,2]<1])]
bestSdFct4<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1])]

if(identical(bestSdFct4,numeric(0))){bestSdFct4<-bestSdFct3}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("sd#4 || min:",round(bestSdFct4,5))
	,ylab="Relative SSE"
	)

abline(h=1)

#################### MEAN4 ########################
gs<-gridSearch(var_dflt,qwe1,"mean",,bestSdFct4)

#bestMeanFct4<-gs$grid[gs$roots[,1]<1 & gs$roots[,2]<1][which.min(apply(gs$aic,1,sum)[gs$roots[,1]<1 & gs$roots[,2]<1])]
bestMeanFct4<-gs$grid[apply(gs$roots,1,max)<1][which.min(gs$sse[apply(gs$roots,1,max)<1)]

if(identical(bestSdFct,numeric(0))){bestMeanFct4<-bestMeanFct3}

kolor=rep("#FF0000",length(gs$sse))
kolor[apply(gs$roots,1,max)<1]<-"#FFFF00"
kolor[apply(gs$roots,1,max)<0.95]]<-"#00FF00"

plot(gs$grid
	,gs$sse/sqrt(sum(resDfltVar^2))
	#,apply(gs$aic,1,sum)
	,col=kolor
	,pch=19
	,main=paste("mean#4 || min:",round(bestMeanFct4,5))
	,ylab="Relative SSE"
	)

abline(h=1)
}

###########################################################################################
gs2<-gridSearch2D(var_dflt,qwe1,"all",0,1,60,250)

dev.new()
filled.contour(gs2$gridS
	,gs2$gridM
	,as.matrix(gs2$sse/sqrt(sum(resDfltVar^2)))
	#,as.matrix(gs2$aic/abs(AIC(var_dflt))/2)
	,xlab="Sd"
	,ylab="Mean"
	,color.palette = matlab.like2
	,main="GridSearch"
	,nlevels=48
	)
 
dev.print(pdf,paste("./PDF/GridAll.pdf",sep=""))
dev.off()


gs2$sse[gs2$roots>1]<-Inf

dev.new()
filled.contour(gs2$gridS
	,gs2$gridM
	,as.matrix(gs2$sse/sqrt(sum(resDfltVar^2)))
	#,as.matrix(gs2$aic/abs(AIC(var_dflt))/2)
	,xlab="Sd"
	,ylab="Mean"
	,color.palette = matlab.like2
	,main="GridSearch"
	,nlevels=24
	)

cord<-which(gs2$sse == min(gs2$sse), arr.ind = TRUE)
points(gs2$gridS[cord[1]],gs2$gridM[cord[2]],pch=19,col="orange")

#################### PRINT ########################
dev.print(pdf,paste("./PDF/GridDomain.pdf",sep=""))
dev.off()

#################################################################################################################

######################################################################
###### Model VAR 

#prob1<-abs(pnorm(qwe1,mean(qwe1),sd(qwe1)*bestSdFct)-0.5)*2

#prob<-pnorm(qwe1,bestMeanFct4,bestSdFct4)
prob<-pnorm(qwe1,gs2$gridM[cord[2]],gs2$gridS[cord[1]])

dev.new()
par(mfrow=c(1,1))
plot(prob, type="l")
dev.print(pdf,paste("./PDF/Probabilities.pdf",sep=""))
dev.off()


### SVAR MATRICES ###
u_mon_xS<-matrix(0,4,4)
u_mon_xS[lower.tri(u_mon_xS)]<-NA
u_mon_xS[1,1]<-NA
u_mon_xS[2,2]<-NA
u_mon_xS[3,3]<-NA
u_mon_xS[4,4]<-NA
diag(u_mon_xS)<-1
u_mon_xS

var1<-var2<-var_dflt

dane_monet1<-dane_monet*prob
dane_monet2<-dane_monet*(1-prob)

dane_all<-dane_monet
for (i in 1:var_dflt$p){
  dane_all<-cbind(dane_all,lag(dane_monet1,-i))
}
for (i in 1:var_dflt$p){
  dane_all<-cbind(dane_all,lag(dane_monet2,-i))
}

colnames(dane_monet1)<-c("GDP","CPI","WIBOR","REER")
colnames(dane_monet2)<-c("GDP","CPI","WIBOR","REER")
colnames(dane_monet)<-c("GDP","CPI","WIBOR","REER")


r1<-lm(dane_all[,1]~dane_all[,-c(1:4)]-1)
r2<-lm(dane_all[,2]~dane_all[,-c(1:4)]-1)
r3<-lm(dane_all[,3]~dane_all[,-c(1:4)]-1)
r4<-lm(dane_all[,4]~dane_all[,-c(1:4)]-1)

comat1<-rbind(r1$coef[c(1:(4*var_dflt$p))],r2$coef[c(1:(4*var_dflt$p))],r3$coef[c(1:(4*var_dflt$p))],r4$coef[c(1:(4*var_dflt$p))])
comat2<-rbind(r1$coef[c((4*var_dflt$p+1):length(r1$coef))],r2$coef[c((4*var_dflt$p+1):length(r1$coef))],r3$coef[c((4*var_dflt$p+1):length(r1$coef))],r4$coef[c((4*var_dflt$p+1):length(r1$coef))])

for (i in 1:4){
  nazwy<-names(var1$varresult[[i]]$coefficients)
  var1$varresult[[i]]$coefficients<-c(comat1[i,],0)  
  names(var1$varresult[[i]]$coefficients)<-nazwy
  
  nazwy<-names(var2$varresult[[i]]$coefficients)
  var2$varresult[[i]]$coefficients<-c(comat2[i,],0)
  names(var2$varresult[[i]]$coefficients)<-nazwy
  
}

### DODAJ CONSTANT
# wyci¹gamy zmienne 
#dane_all<-na.omit(dane_all)
#dane<-dane_all[,-c(1:4)]
dane<-var_dflt$datamat[,-c(1:4,(var_dflt$p+1)*var_dflt$K+1)]
#dane<-dane[-c(1:var_dflt$p),]
#dane<-na.omit(dane)

#przeliczamy reszty

  for (i in 1:4){
    var1$varresult[[i]]$fitted.values<-var1$varresult[[i]]$fitted.values+tail(var1$varresult[[i]]$coefficients,1)
    var1$varresult[[i]]$residuals<-tail(dane_monet,-var_dflt$p)[,i]-var1$varresult[[i]]$fitted.values
    }
  for (i in 1:4){
    var2$varresult[[i]]$fitted.values<-var2$varresult[[i]]$fitted.values+tail(var2$varresult[[i]]$coefficients,1)
    var2$varresult[[i]]$residuals<-tail(dane_monet,-var_dflt$p)[,i]-var2$varresult[[i]]$fitted.values
  }



#Obliczamy stale
c1<-residuals(var1)
const1<-apply(c1,2,mean)
const1

c2<-residuals(var2)
const2<-apply(c2,2,mean)
const2


### podmieniamy sta³e
for (i in 1:4){
  var1$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const1[i]
}

for (i in 1:4){
  var2$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const2[i]
}



  #przeliczamy reszty
  for (i in 1:4){
    var1$varresult[[i]]$fitted.values<-var1$varresult[[i]]$fitted.values+tail(var1$varresult[[i]]$coefficients,1)
    var1$varresult[[i]]$residuals<-tail(dane_monet,-var_dflt$p)[,i]-var1$varresult[[i]]$fitted.values
    }
  for (i in 1:4){
    var2$varresult[[i]]$fitted.values<- var2$varresult[[i]]$fitted.values+tail(var2$varresult[[i]]$coefficients,1)
    var2$varresult[[i]]$residuals<-tail(dane_monet,-var_dflt$p)[,i]-var2$varresult[[i]]$fitted.values
  }


summary(var1)
summary(var2)

#to jest tylko sprawdzenie, ze sie nie pomylilismy
c1<-residuals(var1)
const1<-apply(c1,2,mean)
const1

c2<-residuals(var2)
const2<-apply(c2,2,mean)
const2

#Tworzymy nowe ladne modele
model_var_1<-var1
model_var_2<-var2

MacB<-matrix(0,nrow=4,ncol=4)
diag(MacB)<-NA

model_svar_1<-SVAR(model_var_1,Amat=u_mon_xS,Bmat=MacB,lrtest = FALSE, max.iter=1000)
model_svar_2<-SVAR(model_var_2,Amat=u_mon_xS,Bmat=MacB,lrtest = FALSE, max.iter=1000)
model_svar_dflt<-SVAR(var_dflt,Amat=u_mon_xS,Bmat=MacB,lrtest = FALSE, max.iter=1000)


### IRFY ###
h<-24 	#ustawiamy horyzont zabawy
ci<-0.95 	#szerokoœæ CI
r<-100 	#liczba iteracji do bootstrappu

##### Zmiany CPI
for(i in 1:4){

	b<-c("GDP","CPI","WIBOR","REER")[i]
	a<-c("GDP","CPI","WIBOR","REER")[-i]

	dev.new()
		par(mfcol=c(3,3))
		for (i in 1:3){
			#High regime
    			var_plot=irf(model_svar_1
				,impulse = a[i]
				,response=b	
				,n.ahead = h
				,ortho=TRUE
				,boot=TRUE
				,runs=r
				,ci=ci
			)
    			hc<-h+1
    			plot(x=c(1:hc)
				,y=unlist(var_plot$Lower)
				,type="l"	
				,lwd = 3
				,lty=2
				,col="red"
				,ylab=paste("Response of ",b)
				,xlab=a[i]
				,ylim=range(c(unlist(var_plot$Lower)
				,unlist(var_plot$Upper)))
				,main="High Regime"
			)
    			lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3,lty=2,col="red")
    			lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
    			abline(a = NULL, h = 0)

			#Low regime
    			var_plot=irf(model_svar_2
				,impulse = a[i]
				,response=b	
				,n.ahead = h
				,ortho=TRUE
				,boot=TRUE
				,runs=r
				,ci=ci
			)
    			hc<-h+1
    			plot(x=c(1:hc)
				,y=unlist(var_plot$Lower)
				,type="l"	
				,lwd = 3
				,lty=2
				,col="red"
				,ylab=paste("Response of ",b)
				,xlab=a[i]
				,ylim=range(c(unlist(var_plot$Lower)
				,unlist(var_plot$Upper)))
				,main="Low regime"
			)
    			lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3,lty=2,col="red")
    			lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
    			abline(a = NULL, h = 0)

			#Defalut VAR
    			var_plot=irf(model_svar_dflt
				,impulse = a[i]
				,response=b	
				,n.ahead = h
				,ortho=TRUE
				,boot=TRUE
				,runs=r
				,ci=ci
			)
    			hc<-h+1
    			plot(x=c(1:hc)
				,y=unlist(var_plot$Lower)
				,type="l"	
				,lwd = 3
				,lty=2
				,col="red"
				,ylab=paste("Response of ",b)
				,xlab=a[i]
				,ylim=range(c(unlist(var_plot$Lower)
				,unlist(var_plot$Upper)))
				,main="No Regime"
			)
    			lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3,lty=2,col="red")
    			lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
    			abline(a = NULL, h = 0)
		}

	dev.print(pdf,paste("./PDF/",b,".pdf",sep=""))
	dev.off()
}

#################### 3D IRF #########################################
MEGAIRFCPI<-NULL
MEGAIRFGDP<-NULL
MEGAIRFREER<-NULL

for(p in seq(0,1,by=0.01)){
	model_var_irf<-var1

	for (i in 1:4){
		model_var_irf$varresult[[i]]$coefficients<-p*var1$varresult[[i]]$coefficients+(1-p)*var2$varresult[[i]]$coefficients
    		model_var_irf$varresult[[i]]$fitted.values<-as.matrix(dane)%*%as.matrix(head(model_var_irf$varresult[[i]]$coefficients,-1))+tail(model_var_irf$varresult[[i]]$coefficients,1)
    		model_var_irf$varresult[[i]]$residuals<-tail(dane_monet,-var_dflt$p)[,i]-model_var_irf$varresult[[i]]$fitted.values
    	}

	model_svar_irf<-SVAR(model_var_irf,Amat=u_mon_xS,Bmat=MacB,lrtest = FALSE, max.iter=1000)
	var_plot<-irf(model_svar_irf, impulse = "WIBOR", response=c("CPI","GDP","REER"), n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
   	MEGAIRFCPI<-cbind(MEGAIRFCPI,var_plot$irf$WIBOR[,"CPI"])
	MEGAIRFGDP<-cbind(MEGAIRFGDP,var_plot$irf$WIBOR[,"GDP"])
	MEGAIRFREER<-cbind(MEGAIRFREER,var_plot$irf$WIBOR[,"REER"])
}

poziomy<-48

############### CPI ############################
zakres<-seq(range(as.matrix(MEGAIRFCPI))[1],range(as.matrix(MEGAIRFCPI))[2],length.out=1000)
paleta<-c(colorRampPalette(colors = c("green", "white"), space="Lab")(round(length(zakres[zakres>0])/1000*poziomy))
	,colorRampPalette(colors = c("white", "red"), space="Lab")(poziomy-round(length(zakres[zakres>0])/1000*poziomy))   
	)   

dev.new()
filled.contour(0:h
	,seq(0,1,by=0.01)
	,as.matrix(MEGAIRFCPI)
	,xlab="Horizon"
	,ylab="Probability of higher regime"
	,col = rev(paleta)
	,main="WIBOR shock on CPI"
	,levels = seq(range(as.matrix(MEGAIRFCPI))[1],range(as.matrix(MEGAIRFCPI))[2], length.out=(poziomy+1))
	,nlevels = poziom
	)
dev.print(pdf,"./PDF/WIBORonCPI3D.pdf")
dev.off()

############### GDP ############################
zakres<-seq(range(as.matrix(MEGAIRFGDP))[1],range(as.matrix(MEGAIRFGDP))[2],length.out=1000)
paleta<-c(colorRampPalette(colors = c("green", "white"), space="Lab")(round(length(zakres[zakres>0])/1000*poziomy))
	,colorRampPalette(colors = c("white", "red"), space="Lab")(poziomy-round(length(zakres[zakres>0])/1000*poziomy))   
	)   

dev.new()
filled.contour(0:h
	,seq(0,1,by=0.01)
	,as.matrix(MEGAIRFGDP)
	,xlab="Horizon"
	,ylab="Probability of higher regime"
	,col = rev(paleta)
	,main="WIBOR shock on GDP"
	,levels = seq(range(as.matrix(MEGAIRFGDP))[1],range(as.matrix(MEGAIRFGDP))[2], length.out=(poziomy+1))
	,nlevels = poziom
	)
dev.print(pdf,"./PDF/WIBORonGDP3D.pdf")
dev.off()

############### REER ############################
zakres<-seq(range(as.matrix(MEGAIRFREER))[1],range(as.matrix(MEGAIRFREER))[2],length.out=1000)
paleta<-c(colorRampPalette(colors = c("green", "white"), space="Lab")(round(length(zakres[zakres>0])/1000*poziomy))
	,colorRampPalette(colors = c("white", "red"), space="Lab")(poziomy-round(length(zakres[zakres>0])/1000*poziomy))   
	)   

dev.new()
filled.contour(0:h
	,seq(0,1,by=0.01)
	,as.matrix(MEGAIRFREER)
	,xlab="Horizon"
	,ylab="Probability of higher regime"
	,col = rev(paleta)
	,main="WIBOR shock on REER"
	,levels = seq(range(as.matrix(MEGAIRFREER))[1],range(as.matrix(MEGAIRFREER))[2], length.out=(poziomy+1))
	,nlevels = poziom
	)
dev.print(pdf,"./PDF/WIBORonREER3D.pdf")
dev.off()

