 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/pd50409/Desktop/Monetary_policy")

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
library(matlib)



source("./holidays.R")
source("./PrivateLib.R")
source("./GridSearch.R")

##### PARAMETRY #####
p_lag_high<-3
p_lag_low<-4
assign("p_lag_high",p_lag_high,envir = .GlobalEnv)
assign("p_lag_low",p_lag_low,envir = .GlobalEnv)

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
for (i in 1:32){
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

###### liczymy r�nice
qwe2<-fit-qwe1
#plot(qwe2, type="l")

sdp<-sd(qwe2)

###### liczymy prawdopodobie�stwa
#prob1<-pnorm(qwe2,  mean = 0, sd = sdp)

#z<-length(wynagrodzenie)-length(prob)

start.dane<-c(1993,1)
dane_monet<-window(dane_monet,start=start.dane)
qwe1<-window(qwe1,start=start.dane)

var_dflt<-VAR(dane_monet,p=1)
var_dflt_high<-VAR(dane_monet,p=p_lag_high)
var_dflt_low<-VAR(dane_monet,p=p_lag_low)

#####################################################################################
###### Przekszta�cenie danych do modelu VAR     #####################################
################################################################################################################################

#################### PreGrid ######################
resDfltVar<-residuals(var_dflt)
	resDfltVar[,1]<-resDfltVar[,1]/apply(dane_monet,2,sd)[1]
	resDfltVar[,2]<-resDfltVar[,2]/apply(dane_monet,2,sd)[2]
	resDfltVar[,3]<-resDfltVar[,3]/apply(dane_monet,2,sd)[3]
	resDfltVar[,4]<-resDfltVar[,4]/apply(dane_monet,2,sd)[4]

###########################################################################################
## Liczba [0,1] - % bazy danych; Liczba >1 - liczba rekord�w

gs2<-gridSearch2D(
		var_dflt_high 	# Podstawka pod VAR w wyzszym rezimie
		,var_dflt_low 	# Podstawka pod VAR w nizszym rezimie
		,qwe1			# Sentymenty
		,"all"		# Typ szukania (zawsze all)
		,0			# srednia sentymentu
		,1			# sd sentymentu
		,0.75			# Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
		,50			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
		,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
	)

dev.new()
filled.contour(gs2$gridS
	,gs2$gridM
	,as.matrix(gs2$sse/sqrt(sum(resDfltVar^2)))
	#,as.matrix(gs2$aic/abs(AIC(var_dflt))/2)
	#,as.matrix(gs2$roots)
	#,as.matrix(gs2$irf[[1]])
	,xlab="Sd"
	,ylab="Mean"
	,color.palette = matlab.like2
	#,color.palette=terrain.colors
	,main="GridSearch"
	,nlevels=96
	)
 
dev.print(pdf,paste("./PDF/GridAll.pdf",sep=""))
dev.off()

gs2$sse[gs2$roots>=1]<-Inf
#gs2$sse[gs2$gridM<(-1),]<-Inf
#gs2$roots[gs2$roots>=1]<-Inf
#gs2$sse[gs2$sse>=1.1*sqrt(sum(resDfltVar^2))]<-Inf
cord<-which(gs2$sse == min(gs2$sse), arr.ind = TRUE)

dev.new()
filled.contour(gs2$gridS
	,gs2$gridM
	,as.matrix(gs2$sse/sqrt(sum(resDfltVar^2)))
	#,as.matrix(gs2$aic/abs(AIC(var_dflt))/2)
	#,as.matrix(gs2$roots)
	,xlab="Sd"
	,ylab="Mean"
	,color.palette = matlab.like2
	#,color.palette=terrain.colors
	,main="GridSearch"
	,nlevels=24
	,plot.axes={points(gs2$gridS[cord[1]],gs2$gridM[cord[2]],pch=19,col="orange");axis(1);axis(2)}
	)

gs2$sse[cord]/sqrt(sum(resDfltVar^2))

#cord<-which(gs2$roots == min(gs2$roots), arr.ind = TRUE)


#################### PRINT ########################
dev.print(pdf,paste("./PDF/GridDomain.pdf",sep=""))
dev.off()


dev.new()
	plot(qwe1,ylim=c(min(qwe1,qnorm(0.01,gs2$gridM[cord[2]],gs2$gridS[cord[1]])),max(qwe1,qnorm(0.99,gs2$gridM[cord[2]],gs2$gridS[cord[1]]))))
	abline(h=gs2$gridM[cord[2]],col="blue")
	abline(h=qnorm(0.99,gs2$gridM[cord[2]],gs2$gridS[cord[1]]),col="red",lty=4)
	abline(h=qnorm(0.01,gs2$gridM[cord[2]],gs2$gridS[cord[1]]),col="red",lty=4)
	dev.print(pdf,paste("./PDF/SentReg.pdf",sep=""))
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

var1<-var_dflt_high
var2<-var_dflt_low

# Uy+Zy (Z=U-D -> D=U-Z)
		dane_monet1<-dane_monet*prob
		dane_monet2<-dane_monet*(1-prob)

		colnames(dane_monet1)<-c("GDP","CPI","WIBOR","REER")
		colnames(dane_monet2)<-c("GDP","CPI","WIBOR","REER")

		dane_all<-var_dflt_high$y
		for (i in 1:var_dflt_high$p){
  			dane_all<-cbind(dane_all,lag(dane_monet1,-i))
		}
		for (i in 1:var_dflt_low$p){
 			dane_all<-cbind(dane_all,lag(dane_monet2,-i))
		}

		r1<-lm(dane_all[,1]~dane_all[,-c(1:4)])
		r2<-lm(dane_all[,2]~dane_all[,-c(1:4)])
		r3<-lm(dane_all[,3]~dane_all[,-c(1:4)])
		r4<-lm(dane_all[,4]~dane_all[,-c(1:4)])

		
		comat1<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c(2:(4*var_dflt_high$p+1))]
		comat2<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c((4*var_dflt_high$p+2):length(r1$coef))]

		#comat1<-cbind(comat1,matrix(0,nrow=4,ncol=4*max(var_dflt_high$p,var_dflt_low$p)))[,1:(4*max(var_dflt_high$p,var_dflt_low$p))]
		#comat2<-cbind(comat2,matrix(0,nrow=4,ncol=4*max(var_dflt_high$p,var_dflt_low$p)))[,1:(4*max(var_dflt_high$p,var_dflt_low$p))]

		#comat2<-comat1-comat2
		#comat2<-comat2[,1:(4*var_dflt_low$p)]

		#comat1<-comat1[,apply(comat1,2,sd)!=0]
		#comat2<-comat2[,apply(comat2,2,sd)!=0]

		for (i in 1:4){
  			nazwy<-names(var1$varresult[[i]]$coefficients)
  			var1$varresult[[i]]$coefficients<-c(comat1[i,],0)  
  			names(var1$varresult[[i]]$coefficients)<-nazwy
  
  			nazwy<-names(var2$varresult[[i]]$coefficients)
  			var2$varresult[[i]]$coefficients<-c(comat2[i,],0)
  			names(var2$varresult[[i]]$coefficients)<-nazwy
  		}

		### DODAJ CONSTANT
		dane_h<-var_dflt_high$datamat[,-c(1:4,ncol(var_dflt_high$datamat))]
		dane_l<-var_dflt_low$datamat[,-c(1:4,ncol(var_dflt_low$datamat))]

		#przeliczamy reszty
 		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-as.matrix(dane_h)%*%as.matrix(head(var1$varresult[[i]]$coefficients,-1))
    			var1$varresult[[i]]$residuals<-var_dflt_high$datamat[,i]-var1$varresult[[i]]$fitted.values

    			var2$varresult[[i]]$fitted.values<-as.matrix(dane_l)%*%as.matrix(head(var2$varresult[[i]]$coefficients,-1))
    			var2$varresult[[i]]$residuals<-var_dflt_low$datamat[,i]-var2$varresult[[i]]$fitted.values
  		}

		#Obliczamy stale
		c1<-residuals(var1)
		const1<-apply(c1,2,mean)
		#const1

		c2<-residuals(var2)
		const2<-apply(c2,2,mean)
		#const2


		### podmieniamy sta�e
		for (i in 1:4){
  			var1$varresult[[i]]$coefficients[4*var_dflt_high$p+1]<-const1[i]
  			var2$varresult[[i]]$coefficients[4*var_dflt_low$p+1]<-const2[i]
		}

  		#przeliczamy reszty uwzgledniajac stala
  		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-var1$varresult[[i]]$fitted.values+tail(var1$varresult[[i]]$coefficients,1)
    			var1$varresult[[i]]$residuals<-var1$varresult[[i]]$residuals-tail(var1$varresult[[i]]$coefficients,1)
    			
			var2$varresult[[i]]$fitted.values<-var2$varresult[[i]]$fitted.values+tail(var2$varresult[[i]]$coefficients,1)
    			var2$varresult[[i]]$residuals<-var2$varresult[[i]]$residuals-tail(var2$varresult[[i]]$coefficients,1)
  		}

		checkpoint.reszta<-list(apply(residuals(var1),2,mean),apply(residuals(var2),2,mean))



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


#ortogonalny
Kmax<-24
IRF<-IRF.noGDP<-rep(NaN,Kmax)
impuls<-1

F<-matrix(unlist(Acoef(model_var_1)),nrow=4)
F<- rbind(F
          ,cbind(diag(model_var_1$p-1) %x% diag(4),matrix(0,nrow=4*(model_var_1$p-1),ncol=4))
  )
  
F.noGDP<-F
F.noGDP[2,seq(0,4*(model_var_1$p-1),by=4)+impuls]<-0

A<-model_svar_1$A
A.noGDP<-model_svar_1$A
A.noGDP[2,implus]<-0

temp <- diag(1,4)
temp<-solve(A)%*%model_svar_1$B
temp<- rbind(cbind(temp,matrix(0,nrow=4,ncol=4*(model_var_1$p-1)))
            ,matrix(0,nrow=4*(model_var_1$p-1),ncol=4*(model_var_1$p))
)

temp.noGDP<-temp

for(k in 1:Kmax){
  IRF[k] <- temp[2,3]
  IRF.noGDP[k] <- temp.noGDP[2,3]
  temp   <- F%*%temp
  temp.noGDP   <- F.noGDP%*%temp.noGDP
}

dev.new()

  plot.ts(irf(model_var_1,impulse="WIBOR",response="CPI",n.ahead=24)$irf[[1]],lwd=2)
  lines(IRF,col="blue")
  lines(IRF.noGDP,col="green")
  lines(Psi(model_svar_1$var,nstep=23)[2,3,],col="red",lty=4)

  dev.print(pdf,paste("./PDF/WIBORonCPI_GDPoff.pdf",sep=""))
dev.off()


if(FALSE){
#nieortogonalny
Kmax<-24
IRF    <- rep(NaN,Kmax)
F<-Acoef(model_var_1)[[1]]
temp <- diag(1,4)
#temp<-solve(model_svar_1$A)

for(k in 1:Kmax){
  IRF[k] <- temp[2,3]
  temp   <- F%*%temp
}

plot.ts(IRF,col="blue")
lines(Phi(model_svar_1$var,nstep=23)[2,3,],col="red",lty=4)
}

dane<-var_dflt_high$y

### IRFY ###
h<-24 	#ustawiamy horyzont zabawy
ci<-0.95 	#szeroko�� CI
r<-1000	#liczba iteracji do bootstrappu

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
				,cumulative=TRUE
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

	dev.print(pdf,paste("./PDF/",b,"_short.pdf",sep=""))
	dev.off()
}

### IRFY ###
h<-12*4 	#ustawiamy horyzont zabawy
ci<-0.95 	#szeroko�� CI
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
				,cumulative=TRUE
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

	dev.print(pdf,paste("./PDF/",b,"_long.pdf",sep=""))
	dev.off()
}


#################### 3D IRF #########################################
MEGAIRFCPI<-NULL
MEGAIRFGDP<-NULL
MEGAIRFREER<-NULL

h<-24


for(p in seq(0,1,by=0.01)){
	
	model_var_irf<-var2
	if(var1$p>var2$p){
		model_var_irf<-var1
	}

	for (i in 1:4){
		vc1<-c(head(var1$varresult[[i]]$coefficients,-1),rep(0,4*(max(var1$p,var2$p)-var1$p)))
		vc2<-c(head(var2$varresult[[i]]$coefficients,-1),rep(0,4*(max(var1$p,var2$p)-var2$p)))
		const<-p*tail(var1$varresult[[i]]$coefficients,1)+(1-p)*tail(var2$varresult[[i]]$coefficients,1)
  			
  		nazwy<-names(model_var_irf$varresult[[i]]$coefficients)
		model_var_irf$varresult[[i]]$coefficients<-c(p*vc1+(1-p)*vc2,const)
		names(model_var_irf$varresult[[i]]$coefficients)<-nazwy

    		#model_var_irf$varresult[[i]]$fitted.values<-as.matrix(model_var_irf$datamat)%*%as.matrix(head(model_var_irf$varresult[[i]]$coefficients,-1))+tail(model_var_irf$varresult[[i]]$coefficients,1)
    		
		model_var_irf$varresult[[i]]$fitted.values<-as.matrix(model_var_irf$datamat[,-c(1:4)])%*%as.matrix(model_var_irf$varresult[[i]]$coefficients)
    		model_var_irf$varresult[[i]]$residuals<-model_var_irf$datamat[,i]-model_var_irf$varresult[[i]]$fitted.values
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

