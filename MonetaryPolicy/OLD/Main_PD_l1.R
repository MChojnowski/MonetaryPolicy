##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/piotr.dybka/Desktop/monetary policy/Monetary_policy")

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
#?szbsvar()


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

qwe1<-cumsum(qwe)/sd(cumsum(qwe))

######### for external value of sentiments
#qwe1<-read.csv("sentfull.csv")
#qwe1<-qwe1[,2]

plot(qwe1, type="l")

###### Exponential smoothing
qwe_sm <- holt(qwe1, alpha=0.1,  beta=FALSE, initial="simple", h=12) 
fit<-fitted(qwe_sm)
plot(fit)

###### liczymy ró¿nice
qwe2<-fit-qwe1
plot(qwe2, type="l")

sdp<-sd(qwe2)

###### liczymy prawdopodobieñstwa
prob1<-pnorm(qwe2,  mean = 0, sd = sdp)
plot(prob1, type="l")


#####################################################################################
###### Przekszta³cenie danych do modelu VAR

prob<-prob1
prob_prim<-1-prob

### DANE PREP ###
wynagrodzenie<-log(dane_m[,"wag"]/dane_m[,"cpi"])
zatrudnienie<-log(dane_m[,"emp"])
wspol_kapital<-log(dane_m[,"car"])
konsumpcja<-log(dane_m[,"rt"])
produkcja<-log(dane_m[,"pi"])
kredyty<-log(dane_m[,"loa"]/dane_m[,"cpi"])

pkb<-log(dane_m[,"pi"])
cpi<-log(dane_m[,"cpi"])
wibor<-log(dane_m[,"wibor_on"])
reer<-log(dane_m[,"reer"])

z<-length(wynagrodzenie)-length(prob)


pkb1<-prob*pkb[-c(1:z)]
cpi1<-prob*cpi[-c(1:z)]
wibor1<-prob*wibor[-c(1:z)]
reer1<-prob*reer[-c(1:z)]

pkb2<-prob_prim*pkb[-c(1:z)]
cpi2<-prob_prim*cpi[-c(1:z)]
wibor2<-prob_prim*wibor[-c(1:z)]
reer2<-prob_prim*reer[-c(1:z)]

dane_monet<-cbind(
  diff(pkb)+1,
  diff(cpi)+1,
  diff(wibor)+1,
  diff(reer)+1
)

dane_monet1<-cbind(
  diff(pkb1)+1,
  diff(cpi1)+1,
  diff(wibor1)+1,
  diff(reer1)+1
)

dane_monet2<-cbind(
  diff(pkb2)+1,
  diff(cpi2)+1,
  diff(wibor2)+1,
  diff(reer2)+1
)

colnames(dane_monet1)<-c("GDP","CPI","WIBOR","REER")
colnames(dane_monet2)<-c("GDP","CPI","WIBOR","REER")
colnames(dane_monet)<-c("GDP","CPI","WIBOR","REER")

#do sprawdzenia daty (ju¿ nie pamiêtam)
#dane_monet1<-window(dane_monet,start=c(2006,1))
#dane_monet2<-window(dane_monet,start=c(2006,1))

###### Model VAR 

### SVAR MATRICES ###
u_mon_xS<-matrix(0,4,4)
u_mon_xS[lower.tri(u_mon_xS)]<-NA
u_mon_xS[1,1]<-NA
u_mon_xS[2,2]<-NA
u_mon_xS[3,3]<-NA
u_mon_xS[4,4]<-NA
u_mon_xS

dane_monet1<-ts(dane_monet1,start=c(1993,2),freq=12)
dane_monet2<-ts(dane_monet2,start=c(1993,2),freq=12)

dane_all<-cbind(dane_monet,lag(dane_monet1,-1),lag(dane_monet2,-1))

r1<-lm(dane_all[,1]~dane_all[,-c(1:4)])
r2<-lm(dane_all[,2]~dane_all[,-c(1:4)])
r3<-lm(dane_all[,3]~dane_all[,-c(1:4)])
r4<-lm(dane_all[,4]~dane_all[,-c(1:4)])

comat1<-rbind(r1$coef[2:5],r2$coef[2:5],r3$coef[2:5],r4$coef[2:5])
comat2<-rbind(r1$coef[6:9],r2$coef[6:9],r3$coef[6:9],r4$coef[6:9])

colnames(comat1)<-NULL
colnames(comat2)<-NULL

var_dftl<-VAR(dane_monet,p=1)

var1<-var2<-var_dftl

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
dane<-dane_all[,1:4]
dane<-lag(dane, -1)
dane<-dane[-1,]

#przeliczamy reszty
for (i in 1:4){
var1$varresult[[i]]$residuals<-dane_all[-2,i]-dane[,1]*var1$varresult[[i]]$coefficients[1]+dane[,2]*var1$varresult[[i]]$coefficients[2]+dane[,3]*var1$varresult[[i]]$coefficients[3]+dane[,4]*var1$varresult[[i]]$coefficients[4]
}

for (i in 1:4){
  var2$varresult[[i]]$residuals<-dane_all[-2,i]-dane[,1]*var2$varresult[[i]]$coefficients[1]+dane[,2]*var2$varresult[[i]]$coefficients[2]+dane[,3]*var2$varresult[[i]]$coefficients[3]+dane[,4]*var2$varresult[[i]]$coefficients[4]
}

#Obliczamy stale
c1<-residuals(var1)
const1<-c(mean(c1[,1]),mean(c1[,2]),mean(c1[,3]),mean(c1[,4]))
const1

c2<-residuals(var2)
const2<-c(mean(c2[,1]),mean(c2[,2]),mean(c2[,3]),mean(c2[,4]))
const2

### podmieniamy sta³e
for (i in 1:4){
var1$varresult[[i]]$coefficients[5]<-const1[i]
}

for (i in 1:4){
  var2$varresult[[i]]$coefficients[5]<-const2[i]
}

#var1
#var2

#przeliczamy reszty jeszcze raz, ale tym razem ze sta³ymi!
for (i in 1:4){
  var1$varresult[[i]]$residuals<-dane_all[-2,i]-dane[,1]*var1$varresult[[i]]$coefficients[1]+dane[,2]*var1$varresult[[i]]$coefficients[2]+dane[,3]*var1$varresult[[i]]$coefficients[3]+dane[,4]*var1$varresult[[i]]$coefficients[4]-var1$varresult[[i]]$coefficients[5]
}

for (i in 1:4){
  var2$varresult[[i]]$residuals<-dane_all[-2,i]-dane[,1]*var2$varresult[[i]]$coefficients[1]+dane[,2]*var2$varresult[[i]]$coefficients[2]+dane[,3]*var2$varresult[[i]]$coefficients[3]+dane[,4]*var2$varresult[[i]]$coefficients[4]-var2$varresult[[i]]$coefficients[5]
}


#to jest tylko sprawdzenie, ze sie nie pomylilismy
c1<-residuals(var1)
const1<-c(mean(c1[,1]),mean(c1[,2]),mean(c1[,3]),mean(c1[,4]))
const1

c2<-residuals(var2)
const2<-c(mean(c2[,1]),mean(c2[,2]),mean(c2[,3]),mean(c2[,4]))
const2

#Tworzymy nowe ladne modele
model_var_1<-var1
model_var_2<-var2


model_svar_1<-SVAR(model_var_1,Amat=u_mon_xS)
model_svar_2<-SVAR(model_var_2,Amat=u_mon_xS)



### IRFY ###
h<-24 #ustawiamy horyzont zabawy
ci<-0.9 #szerokoœæ CI
r<-100 #liczba iteracji do bootstrappu

##### Zmiany CPI
b<-c("CPI") #ktora zmienna badamy
a<-c("WIBOR","REER","GDP") #definiujemy które zmienne wp³ywaj¹

par(mfcol=c(3,2))
for (i in 1:3){
    var_plot=irf(model_svar_1, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
    hc<-h+1
    plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
    lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
    lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
    abline(a = NULL, h = 0)
}
for (i in 1:3){
  var_plot=irf(model_svar_2, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b ), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}


##### Zmiany WIBOR
b<-c("WIBOR") #ktora zmienna badamy
a<-c("CPI","REER","GDP") #definiujemy które zmienne wp³ywaj¹

par(mfcol=c(3,2))
for (i in 1:3){
  var_plot=irf(model_svar_1, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}
for (i in 1:3){
  var_plot=irf(model_svar_2, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}

##### Zmiany REER
b<-c("REER") #ktora zmienna badamy
a<-c("CPI","WIBOR","GDP") #definiujemy które zmienne wp³ywaj¹

par(mfcol=c(3,2))
for (i in 1:3){
  var_plot=irf(model_svar_1, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}
for (i in 1:3){
  var_plot=irf(model_svar_2, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}

##### Zmiany GDP
b<-c("GDP") #ktora zmienna badamy
a<-c("CPI","WIBOR","REER") #definiujemy które zmienne wp³ywaj¹

par(mfcol=c(3,2))
for (i in 1:3){
  var_plot=irf(model_svar_1, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}
for (i in 1:3){
  var_plot=irf(model_svar_2, impulse = a[i], response=b, n.ahead = h, ortho=TRUE, boot=TRUE, runs=r, ci=ci)
  hc<-h+1
  plot(x=c(1:hc), y=unlist(var_plot$Lower), type="l", lwd = 3, lty=2,col="red", ylab=paste("Response of ",b), xlab=a[i], ylim=range(c(unlist(var_plot$Lower),unlist(var_plot$Upper))) )
  lines(x=c(1:hc),y=unlist(var_plot$Upper),type="l",lwd = 3, lty=2,col="red")
  lines(x=c(1:hc),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(a = NULL, h = 0)
}

