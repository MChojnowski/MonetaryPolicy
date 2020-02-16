 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/42874/OneDrive - Bain/PhD/Workspace/MonetaryPolicy/MonetaryPolicy")
library(devtools)
#install_github("Mchojnowski/SentR")
set.seed(2019)

# source
source("./bin/config.R")
source("./bin/functions.R")
source("./bin/holidays.R")
source("./bin/PrivateLib.R")
source("./bin/DataPrepUS.R")
source("./bin/gridSearch2.R")

#szbsvar
#?szbsvar()
#c(0.7,0.2,2,1000,1,0.1,100)

##### SENTIMENTS #####
sentymenty <- sent_US

dane_monet<-window(dane_monet,start=poczatek_proby_US)
sentymenty <- window(sentymenty,start=poczatek_proby_US)
mci <- (mci - mean(mci))/sd(mci)
mci <- window(mci,start=poczatek_proby_US)


plot(sentymenty)
plot(mci)

##### MODEL SEARCH #####
grid_sent <- gridSearch2(
    dane_monet # Data
    ,1 # P Higher
    ,1 # P Lower
		,sentymenty		# Sentymenty
		,0			# srednia sentymentu
		,1		# sd sentymentu
		,0.1 # Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
		,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
		,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
		,stable_param = 1
	)

##### MODEL SEARCH #####
grid_mci <- gridSearch2(
  dane_monet # Data
  ,1 # P Higher
  ,1 # P Lower
  ,mci		# Sentymenty
  ,0		# srednia sentymentu
  ,1		# sd sentymentu
  ,0.1 # Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
  ,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
  ,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
  ,stable_param = 1
)


plotGrid(grid_sent)	
plotGrid(grid_mci)	

prob_sent <- pnorm(sentymenty,grid_sent$optparam[2],grid_sent$optparam[1])
prob_mci <- pnorm(mci,grid_mci$optparam[2],grid_mci$optparam[1])

plot(prob_sent)
plot(prob_mci)

var_sent <- STSVAR(dane_monet,1,1,prob_sent)
var_mci <- STSVAR(dane_monet,1,1,as.vector(prob_mci))

roots(var_sent$varH)
roots(var_sent$varL)

AmatMP <- diag(4)
AmatMP[lower.tri(AmatMP)] <- NA
BmatMP <- diag(4)
diag(BmatMP) <- NA

model_D <- VAR(dane_monet,p=1)

plot(irf(SVAR(var_sent$varH,Amat=AmatMP,Bmat=BmatMP),impulse = "pce", response="mon_pol_RBNZ"))
plot(irf(SVAR(var_sent$varL,Amat=AmatMP,Bmat=BmatMP),impulse = "pce", response="mon_pol_RBNZ"))


##### PLOTS IRF #####
zmienne <- colnames(dane_monet)
# Extracted sentiments
for (i in 3){
#for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_sent,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  #title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/US/",zmienne[i],"Shock_SENT.pdf",sep=""))
  dev.off()
}

# Extracted sentiments
for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_mci,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  #title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/US/",zmienne[i],"US_Shock_MCI.pdf",sep=""))
  dev.off()
}

## IRF 3D for sentiments
IRF3D <- irf3D(var_sent,c("mon_pol_RBNZ","pce"),c("pi","pce","reer_bis"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"mon_pol_RBNZ","pce",48,100)
dev.print(pdf,paste("./PDF/US/IRF3D_SENT.pdf",sep=""))
dev.off()

## IRF 3D for mci
IRF3D <- irf3D(var_mci,c("mon_pol_RBNZ","pce"),c("pi","pce","reer_bis"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"mon_pol_RBNZ","pce",48,100)
dev.print(pdf,paste("./PDF/US/IRF3D_MCI.pdf",sep=""))
dev.off()

##### OTHER PLOTS #####
## Probability SENT
dev.new()
plot(prob_mci, ylab="Probability", main="Probability of higher regime - SENT")
dev.print(pdf,paste("./PDF/US/Probability_SENT.pdf",sep=""))
dev.off()

## Probability mci
dev.new()
plot(prob_mci, ylab="Probability", main="Probability of higher regime - MCI")
dev.print(pdf,paste("./PDF/US/Probability_MCI.pdf",sep=""))
dev.off()

## Plot MCI, plot sent

## extracted sentiments
dev.new()
plot(sentymenty, ylab="", main="Extracted sentiments")
dev.print(pdf,paste("./PDF/US/SENT.pdf",sep=""))
dev.off()

## Probability SENT
dev.new()
plot(mci, ylab="Probability", main="Probability of higher regime - MCI")
dev.print(pdf,paste("./PDF/US/MCI.pdf",sep=""))
dev.off()