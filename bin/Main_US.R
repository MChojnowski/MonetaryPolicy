 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/42874/OneDrive - Bain/PhD/Workspace/MonetaryPolicy")
library(devtools)
install_github("Mchojnowski/SentR")
set.seed(2019)

# source
source("./bin/config.R")
source("./bin/functions.R")
source("./bin/holidays.R")
source("./bin/PrivateLib.R")
source("./bin/DataPrepUS.R")

#szbsvar
#?szbsvar()
#c(0.7,0.2,2,1000,1,0.1,100)

##### SENTIMENTS #####
sentymenty <- sent_US

#dane_monet<-window(dane_monet,start=poczatek_proby)
#sentymenty <- window(sentymenty,start=poczatek_proby)
#esi_sent <- window(esi_norm,start=poczatek_proby,end=end(dane_monet))

plot(sentymenty)

##### MODEL SEARCH #####
grid_sent <- gridSearch(
    dane_monet # Data
    ,1 # P Higher
    ,1 # P Lower
		,sentymenty		# Sentymenty
		,mean(sentymenty)			# srednia sentymentu
		,sd(sentymenty)			# sd sentymentu
		,0.1 # Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
		,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
		,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
	)

plotGrid(grid_sent)	

prob_sent <- pnorm(sentymenty,grid_sent$optparam[2],grid_sent$optparam[1])

var_sent <- STSVAR(dane_monet,1,1,prob_sent)

AmatMP <- diag(4)
AmatMP[lower.tri(AmatMP)] <- NA
BmatMP <- diag(4)
diag(BmatMP) <- NA

model_D <- VAR(dane_monet,p=1)

##### PLOTS IRF #####
zmienne <- colnames(dane_monet)
# Extracted sentiments
for (i in 3){
#for(i in 1:ncol(dane_monet)){
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