 ##### HEAD #####
rm(list=ls(all=TRUE))
setwd("C:/Users/pdybka-las/Desktop/MonPol/MonetaryPolicy")
library(devtools)
install_github("Mchojnowski/SentR")
set.seed(2019)

# source
source("./bin/config.R")
source("./bin/functions.R")
source("./bin/holidays.R")
source("./bin/PrivateLib.R")
source("./bin/DataPrep.R")

#szbsvar
#?szbsvar()
#c(0.7,0.2,2,1000,1,0.1,100)

##### SENTIMENTS #####
# Oblicza sentymenyty
# Zmienna w if - TRUE - oblicza; FALSE - wczytuje z pliku

lSENTY<-list(NULL)
lIS<-list(NULL)

if (TRUE){
# 0.95,0.15,1,0.1584893,1,6.3095734,1584.8931925
# 10/20/2020 : c(0.7,0.2,1.6,1,1,10,10)
# 11/20/2020 : c(0.75,0.2,1,1,1,10,1000)
  
  
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
    ,c(0.9,0.2,1,1,1,1,100) # SZBSVAR
  )
  
  sentymenty_diff <- ts(sent_org$sent,end=end(dane_var),freq=12)
  
  sentymenty <- ts(cumsum(sentymenty_diff)/sd(cumsum(sentymenty_diff)),freq=12,end=end(dane_var))
  sentymenty <- (sentymenty-mean(sentymenty))/sd(sentymenty)
}

sentymenty_all <- sentymenty

dane_monet<-window(dane_monet,start=poczatek_proby)
sentymenty <- window(sentymenty_all,start=poczatek_proby)
#sentymenty <- (sentymenty-mean(sentymenty))/sd(sentymenty)
esi_sent <- window(esi_norm,start=poczatek_proby,end=end(dane_monet))

plot(sentymenty)

##### MODEL SEARCH #####
grid_sent <- gridSearch2(
    dane_monet # Data
    ,1 # P Higher
    ,1 # P Lower
		,sentymenty		# Sentymenty
		,0			# srednia sentymentu
		,1		# sd sentymentu
		,0.25 # Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
		,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
		,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
		,stable_param = 0.999
	)

grid_esi <- gridSearch2(
  dane_monet # Data
  ,1 # P Higer
  ,1 # P Lower
  ,esi_sent			# Sentymenty
  ,0			# srednia sentymentu
  ,1			# sd sentymentu
  ,0.7			# Minimalna liczba obserwacji w przedziale prawdopodobienstw podanych nizej (im mniejsza liczba, tym ekstrema bardziej dostepne)
  ,100			# Liczba obserwacji w 'gridzie' - liczba kalkulacji rosnie kwadratowo !!!
  ,c(0.05,0.95)		# Przedzial dla minimalnej liczbe obserwacji
)

plotGrid(grid_sent)	
plotGrid(grid_esi)		

prob_sent <- pnorm(sentymenty,grid_sent$optparam[2],grid_sent$optparam[1])
prob_esi <- pnorm(esi_sent,grid_esi$optparam[2],grid_esi$optparam[1])

var_sent <- STSVAR(dane_monet,1,1,prob_sent)
roots(var_sent$varH)
roots(var_sent$varL)

var_esi <- STSVAR(dane_monet,1,1,prob_esi)

AmatMP <- diag(4)
AmatMP[lower.tri(AmatMP)] <- NA
BmatMP <- diag(4)
diag(BmatMP) <- NA

model_D <- VAR(dane_monet,p=1)

plot(prob_sent)
plot(irf(SVAR(var_sent$varH,Amat=AmatMP,Bmat=BmatMP),impulse = "WIBOR", response="CPI"))
plot(irf(SVAR(var_sent$varL,Amat=AmatMP,Bmat=BmatMP),impulse = "WIBOR", response="CPI"))

##### PLOTS IRF #####
zmienne <- colnames(dane_monet)
# Extracted sentiments

for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_sent,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/PL/",zmienne[i],"Shock_SENT.pdf",sep=""))
  dev.off()
}

# Extracted sentiments
for(i in 1:ncol(dane_monet)){
  dev.new()
  par(mfrow=c(3,3),oma=c(0,0,2,0))
  plotIRF(var_esi,model_D,AmatMP,BmatMP,zmienne[i],zmienne[-i],24,0.95,1000,c("High","Low","Control"))
  title(paste("Response on ",zmienne[i]," shock - extracted sentiments"),outer=TRUE)
  dev.print(pdf,paste("./PDF/PL/",zmienne[i],"Shock_ESI.pdf",sep=""))
  dev.off()
}

## IRF 3D for sentiments
IRF3D <- irf3D(var_sent,c("WIBOR","CPI"),c("GDP","CPI","REER"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"WIBOR","CPI",48,100)
dev.print(pdf,paste("./PDF/PL/IRF3D_SENT.pdf",sep=""))
dev.off()

## IRF 3D for ESI
IRF3D <- irf3D(var_esi,c("WIBOR","CPI"),c("GDP","CPI","REER"),AmatMP,BmatMP,24,100)

dev.new()
plotIRF3D(IRF3D,"WIBOR","CPI",48,100)
dev.print(pdf,paste("./PDF/PL/IRF3D_ESI.pdf",sep=""))
dev.off()

##### OTHER PLOTS #####
## Probability SENT
dev.new()
plot(prob_esi, ylab="Probability", main="Probability of higher regime - SENT")
dev.print(pdf,paste("./PDF/PL/Probability_SENT.pdf",sep=""))
dev.off()

## Probability ESI
dev.new()
plot(prob_esi, ylab="Probability", main="Probability of higher regime - ESI")
dev.print(pdf,paste("./PDF/PL/Probability_ESI.pdf",sep=""))
dev.off()

## Plot ESI, plot sent

## extracted sentiments
dev.new()
plot(sentymenty, ylab="", main="Extracted sentiments")
dev.print(pdf,paste("./PDF/PL/SENT.pdf",sep=""))
dev.off()

## Probability SENT
dev.new()
plot(esi_sent, ylab="Probability", main="Probability of higher regime - ESI")
dev.print(pdf,paste("./PDF/PL/ESI.pdf",sep=""))
dev.off()