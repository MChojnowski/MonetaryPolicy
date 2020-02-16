##### LIBRARIES #####
#library(MSBVAR)
#source("./bin/MSBVARfunctions.R")

for (file in dir("../MSBVAR/R")){
  source(paste0("../MSBVAR/R/",file))
  
}

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
library(readxl)

lag <- stats::lag
irf <- vars::irf

##### PARAMETRY PL #####
poczatek_proby <- c(1999,1)
poczatek_model <- c(1999,1)

u_oS<-matrix(0,6,6)
u_oS[1,c(2)]<-NA
u_oS[2,c(1)]<-NA
u_oS[3,c(1,2)]<-NA
u_oS[4,c(1,5)]<-NA
u_oS[5,c(2,4)]<-NA
u_oS[6,c(3,4,5)]<-NA
diag(u_oS)<-1

##### PARAMETRY US #####
poczatek_proby_US <- c(1991,1)
poczatek_model_US <- c(1991,1)

u_oS_US<-matrix(0,6,6)
u_oS_US[1,c(2)]<-NA
u_oS_US[2,c(1)]<-NA
u_oS_US[3,c(1,2)]<-NA
u_oS_US[4,c(1,5)]<-NA
u_oS_US[5,c(2,4)]<-NA
u_oS_US[6,c(3,4,5)]<-NA
diag(u_oS_US)<-1