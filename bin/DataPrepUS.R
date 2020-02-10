dane_m<-ts(read.csv("./Dane/dane_USA.csv",sep=";",dec=".")[,-1],start=c(1964,1),freq=12)
dane_m <- dane_m[,c("pi","cpi","pce","reer_bis")]

sent_US <- ts(read.csv("./Dane/USSent.csv",sep=",",dec=".")[,3],end=c(2018,3),freq=12)

###################################


dane_monet <- diff(dane_m)
dane_monet[,4] <- tail(dane_m[,4],-1)

dane_monet <- window(dane_monet, start=c(1964,2), end=c(2018,3))
sent_US <- window(sent_US, start=c(1964,2), end=c(2018,3))
