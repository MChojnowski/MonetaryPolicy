dane_m<-ts(read.csv("./Dane/dane_ctrl.csv",sep=",",dec=".")[,-1],start=c(1993,1),freq=12)

###################################

wynagrodzenie<-log(dane_m[,"wag"]/dane_m[,"cpi"])
#lm(wynagrodzenie[73:length(wynagrodzenie)]~c(73:length(wynagrodzenie)))
#lm(wynagrodzenie[1:72]~c(1:72))
#lin<-2.669539+0.001928*73
#lin2<-2.331842+0.004056*73
wynagrodzenie[1:72]<-wynagrodzenie[1:72]+0.182353

#lm(wynagrodzenie[230:240]~c(230:240))
#lm(wynagrodzenie[200:228]~c(200:228))
#lin<-2.874277+0.000936*230
#lin2<-2.626849+0.001981*230
wynagrodzenie[230:length(wynagrodzenie)]<-wynagrodzenie[230:length(wynagrodzenie)]+0.007078

wynagrodzenie<-wynagrodzenie-stl(wynagrodzenie,s.window="periodic")$time.series[,1]

ev1<-rep(0,length(wynagrodzenie))
ev1[c(73,192)]<-1
ev1[180]<-0.5
arima(wynagrodzenie,order=c(2,0,0),xreg=ev1)$coef[4]
wynagrodzenie<-wynagrodzenie - (ev1 * arima(wynagrodzenie,order=c(2,0,0),xreg=ev1)$coef[4])

zatrudnienie<-log(dane_m[,"emp"])

Ev<-ts(rep(0,length(zatrudnienie)),start=start(zatrudnienie),freq=12)
Ev1<-ts(rep(1,24),start=c(1994,1),freq=12)
Ev2<-ts(rep(1,24),start=c(1996,1),freq=12)
Ev3<-ts(rep(1,12),start=c(1998,1),freq=12)
Ev32<-ts(rep(1,12),start=c(1999,1),freq=12)
Ev4<-ts(rep(1,36),start=c(2000,1),freq=12)
Ev5<-ts(rep(1,24),start=c(2003,1),freq=12)
Ev6<-ts(rep(1,12),start=c(2005,1),freq=12)
Ev62<-ts(rep(1,24),start=c(2006,1),freq=12)
Ev7<-ts(rep(1,36),start=c(2008,1),freq=12)
Ev8<-ts(rep(1,60),start=c(2011,1),freq=12)
Ev9<-ts(rep(1,12),start=c(2016,1),freq=12)
Ev10<-ts(rep(1,12),start=c(2017,1),freq=12)

Ev<-cbind(Ev
          ,Ev1
          ,Ev2
          ,Ev3
          ,Ev32
          ,Ev4
          ,Ev5
          ,Ev6
          ,Ev62
          ,Ev7
          ,Ev8
          ,Ev9
          ,Ev10
          )[,-1]
Ev[is.na(Ev)]<-0

model.zat<-arima(zatrudnienie,order=c(2,0,1),xreg=Ev)
zatrudnienie<-zatrudnienie-ts(Ev%*%model.zat$coef[-c(1:4)],start=start(zatrudnienie),freq=12)

wspol_kapital<-log(dane_m[,"car"])
Ev<-ts(rep(0,length(wspol_kapital)),start=start(wspol_kapital),freq=12)
a<-0.61
Ev22<-ts(c(a^0,a^3,a^5),start=c(1997,1),freq=12)
Ev11<-ts(c(a^5,a^3,a^1,a^0,a^1,a^3),start=c(1998,10),freq=12)

Ev1<-ts(apply(cbind(Ev,Ev11),1,sum,na.rm=TRUE),start=start(Ev),freq=12)
Ev2<-ts(apply(cbind(Ev,Ev22),1,sum,na.rm=TRUE),start=start(Ev),freq=12)

model.kap<-arima(wspol_kapital,order=c(2,0,3),xreg=data.frame(Ev1,Ev2),method="CSS")

wspol_kapital<-wspol_kapital-Ev1*model.kap$coef["Ev1"]-Ev2*model.kap$coef["Ev2"]

konsumpcja<-log(dane_m[,"rt"])

#lm(konsumpcja[53:73]~c(53:73))
#lm(konsumpcja[74:94]~c(74:94))
#lin<-4.45770+0.00847*73
#lin2<-4.991324+0.002609*73
konsumpcja[1:73]<-konsumpcja[1:73]+0.105


a<-0.3
ev<-ts(c(a^c(1,0,1)),start=c(2004,3),freq=12)
ev<-cbind(
	ts(0,start=start(konsumpcja),end=end(konsumpcja),freq=12),
	ev
	)
ev<-apply(ev,1,sum,na.rm=TRUE) %>%
	ts(.,start(konsumpcja),end=end(konsumpcja),freq=12)


konsumpcja<-konsumpcja - (ev * arima(konsumpcja,order=c(2,0,0),xreg=ev)$coef[4])



produkcja<-log(dane_m[,"pi"])
Ev<-ts(rep(0,length(produkcja)),start=start(produkcja),freq=12)
Ev1<-ts(c(1,1),start=c(2002,4),freq=12)
Ev<-ts(apply(cbind(Ev,Ev1),1,sum,na.rm=TRUE),start=start(produkcja),freq=12)
model.prod<-arima(produkcja,order=c(0,0,5),xreg=Ev)

produkcja<-produkcja-model.prod$coef["Ev"]*Ev

kredyty<-log(dane_m[,"loa"]/dane_m[,"cpi"])


ev1<-rep(0,length(kredyty))
ev1[c(90,142)]<-1
ev1[149]<-0.5
arima(kredyty,order=c(2,0,0),xreg=ev1)$coef[4]
kredyty<-kredyty - (ev1 * arima(kredyty,order=c(2,0,0),xreg=ev1)$coef[4])


pkb<-produkcja
cpi<-log(dane_m[,"cpi"])
wibor<-log(dane_m[,"wibor_1m"])
reer<-log(dane_m[,"reer"])
commodity<-log(dane_m[,"all_com_pri_ind"])
m3<-log(dane_m[,"m3"])
ea<-log(dane_m[,"pi_ea"])

dane_var<-na.omit(cbind(
		diff(wynagrodzenie),
		diff(zatrudnienie),
		diff(wspol_kapital),
		diff(konsumpcja),
		diff(produkcja),
		diff(kredyty)
		))

colnames(dane_var)<-c("Wynagrodzenie","Zatrudnienie","Wspol Kapitalowy","Konsumpcja","Produkcja","Kredyty")

sent<-get_sents(dane_var,u_oS,12,2,12)$sent
dane_var<-na.omit(cbind(dane_var,sent))


######################################################################


if(FALSE){
	start.dane<-c(1999,1)
	pkb<-window(pkb,start=start.dane)
	cpi<-window(cpi,start=start.dane)
	wibor<-window(wibor,start=start.dane)
	reer<-window(reer,start=start.dane)
}

#cpi.d<-diff(cpi)
ev1<-ts(rep(0,length(cpi)),start=start(cpi),freq=12)
ev1.1<-ts(1,start=c(1997,1),freq=12)
ev1<-ts(apply(cbind(ev1,ev1.1),1,sum,na.rm=TRUE),start=start(cpi),freq=12)

model.cpi<-arima(cpi,order=c(3,0,2),xreg=ev1,method="CSS")
cpi<-cpi-ev1*model.cpi$coef["ev1"]

###

ev1<-ts(rep(0,length(wibor)),start=start(wibor),freq=12)
ev1.1<-ts(c(0.5,1),start=c(1999,11),freq=12)
ev1<-ts(apply(cbind(ev1,ev1.1),1,sum,na.rm=TRUE),start=start(wibor),freq=12)

model.wibor<-arima(wibor,order=c(0,0,5),xreg=ev1)
wibor<-wibor-ev1*model.wibor$coef["ev1"]




dane_monet<-na.omit(
	cbind(
  		diff(pkb),
  		diff(cpi),
  		diff(wibor),
  		diff(reer)
	)
)


colnames(dane_monet)<-c("GDP","CPI","WIBOR","REER")

esi <- read.csv("./Dane/ESIPL.csv")[,2]
esi_norm <- ts((esi-mean(esi))/sd(esi),start=c(1993,10),freq=12)
