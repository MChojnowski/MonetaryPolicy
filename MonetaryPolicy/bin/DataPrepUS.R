dane_m <- ts(read.csv("./Dane/dane_USA.csv",sep=";",dec=".")[,-1],start=c(1964,1),freq=12)
dane_m <- dane_m[,c("pi","pce","mon_pol_RBNZ","reer_bis")]
# dane_m <- dane_m[,c("pi","mon_pol_RBNZ","pce","reer_bis")]

sent_US <- ts(read.csv("./Dane/USSent.csv",sep=",",dec=".")[,3],end=c(2018,3),freq=12)
mci <- read_xls('./Dane/tbmics.xls', sheet = "Dane")

###################################
mci <- mci %>% 
          filter(Time >= as.Date("1978-01-01")) %>%
          as.data.frame(.)
          

mci <- ts(mci[,2],start=c(1978,1),freq=12)

dane_monet <- cbind(
  diff(dane_m[,1])
  ,diff(dane_m[,2])
  ,diff(dane_m[,3])
  ,diff(dane_m[,4])
)

colnames(dane_monet) <- colnames(dane_m)

dane_monet <- window(dane_monet, start=c(1978,1), end=c(2018,3))
sent_US <- window(sent_US, start=c(1978,1), end=c(2018,3))
mci <- window(mci, start=c(1978,1), end=c(2018,3))
