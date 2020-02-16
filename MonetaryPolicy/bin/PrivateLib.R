### RMSE
rmse<-function(x,type="rmse"){
	x<-na.omit(x)
	x<-x[x<Inf]
	x<-x[x>-Inf]
	if(type=="rmse"){
		return(sqrt(mean(x^2)))
	} else if(type=="mse"){
		return(mean(abs(x)))
	}
}

### SENTIMENT EXTRACTION
get_sents<-function(f_dane,Amat,opoznienie,szlif,frek=12){
	f_poczatek<-c(1,1,1)
	if(sum(class(dane_var)=="ts")>0){
		f_poczatek[1]<-start(f_dane)[1]	
		f_poczatek[2]<-start(f_dane)[2]	
		f_poczatek[3]<-start(f_dane)[3]	
		f_poczatek[is.na(f_poczatek)]<-1
		f_poczatek<-as.Date(paste(f_poczatek[1],f_poczatek[2],f_poczatek[3]),"%Y %m %d")
		f_poczatek<-f_poczatek+
				years(
					if(frek==1){opoznienie}else{0}
				)+
				months(
					if(frek==12){opoznienie}else{0}
				)+
				days(
					if(frek==365){opoznienie}else{0}
				)
		f_poczatek<-c(year(f_poczatek),month(f_poczatek),day(f_poczatek))
	}

	Bmat<-matrix(0,ncol(f_dane),ncol(f_dane))
	diag(Bmat)<-NA
	
	model_var<-VAR(f_dane,p=opoznienie)
	model_svar<-SVAR(model_var,max.iter=1000,Amat=Amat,Bmat=Bmat,estmethod="direct")

	t(residuals(model_var))->res

	signal<-(solve(model_svar$A)%*%model_svar$B%*%res)
	signal<-apply(signal,1,function(x){(x)/sd(x)})

	x2<-NULL
	for(fi in 1:(ncol(f_dane)/2)){
		x2<-cbind(x2,
			signal[,(ncol(f_dane)/2)+fi]-apply(as.matrix(signal[,c(1:(ncol(f_dane)/2))[-fi]]),1,sum)
		)
	}
	#x2<-apply(x2,2,function(x){(x-mean(x))/sd(x)})
	x2<-apply(x2,1,sum)
	x2<-c(
		forecast::forecast(arima(rev(x2),order=c(3,0,0)),h=3)$mean,
		x2,
		forecast::forecast(arima(x2,order=c(3,0,0)),h=3)$mean
	)
	return(list(sent=ts(tail(head(hpfilter(x2,type="lambda",freq=szlif)$trend,-3),-3),start=f_poczatek,freq=frek),media=ts(tail(head(hpfilter(x2,type="lambda",freq=szlif)$cycle,-3),-3),start=f_poczatek,freq=frek)))
}

### Bayesian Sim-Zhao
getSentsB<-function(f_dane,Amat,opoznienie,opoz_B,szlif,bvec=c(1,1,1,1,1,1,1),frek=12){
	f_poczatek<-c(1,1,1)
	if(sum(class(dane_var)=="ts")>0){
		f_poczatek[1]<-start(f_dane)[1]	
		f_poczatek[2]<-start(f_dane)[2]	
		f_poczatek[3]<-start(f_dane)[3]	
		f_poczatek[is.na(f_poczatek)]<-1
		f_poczatek<-as.Date(paste(f_poczatek[1],f_poczatek[2],f_poczatek[3]),"%Y %m %d")
		f_poczatek<-f_poczatek+
				years(
					if(frek==1){opoznienie}else{0}
				)+
				months(
					if(frek==12){opoznienie}else{0}
				)+
				days(
					if(frek==365){opoznienie}else{0}
				)
		f_poczatek<-c(year(f_poczatek),month(f_poczatek),day(f_poczatek))
	}

	Bmat<-matrix(0,ncol(f_dane),ncol(f_dane))
	diag(Bmat)<-NA
	
	model_var<-VAR(f_dane,p=opoznienie)
	model_svar<-SVAR(model_var,max.iter=1000,Amat=Amat,Bmat=Bmat,estmethod="direct")
	model_bsvar<-szbsvar(f_dane,opoz_B,,bvec[1],bvec[2],bvec[3],bvec[4],bvec[5],bvec[6],bvec[7],model_svar$A,12)
	res<-t(model_bsvar$structural.innovation)

	#t(residuals(model_var))->res
	
	signal<-(solve(model_bsvar$A0.mode)%*%res)	
	signal<-apply(signal,1,function(x){(x)/sd(x)})

	x2<-NULL
	for(fi in 1:(ncol(f_dane)/2)){
		x2<-cbind(x2,
			signal[,(ncol(f_dane)/2)+fi]-apply(as.matrix(signal[,c(1:(ncol(f_dane)/2))[-fi]]),1,sum)
		)
	}
	#x2<-apply(x2,2,function(x){(x-mean(x))/sd(x)})
	x2<-apply(x2,1,sum)
	x2<-c(
		forecast::forecast(arima(rev(x2),order=c(3,0,0),method="CSS"),h=3)$mean,
		x2,
		forecast::forecast(arima(x2,order=c(3,0,0),method="CSS"),h=3)$mean
	)
	return(list(
		sent=ts(tail(head(hpfilter(x2,type="lambda",freq=szlif)$trend,-3),-3),start=f_poczatek,freq=frek),
		media=ts(tail(head(hpfilter(x2,type="lambda",freq=szlif)$cycle,-3),-3),start=f_poczatek,freq=frek))
		)
}

### DATA ARRAY 
inter_daty<-function(pocz,kon){
	dni<-seq(as.Date(pocz,"%d-%m-%Y"),as.Date(kon,"%d-%m-%Y"),by=1)
	dni<-dni[weekdays(dni)!="Sunday"]
	dni<-dni[!(dni %in% holidays)]
	dni<-cbind(dni,1)
	dni[weekdays(as.Date(dni[,1]))=="Saturday",2]<-sob
	dni<-aggregate(dni[,2],by=list(format(as.Date(dni[,1]),"%Y-%m")),sum)
	return(dni)
}

### DATES
daty<-function(pocz,n){

	kon<-pocz
	month(kon)<-month(kon)+nrow(dane_m)
	day(kon)<-day(kon)-1
	dni<-inter_daty(pocz,kon)
	konn<-kon+1
	month(konn)<-month(konn)+15
	day(konn)<-day(konn)-1
	dni_fory<-inter_daty(kon+1,konn)
	dni_tytul<-paste(format(seq(kon+1,konn,by=31),"%b'%y"))

	return(list(kon=kon,konn=konn,dni=dni,fory=dni_fory,tytul=dni_tytul))
}

### DATA REDIFINITION
zrob_zmienna<-function(input,czy.log=FALSE){
	x<-ts(c(1,cumprod(input/100)),start=c(2006,12),freq=12)
	if (czy.log) {
		x<-ts(log(c(1,cumprod(input/100))),start=c(2006,12),freq=12)
	}
	xd<-stl(x,s.window="periodic")
	x_xS<-xd$time.series[,2]+xd$time.series[,3]
	x_tr<-xd$time.series[,2]
	return(list(org=x,no.seas=x_xS,seas=x_tr))
}

### PROGNOZA
### InformationSet
### Dopasuj.Google
### Prognoza.Linear
### Prognoza.NonLinear
### Prognoza.ANN

#Load.InformationSet<-function(){
#}


Save.InformationSet<-function(f_IS,f_FileName){
	write.table(f_IS,f_FileName)		
}


#Dopasuj.Google<-function(){
#}


InformationSet<-function(f_sent,f_dane,f_mc=100,f_sp=1,f_st=10e-04,f_rep=1,opt.start=NULL){	
	
	#Wykryj bledy parametrow

	#Revise
	zzz<-f_sent
	omega_t<-""
	if(!is.null(opt.start)){	
		omega_t<-opt.start
		f_model<-spikeslab(f_sent~f_dane[,omega_t],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)	
		zzz<-f_sent-as.matrix(f_dane[,omega_t])%*%as.matrix(f_model$bma.scale)
	}
	cat("Revised\n")

	#Learn
	check_model<-Inf
	for(model_repet in 1:f_rep){
		tf_model<-spikeslab(zzz~f_dane[,!(colnames(f_dane) %in% omega_t)],n.iter1=5000,n.inter2=1000,ntree=750,bigp.smalln=TRUE)
		if(tf_model$sigma.hat*log(tf_model$phat)<check_model){
			f_model<-tf_model
			check_model<-tf_model$sigma.hat*log(tf_model$phat)
		}
	}	
	cat("Learnt\n")
	
	omega_t1<-rownames(f_model$summary)[abs(f_model$summary$bma)>f_st]
	omega_t1<-substr(omega_t1,regexpr("]",omega_t1)[1]+1,nchar(omega_t1))

	# Forget
	f_model<-spikeslab(f_sent~f_dane[,colnames(f_dane) %in% c(omega_t,omega_t1)],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)
	omega_t1<-rownames(f_model$summary)[abs(f_model$summary$bma)>f_st]
	omega_t1<-substr(omega_t1,regexpr("]",omega_t1)[1]+1,nchar(omega_t1))[1:min(length(omega_t1),f_mc)]
	f_model<-spikeslab(f_sent~f_dane[,omega_t1],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)

	check_model<-Inf
	for(model_repet in 1:f_rep){
		tf_model<-spikeslab(f_sent~f_dane,n.iter1=5000,n.inter2=1000,ntree=750,intercept=TRUE,bigp.smalln=TRUE)
		if(tf_model$sigma.hat*log(tf_model$phat)<check_model){
			f_model_exo<-tf_model
			check_model<-tf_model$sigma.hat*log(tf_model$phat)
		}
	}
	cat("Forgot\n")

	#Adapt
	omega_t1_exo<-rownames(f_model_exo$summary)[abs(f_model_exo$summary$bma)>f_st]
	omega_t1_exo<-gsub("f_dane","",omega_t1_exo)
	omega_t1_exo<-omega_t1_exo[1:min(length(omega_t1_exo),f_mc)]

	if (identical(character(0),omega_t1_exo)){
		f_model_exo<-f_model
		omega_t1_exo<-omega_t1
	} else {
		f_model_exo<-spikeslab(f_sent~f_dane[,omega_t1_exo],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)
	}
	cat("Adapted\n")		

	if(f_sp*f_model_exo$sigma.hat<f_model$sigma.hat){
		omega_t1<-omega_t1_exo
		f_model<-f_model_exo
	}
	rm(f_model_exo)

	#Info print
	persy_temp<-NA
	if(!is.null(opt.start)){
		persy_temp<-sum(table(c(omega_t1,omega_t))==2)/length(unique(c(omega_t,omega_t1)))
	}
	return(list(
		IS=omega_t1,
		persistency=persy_temp,
		model=f_model
	))
}

Prognoza.Linear<-function(f_x,f_xreg=NULL,f_exoreg=NULL,pq=c(1,1,1,1),fhor=12){

	ARIMA<-NULL
	ordery<-NULL

	assign("f_xreg",f_xreg,envir=parent.frame())
	assign("f_exoreg",f_exoreg,envir=parent.frame())
	assign("fhor",fhor,envir=parent.frame())

	for(fp in pq[1]:pq[3]){
		for(fq in pq[2]:pq[4]){
			ordery<-c(ordery,paste("(",fp,"0",fq,")",sep=""))
			ARIMA<-cbind(ARIMA,forecast(arima(f_x,order=c(fp,0,fq),method="CSS",xreg=f_xreg),xreg=f_exoreg,h=fhor)$mean)
		}
	}
	
	if (length(ordery)>1){colnames(ARIMA)<-ordery}
	rm(f_xreg,f_exoreg,fhor,envir=parent.frame())
	return(ARIMA)
}

Prognoza.NonLinear<-function(f_x,f_xreg=NULL,f_exoreg=NULL,pq=c(1,1,1,1),fhor=12,fm=max(pq[3],pq[4])+1){
	LSTAR<-NULL
	ordery<-NULL

	assign("f_x",f_x,envir=parent.frame())
	assign("f_xreg",f_xreg,envir=parent.frame())
	assign("f_exoreg",f_exoreg,envir=parent.frame())

	controls <- list(trim=0.3)

	if(is.null(f_xreg)){
		for(fp in pq[1]:pq[3]){
			for(fq in pq[2]:pq[4]){
				ordery<-c(ordery,paste("(",fp,"0",fq,")",sep=""))
				LSTAR<-cbind(LSTAR,
						tryCatch({
							predict(lstar(f_x,m=max(fp,fq),mL=fp,mH=fq,thDelay=1,trace=FALSE,include="const",control=list(maxit=4500)),n.ahead=fhor)
							},error=function(a){return(forecast(arima(f_x,order=c(fp,0,0)),h=12)$mean)})
					)
			}
		}
	} else {
		for(fp in pq[1]:pq[3]){
			for(fq in pq[2]:pq[4]){
				ordery<-c(ordery,paste("(",fp,"0",fq,")",sep=""))
				LSTAR<-cbind(LSTAR,
						tryCatch({	
							predict(lstar(f_x,m=max(fp,fq),mL=fp,mH=fq,thVar=f_xreg,trace=FALSE,include="const",th=0,control=list(maxit=4500), starting.control=controls),thVar=f_exoreg,n.ahead=fhor)
							},error=function(a){return(rep(NA,12))})
					)
			}
		}

	}
#forecast(arima(f_x,order=c(fp,0,0)),h=12)$mean)})
	if (length(ordery)>1){colnames(LSTAR)<-ordery}
	rm(f_x,f_xreg,f_exoreg,envir=parent.frame())
	return(LSTAR)
}

Prognoza.ANNARMA<-function(f_x,f_sent,f_gt,f_for.gt,f_hid=1,f_n=16,f_which=NULL,f_prewagi=NULL,saveWts=FALSE,okres=NULL){
	RESULT<-NULL
	ordery<-NULL
	f_for.sent<-rep(0,12)
	assign("f_sent",f_sent,envir=parent.frame())

	final.wagi<-NULL

	val_per<-1:12
	lear_per<-13:57
	
	fsent.diffex<-f_sent-ets(f_sent,alpha=0.4)$fitted

	for(k in 1:9){	
		ordery<-c(ordery,paste("(",k %% 3,"0",floor(k/3),")",sep=""))		
		modar<-arima(f_x,order=c(k %% 3,0,floor(k/3)),xreg=fsent.diffex,method="CSS")
		res_nnet<-modar$res+modar$coef["fsent.diffex"]*fsent.diffex
	
		restarget<-cbind(
			as.numeric(res_nnet<qnorm(1/6,0,sd(res_nnet))),
			as.numeric(qnorm(1/6,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(2/6,0,sd(res_nnet))),
			as.numeric(qnorm(2/6,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(3/6,0,sd(res_nnet))),
			as.numeric(qnorm(3/6,0,sd(res_nnet))<res_nnet & res_nnet<qnorm(4/6,0,sd(res_nnet))),
			as.numeric(qnorm(4/6,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(5/6,0,sd(res_nnet))),
			as.numeric(qnorm(5/6,0,sd(res_nnet))<=res_nnet)
		)
		
		restarget<-restarget[,apply(restarget[lear_per,],2,sum)>0]

	
		wagi<-foreach(o=1:f_n,.combine=rbind)%dopar%{
			nnet(restarget[lear_per,]~.,size=f_hid,data=f_gt[lear_per,],trace=FALSE,softmax=TRUE)$wts
		}

		if(!is.null(f_prewagi)){
		  wagi<-tryCatch({rbind(wagi,f_prewagi[k,])},error=function(e){wagi})
		}
		
		val<-foreach(o=1:nrow(wagi),.combine=rbind)%dopar%{
			rmse(
				predict(
					nnet(restarget[lear_per,]~.,size=f_hid,data=f_gt[lear_per,],trace=FALSE,softmax=TRUE,Wts=wagi[o,])
				,f_gt[val_per,]) %*% apply((restarget[lear_per,]*c(res_nnet[lear_per])),2,function(x){mean(x[x!=0])})
				-res_nnet[val_per]
			)
		}
			
		siec<-nnet(restarget[val_per,]~.,size=f_hid,data=f_gt[val_per,],trace=FALSE,softmax=TRUE,Wts=wagi[which.min(val),],maxit=250)
		
		final.wagi<-cbind(final.wagi,ts(wagi[which.min(val),]))
		
		RESULT<-cbind(RESULT,
				forecast(modar,xreg=data.frame(fsent.diffex=f_for.sent),h=12)$mean+
				ts((predict(siec,f_for.gt)%*%apply((restarget*c(res_nnet)),2,function(x){mean(x[x!=0])})),freq=12,start=start(f_for.gt))
			)
	}

	if(saveWts){write.csv(final.wagi,paste("./ANN Weights/ARMAANN_",f_which,"_",okres,".csv",sep=""))}


	if (length(ordery)>1){colnames(RESULT)<-ordery}
	return(RESULT)
}


Prognoza.ANNLSTAR<-function(f_x,f_sent,f_gt,f_for.gt,pq=c(1,1,1,1),f_hid=1,f_n=16,f_which=NULL,f_prewagi=NULL,saveWts=FALSE,okres=NULL){
	RESULT<-NULL
	ordery<-NULL
	fm<-max(pq[3],pq[4])+1
	
	final.wagi<-NULL

	

#	if(!is.null(f_w_pre)){
#   prewagi<-tryCatch({read.csv(paste("./ANN Weights/LSTARANN_",okres,".csv",sep=""))[f_which,-1]},error=function(e){return(NULL)})
#		presiec<-tryCatch({nnet(f_sent~.,size=f_hid,data=f_gt,trace=FALSE,linout=TRUE,Wts=f_w_pre)},error=function(e){return(NULL)})
#	}
	
	val_per<-1:12
	lear_per<-13:57
			
	wagi<-foreach(o=1:f_n,.combine=rbind)%dopar%{
			nnet(f_sent[lear_per]~.,size=f_hid,data=f_gt[lear_per,],trace=FALSE,linout=TRUE)$wts
	}
    
	if(!is.null(f_prewagi)){
		  wagi<-tryCatch({rbind(wagi,f_prewagi)},error=function(e){wagi})
	}
		
	val<-foreach(o=1:nrow(wagi),.combine=rbind)%dopar%{
  		rmse(
				predict(nnet(f_sent[lear_per]~.,size=f_hid,data=f_gt[lear_per,],trace=FALSE,linout=TRUE,Wts=wagi[o,]),f_gt[val_per,])-
				f_sent[val_per]
			)
	}
			
	siec<-nnet(f_sent[lear_per]~.,size=f_hid,data=f_gt[lear_per,],trace=FALSE,linout=TRUE,Wts=wagi[which.min(val),],maxit=250)
	final.wagi<-rbind(final.wagi,wagi[which.min(val),])

	#new_siec<-tail(cumsum(predict(siec)),1)+cumsum(predict(siec,f_for.gt))

	fsent.diffex<-ts(cumsum(c(f_sent,predict(siec,f_for.gt)))-ets(cumsum(c(f_sent,predict(siec,f_for.gt))),alpha=0.4)$fitted,freq=12,start=start(f_sent))

	for(fp in pq[1]:pq[3]){
	  for(fq in pq[2]:pq[4]){
	    ordery<-c(ordery,paste("(",fp,"0",fq,")",sep=""))
	    
    		RESULT<-cbind(RESULT,
		  		ts(tryCatch({	
			  		predict(lstar(f_x,m=max(fp,fq),mL=fp,mH=fq,thVar=head(fsent.diffex,-12),trace=FALSE,include="const",th=0,control=list(maxit=9000),starting.control=list(trim=0.3)),thVar=tail(fsent.diffex,12),n.ahead=12)
					},error=function(a){return(rep(NA,12))})
				,freq=12,start=start(forecast(arima(f_x,order=c(1,0,0)),h=3)$mean))
		  	)
    }
	}

	if(saveWts){write.csv(final.wagi,paste("./ANN Weights/LSTARANN_",f_which,"_",okres,".csv",sep=""))}


	if (length(ordery)>1){colnames(RESULT)<-ordery}
	return(RESULT)
}


##########################################

mc_prognoza2<-function(f_dane_var_sent,f_dane_gt,okres,dl=60,f_mc=100,f_sp=1,szukaj=TRUE){
	
	source("./ord_for_gg.R")

	hor<-c(1,3,6,12)
	OUT_PROGNOZA<-NULL
	KWERY<-NULL
	PERSY<-NULL
	MEGASENT<-matrix(NA,length(3:57),length(3:57)+dl)	

	kwerendy<-list(NULL)
	length(kwerendy)<-length(okres)
	
	u_oS<-matrix(0,4,4)

	u_oS[1,c(2)]<-NA
	u_oS[2,c(1)]<-NA
	u_oS[3,c(1,4)]<-NA
	u_oS[4,c(2,3)]<-NA
	diag(u_oS)<-1
	
	B_oS<-matrix(0,4,4)
	diag(B_oS)<-NA

	for (i in okres){
		cat("\n")
		cat(paste("START:",i,"\n"))
		cat(paste("PROGRESS: ",round(100*which(okres==i)/length(okres),1),"%","\n",sep=""))
	
		f_dane_var<-ts(f_dane_var_sent[(i-1):(i+dl-1),],start=c(1,1),freq=12)
		f_rynek<-ts(f_dane_var[-(1:p_lag),2],start=c(1,1),freq=12)
		f_wwuk<-ts(f_dane_var[-(1:p_lag),3],start=c(1,1),freq=12)
		f_sent<-ts(c(get_sents(f_dane_var[,-c(3,6,7)])$sent),start=c(1,1),freq=12)
		f_sent_gt<-ts(f_dane_gt[(i-1+p_lag):(i+dl-1),],start=c(1,1),freq=12)
		f_temp_wwuk<-hpfilter(f_wwuk,freq=1)$trend
		hype<-stl(ts(f_sent,freq=9),s.window="periodic")$time.series[,1]
		f_sent_oS<-ts(apply(stl(ts(f_sent,freq=9),s.window="periodic")$time.series[,-1],1,sum),start=c(1,1),freq=12)
		hype_wwuk<-stl(ts(diff(c(f_temp_wwuk)),freq=8),s.window="periodic")$time.series[,1]
		MEGASENT[i-2,(i-2):(i-3+length(f_sent))]<-f_sent
			
		assign("f_sent",f_sent,envir = .GlobalEnv)
		assign("f_wwuk",f_wwuk,envir = .GlobalEnv)
		assign("hype",hype,envir = .GlobalEnv)
		assign("hype_wwuk",hype_wwuk,envir = .GlobalEnv)

		f_gt_mod<-f_sent_gt

		#for (i_mod in 1:4){
		#	model_temp<-spikeslab(f_sent~f_sent_gt,n.iter1=100,n.inter2=50,ntree=300,intercept=TRUE,bigp.smalln = TRUE,max.var=24,bigp.smalln.factor=0.6)
		#	if(f_model$sigma.hat>model_temp$sigma.hat){
		#		f_model<-model_temp
		#	}
		#}
		#rm(model_temp,i_mod)

		ARIMAS<-NULL
		ARSENTIS<-NULL
		ARWUKIS<-NULL
		LSTARY<-NULL
		LSTARYWWUK<-NULL
		GUGLE<-NULL
		ARNEUNET<-NULL
		NEUNET<-NULL
		NEUNETGT<-NULL
		NEUNETGTMEAN<-NULL
		VARY<-NULL

		#c(ceiling((hor-1)/4)+1,0,0)

		fory<-NULL
		for(fi in 1:(ncol(f_dane_gt)/7)){
			fory<-cbind(fory,forecast(arima(f_dane_gt[,fi],order=ord_for_gg[[fi]],method="CSS"),h=12)$mean)
		}

		fory<-cbind(
			fory,
			rbind(
				f_dane_gt[nrow(f_dane_gt),1:(ncol(f_dane_gt)/7)],fory[1:11,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-1):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],fory[1:10,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-2):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],fory[1:9,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-3):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],fory[1:8,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-4):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],fory[1:7,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-5):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],fory[1:6,]
			)
		)
		
		#201
		foryAR<-NULL
		for(fi in 1:(ncol(f_dane_gt)/7)){
			foryAR<-cbind(foryAR,forecast(arima(f_dane_gt[,fi],order=ord_for_gg[[fi]],method="CSS"),h=12)$mean)
		}
		foryAR<-cbind(
			foryAR,
			rbind(
				f_dane_gt[nrow(f_dane_gt),1:(ncol(f_dane_gt)/7)],foryAR[1:11,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-1):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryAR[1:10,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-2):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryAR[1:9,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-3):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryAR[1:8,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-4):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryAR[1:7,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-5):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryAR[1:6,]
			)
		)

		#mean
		foryMEAN<-NULL
		for(fi in 1:(ncol(f_dane_gt)/7)){
			foryMEAN<-cbind(foryMEAN,rep(mean(f_dane_gt[,fi]),12))
		}
		foryMEAN<-cbind(
			foryMEAN,
			rbind(
				f_dane_gt[nrow(f_dane_gt),1:(ncol(f_dane_gt)/7)],foryMEAN[1:11,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-1):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryMEAN[1:10,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-2):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryMEAN[1:9,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-3):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryMEAN[1:8,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-4):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryMEAN[1:7,]
			),
			rbind(
				f_dane_gt[(nrow(f_dane_gt)-5):(nrow(f_dane_gt)),1:(ncol(f_dane_gt)/7)],foryMEAN[1:6,]
			)
		)

		rm(fi)
		colnames(foryMEAN)<-colnames(foryAR)<-colnames(fory)<-colnames(f_dane_gt)
		omega_t1<-try(read.table(paste("./TempTXT/Kwerendy_",i,".txt",sep=""))[,1])
		persy_temp<--1

	if (szukaj){
		if (i==3){
			wut<-Inf
			for(qqq in 1:16){
			tf_model<-spikeslab(f_sent~f_gt_mod,n.iter1=5000,n.inter2=1000,ntree=750,intercept=TRUE,bigp.smalln=TRUE)
			if(tf_model$sigma.hat*log(tf_model$phat)<wut){
				f_model<-tf_model
				wut<-tf_model$sigma.hat*log(tf_model$phat)
			}}			
			omega_t1<-gsub("f_gt_mod","",rownames(f_model$summary)[round(f_model$summary$bma,4)!=0])
			omega_t1<-omega_t1[1:min(length(omega_t1),f_mc)]
			
		} else {
			stare_kwery<-read.table(paste("./TempTXT/Kwerendy_",i-1,".txt",sep=""))
			omega_t<-as.character(stare_kwery[,1])
			f_model<-spikeslab(f_sent~f_gt_mod[,omega_t],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)

			zzz<-f_sent-as.matrix(f_gt_mod[,omega_t])%*%as.matrix(f_model$bma.scale)

			wut<-Inf
			for(qqq in 1:1){
			tf_model<-spikeslab(zzz~f_gt_mod[,!(colnames(f_gt_mod) %in% omega_t)],n.iter1=5000,n.inter2=1000,ntree=750,bigp.smalln=TRUE)
			if(tf_model$sigma.hat*log(tf_model$phat)<wut){
				f_model<-tf_model
				wut<-tf_model$sigma.hat*log(tf_model$phat)
			}}	

			#f_model<-spikeslab(zzz~f_gt_mod[,!(colnames(f_gt_mod) %in% omega_t)],n.iter1=5000,n.inter2=1000,ntree=750,bigp.smalln=TRUE)
	
			omega_t1<-rownames(f_model$summary)[round(f_model$summary$bma,4)!=0]
			omega_t1<-substr(omega_t1,47,nchar(omega_t1))

			f_model<-spikeslab(f_sent~f_gt_mod[,c(omega_t,omega_t1)],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)
			omega_t1<-rownames(f_model$summary)[round(f_model$summary$bma,4)!=0]
			omega_t1<-substr(omega_t1,33,nchar(omega_t1))[1:min(length(omega_t1),f_mc)]
			f_model<-spikeslab(f_sent~f_gt_mod[,omega_t1],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)

			wut<-Inf
			for(qqq in 1:1){
			tf_model<-spikeslab(f_sent~f_gt_mod,n.iter1=5000,n.inter2=1000,ntree=750,intercept=TRUE,bigp.smalln=TRUE)
			if(tf_model$sigma.hat*log(tf_model$phat)<wut){
				f_model_exo<-tf_model
				wut<-tf_model$sigma.hat*log(tf_model$phat)
			}}

			#f_model_exo<-spikeslab(f_sent~f_gt_mod,n.iter1=5000,n.inter2=1000,ntree=750,intercept=TRUE,bigp.smalln=TRUE)
			omega_t1_exo<-rownames(f_model_exo$summary)[round(f_model_exo$summary$bma,4)!=0]
			omega_t1_exo<-substr(omega_t1_exo,9,nchar(omega_t1_exo))
			omega_t1_exo<-omega_t1_exo[1:min(length(omega_t1_exo),f_mc)]

			if (identical(character(0),omega_t1_exo)){
				f_model_exo<-f_model
				omega_t1_exo<-omega_t1
			} else {
				f_model_exo<-spikeslab(f_sent~f_gt_mod[,omega_t1_exo],n.iter1=500,n.inter2=100,ntree=75,intercept=TRUE)
			}
			

			if(f_sp*f_model_exo$sigma.hat<f_model$sigma.hat){
				omega_t1<-omega_t1_exo
				f_model<-f_model_exo
			}
			rm(f_model_exo)

			persy_temp<-sum(table(c(omega_t1,omega_t))==2)/length(unique(c(omega_t,omega_t1)))
		}
	}
#### stare szukanie zbioru
		if (FALSE){
			f_model<-spikeslab(f_sent~f_gt_mod,n.iter1=5000,n.inter2=1000,ntree=750,intercept=TRUE,bigp.smalln = TRUE)
			omega_t1<-gsub("f_gt_mod","",rownames(f_model$summary)[round(f_model$summary$bma,4)!=0])
	
			if(i!=3){
				stare_kwery<-read.table(paste("./TempTXT/Kwerendy_",i-1,".txt",sep=""))
				if(identical(omega_t1,character(0))){omega_t1<-as.character(stare_kwery[,1])}
				Omega<-c(as.character(stare_kwery[,1]),rownames(f_model$summary[round(f_model$summary$bma,4)!=0,]))
				Omega<-gsub("f_gt_mod","",Omega)
				persy_temp<-sum(table(Omega)>1)/length(unique(Omega))
				omega_t1<-gsub("f_gt_mod","",rownames(f_model$summary[round(f_model$summary$bma,4)!=0,]))
				
				if(identical(omega_t1,character(0))){
					omega_t1<-as.character(stare_kwery[,1])
				} else{
				
				if(persy_temp<pers){
				
					req<-ceiling(length(unique(Omega))*pers)
					need<-req-sum(table(Omega)>1)
					ava<-nrow(stare_kwery)-(length(Omega)-length(unique(Omega)))

					if(ava>need){
						omega_t1<-c(
							rownames(f_model$summary[round(f_model$summary$bma,3)!=0,]),
							as.character(stare_kwery[order(-abs(stare_kwery[,2])),][!(gsub("f_gt_mod","",stare_kwery[order(-abs(stare_kwery[,2])),1]) %in% gsub("f_gt_mod","",rownames(f_model$summary[round(f_model$summary$bma,4)!=0,]))),][1:need,1])
						)
					
					} else {
						req<-floor(nrow(stare_kwery)*(pers^(-1)))
						need<-req-nrow(stare_kwery)
						ava<-f_model$summary[round(f_model$summary$bma,4)!=0,][!(gsub("f_gt_mod","",rownames(f_model$summary[round(f_model$summary$bma,4)!=0,])) %in% gsub("f_gt_mod","",stare_kwery[order(-abs(stare_kwery[,2])),1])),]
	
						omega_t1<-c(
							as.character(stare_kwery[order(-abs(stare_kwery[,2])),1]),
							rownames(ava[order(-abs(ava[,2])),][1:min(need,nrow(ava)),])
							)

					}}
				
					omega_t1<-gsub("f_gt_mod","",unique(omega_t1))
					Omega<-c(as.character(stare_kwery[,1]),omega_t1)
					Omega<-gsub("f_gt_mod","",Omega)
					persy_temp<-sum(table(Omega)>1)/length(unique(Omega))
				
					neu_gt<-ts(f_sent_gt[,omega_t1],start=c(1,1),freq=12)
					neu_gt_fory<-fory[,omega_t1]

					#print(persy_temp)
					PERSY<-c(PERSY,persy_temp)
				}

			} 
		}
		
		neu_gt<-ts(f_sent_gt[,omega_t1],start=c(1,1),freq=12)
		neu_gt_fory<-fory[,omega_t1]
		neu_gt_fory_mean<-foryMEAN[,omega_t1]

		f_gt_mod<-f_sent_gt[,omega_t1]
		f_model<-spikeslab(f_sent~f_gt_mod,n.iter1=1000,n.inter2=200,ntree=100,intercept=TRUE)
		
		if(szukaj){
			write.table(
				cbind(
					gsub("f_gt_mod","",rownames(f_model$summary)),
					f_model$summary$bma
				),paste("./TempTXT/Kwerendy_",i,".txt",sep=""))
		}

		#if (length(neu_gt)==0){
		#	f_model<-spikeslab(f_sent~f_sent_gt,n.iter1=1500,n.inter2=750,ntree=1000,intercept=TRUE,bigp.smalln = TRUE)
		#	neu_gt<-ts(f_sent_gt[,round(f_model$gnet,4)!=0],start=c(1,1),freq=12)
		#	neu_gt_fory<-fory[,round(f_model$gnet,4)!=0]	
		#}

		#if (sum(abs(f_model$gnet)>0.001)>5){
		#	neu_gt<-ts(f_sent_gt[,abs(f_model$gnet)>0.001],start=c(1,1),freq=12)
		#	neu_gt_fory<-fory[,abs(f_model$gnet)>0.001]
		#}

		nnet_hype<-ts(hype,freq=12,start=start(neu_gt))
		assign("nnet_hype",nnet_hype,envir = .GlobalEnv)
		exo_hype<-ts(hype[6:18],freq=12)
		exo_hype_wwuk<-ts(hype_wwuk[6:18],freq=12)
		exo_sent<-forecast(arima(f_sent,order=c(2,1,1),method="CSS"),h=12)$mean		
		exo_wwuk<-forecast(arima(f_wwuk,order=c(1,1,0),method="CSS"),h=12)$mean
		f_gt<-as.matrix(f_dane_gt[(i-1+p_lag):(i+dl-1),omega_t1])%*%as.matrix(f_model$gnet.scale)
		assign("f_gt",f_gt,envir = .GlobalEnv)
		exo_gt<-(as.matrix(foryAR[,omega_t1])%*%as.matrix(f_model$gnet.scale))

		#if (is.null(ncol(neu_gt))){
		#	neu_gt<-ts(rep(1,nrow(f_sent_gt)),start=c(1,1),freq=12)
		#	neu_gt_fory<-rep(1,nrow(f_sent_gt))
		#	exo_gt<-rep(1,nrow(f_sent_gt))
		#} else {
		#	if (length(neu_gt)==0){
		#		neu_gt<-ts(rep(1,nrow(f_sent_gt)),start=c(1,1),freq=12)
		#		neu_gt_fory<-rep(1,nrow(f_sent_gt))
		#		exo_gt<-rep(1,nrow(f_sent_gt))
		#	}
		#}
		
		cat("GOOGLE TRENDS: CHECKED \n")
		cat(paste("Number of Google Trends queries:",ncol(neu_gt)," \n"))
		if(i!=3){cat(paste("Persistency:",round(persy_temp,3)," \n"))}
		cat(paste("MSE of spike-and-slab:",round(f_model$sigma.hat,6)," \n"))
		

		#rm(f_model,foryAR)

###############################################################
			ARIMAS<-c(
				forecast(arima(f_rynek,order=c(0,0,0),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,1),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,2),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,0),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,1),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,2),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,0),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,1),method="CSS"),h=12)$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,2),method="CSS"),h=12)$mean[hor]
			)
			cat("ARIMA: CHECKED \n")

###############################################################
			ARSENTIS<-c(
				forecast(arima(f_rynek,order=c(0,0,0),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,1),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,2),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,0),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,1),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,2),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,0),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,1),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,2),xreg=data.frame(f_sent),method="CSS"),h=12,xreg=data.frame(exo_sent))$mean[hor]
			)
			cat("ARIMA with SENTIMENTS: CHECKED \n")
################################################################
			ARWUKIS<-c(
				forecast(arima(f_rynek,order=c(0,0,0),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,1),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(0,0,2),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,0),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,1),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(1,0,2),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,0),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,1),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor],
				forecast(arima(f_rynek,order=c(2,0,2),xreg=data.frame(f_wwuk),method="CSS"),h=12,xreg=data.frame(exo_wwuk))$mean[hor]
			)
			cat("ARIMA with WWUK: CHECKED \n")
################################################################	
#cumsum(f_sent)		
#(tail(cumsum(f_sent),1)+cumsum(exo_sent))
			LSTARY<-c(
				tryCatch({predict(lstar(f_rynek,m=2,d=1,mL=1,mH=1,thVar=hype,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=1,mH=2,thVar=hype,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=2,mH=1,thVar=hype,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=2,mH=2,thVar=hype,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))})		
			)
			cat("LSTAR SENT: CHECKED \n")

			LSTARYWWUK<-c(
				tryCatch({predict(lstar(f_rynek,m=2,d=1,mL=1,mH=1,thVar=hype_wwuk,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype_wwuk,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=1,mH=2,thVar=hype_wwuk,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype_wwuk,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=2,mH=1,thVar=hype_wwuk,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype_wwuk,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({predict(lstar(f_rynek,m=3,d=1,mL=2,mH=2,thVar=hype_wwuk,trace=FALSE,include="const",control=list(maxit=4500)),thVar=exo_hype_wwuk,n.ahead=12)[hor]},error=function(a){return(rep(Inf,4))})		
			)
			cat("LSTAR WWUK: CHECKED \n")
##############################################################
			GUGLE<-c(
				tryCatch({forecast(arima(f_rynek,order=c(0,0,0),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(0,0,1),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(0,0,2),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(1,0,0),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(1,0,1),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(1,0,2),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(2,0,0),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(2,0,1),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))}),
				tryCatch({forecast(arima(f_rynek,order=c(2,0,2),xreg=data.frame(f_gt),method="CSS"),h=12,xreg=data.frame(exo_gt))$mean[hor]},error=function(a){return(rep(Inf,4))})
			)
			cat("ARIMA with SPIKESLAB: CHECKED \n")
##############################################################	
rdzenie<-detectCores()-1
cl<-makeCluster(rdzenie)
WAGI_ARMANNET<-NULL
WAGI_NNET<-NULL
WAGI_NNET_MEAN<-NULL
for(k in 0:8){
###############################################################	

			modar<-arima(f_rynek,order=c(k %% 3,0,floor(k/3)),xreg=nnet_hype,method="CSS")
			res_nnet<-modar$res
	
			target<-cbind(
				as.numeric(res_nnet<qnorm(0.05,0,sd(res_nnet))),
				as.numeric(qnorm(0.05,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(0.15,0,sd(res_nnet))),
				as.numeric(qnorm(0.15,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(0.35,0,sd(res_nnet))),
				as.numeric(qnorm(0.35,0,sd(res_nnet))<=res_nnet & res_nnet<0),
				as.numeric(0<res_nnet & res_nnet<qnorm(0.65,0,sd(res_nnet))),
				as.numeric(qnorm(0.65,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(0.85,0,sd(res_nnet))),
				as.numeric(qnorm(0.85,0,sd(res_nnet))<=res_nnet & res_nnet<qnorm(0.95,0,sd(res_nnet))),
				as.numeric(qnorm(0.95,0,sd(res_nnet))<=res_nnet)
			)
			target<-target[,apply(target[1:47,],2,sum)>0]

			prewagi<-tryCatch({read.csv(paste("./TempNNET/ARMANET_",i,".csv",sep=""))},error=function(e){return(NULL)})
			presiec<-tryCatch({nnet(target[1:47,]~.,size=8,data=neu_gt[1:47,],trace=FALSE,softmax=TRUE,Wts=prewagi[,k+1])},error=function(e){return(NULL)})
			
			test_per<-round(seq(1,59,length.out=12),0)
			lear_per<-(1:59)[-test_per]
			
			wagi<-foreach(o=1:128,.combine=rbind)%dopar%{
				nnet(target[lear_per,]~.,size=8,data=neu_gt[lear_per,],trace=FALSE,softmax=TRUE)$wts
			}

			val<-foreach(o=1:128,.combine=rbind)%dopar%{
				rmse(((round(predict(
					nnet(target[lear_per,]~.,size=8,data=neu_gt[lear_per,],trace=FALSE,softmax=TRUE,Wts=wagi[o,])
				,neu_gt[test_per,]),1))%*%apply((target*c(res_nnet)),2,function(x){mean(x[x!=0])}))-res_nnet[test_per])
			}
			
			siec<-nnet(target~.,size=8,data=neu_gt,trace=FALSE,softmax=TRUE,Wts=wagi[which.min(val),],maxit=250)
			presiecRMSE<-tryCatch({rmse(((round(predict(nnet(target[lear_per,]~.,size=10,data=neu_gt[lear_per,],trace=FALSE,softmax=TRUE,Wts=prewagi[,k+1])
							,neu_gt[test_per,]),1))%*%apply((target*c(res_nnet)),2,function(x){mean(x[x!=0])}))-res_nnet[test_per])},
					error=function(e){return(Inf)})
			rm(prewagi)
			if(presiecRMSE<min(val)){
				siec<-presiec
			}
			
			WAGI_ARMANNET<-cbind(WAGI_ARMANNET,siec$wts)			

			ARNEUNET<-c(ARNEUNET,
					forecast(modar,xreg=data.frame(nnet_hype=exo_hype),h=12)$mean[hor]+(predict(siec,neu_gt_fory)%*%apply((target*c(res_nnet)),2,function(x){mean(x[x!=0])}))[hor]
				)
			rm(res_nnet,siec,modar,wagi,val)
###############################################################
if(FALSE){
			nnetdata<-tail(na.omit(cbind(lag(f_rynek,-1),lag(f_rynek,-2),lag(f_rynek,-3))),-1)
	
			wagi<-foreach(o=1:8,.combine=rbind)%dopar%{
				nnet(tail(f_rynek,-3)[1:45]~.,size=4,data=nnetdata[1:45,],linout=TRUE,trace=FALSE)$wts
			}

			val<-foreach(o=1:8,.combine=rbind)%dopar%{
				rmse(predict(
					nnet(tail(f_rynek,-3)~.,size=4,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[o,])
				,nnetdata[46:56,])-tail(f_rynek,-3)[46:56])
			}

			siec<-nnet(tail(f_rynek,-3)~.,size=4,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[which.max(val),])
			
			rm(temp_siec)
			ff<-forecast(arima(f_rynek,order=c(k %% 3,0,floor(k/3)),method="CSS"),h=12)$mean
			fff<-cbind(c(tail(f_rynek,1),ff[1:11]),c(tail(f_rynek,2),ff[1:10]),c(tail(f_rynek,3),ff[1:9]))
			colnames(fff)<-colnames(nnetdata)
			NEUNET<-c(NEUNET,(predict(siec,fff))[hor])
			rm(ff,fff,nnetdata,siec)
}
################################################################
			
			nnetdata<-na.omit(cbind(lag(f_rynek,-1),lag(f_rynek,-2),lag(f_rynek,-3),neu_gt,nnet_hype))

			prewagi<-tryCatch({read.csv(paste("./TempNNET/ANN_",i,".csv",sep=""))},error=function(e){return(NULL)})
			presiec<-tryCatch({nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=prewagi[,k+1])},error=function(e){return(NULL)})

			test_per<-round(seq(1,56,length.out=12),0)
			lear_per<-(1:56)[-test_per]

			wagi<-foreach(o=1:128,.combine=rbind)%dopar%{
				nnet(tail(f_rynek,-3)[lear_per]~.,size=8,data=nnetdata[lear_per,],linout=TRUE,trace=FALSE)$wts
			}

			val<-foreach(o=1:128,.combine=rbind)%dopar%{
				rmse(predict(
					nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[o,])
				,nnetdata[test_per,])-tail(f_rynek,-3)[test_per])
			}

			siec<-nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[which.min(val),])			
			presiecRMSE<-tryCatch({rmse(predict(
					nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=prewagi[,k+1])
				,nnetdata[test_per,])-tail(f_rynek,-3)[test_per])},
					error=function(e){return(Inf)})
			rm(prewagi)
			if(presiecRMSE<min(val)){
				siec<-presiec
			}
			
			WAGI_NNET<-cbind(WAGI_NNET,siec$wts)	
			
			ff<-forecast(arima(f_rynek,order=c(k %% 3,0,floor(k/3)),method="CSS"),h=12)$mean
			fff<-cbind(c(tail(f_rynek,1),ff[1:11]),c(tail(f_rynek,2),ff[1:10]),c(tail(f_rynek,3),ff[1:9]),neu_gt_fory,c(exo_hype[1:12]))
			colnames(fff)<-colnames(nnetdata)
			NEUNETGT<-c(NEUNETGT,(predict(siec,fff))[hor])
			rm(ff,fff,net_wagi,net_blad,nnetdata,siec)

################################################################
			
			nnetdata<-na.omit(cbind(lag(f_rynek,-1),lag(f_rynek,-2),lag(f_rynek,-3),neu_gt,nnet_hype))

			prewagi<-tryCatch({read.csv(paste("./TempNNET/ANNMEAN_",i,".csv",sep=""))},error=function(e){return(NULL)})
			presiec<-tryCatch({nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=prewagi[,k+1])},error=function(e){return(NULL)})

			test_per<-round(seq(1,56,length.out=12),0)
			lear_per<-(1:56)[-test_per]

			wagi<-foreach(o=1:128,.combine=rbind)%dopar%{
				nnet(tail(f_rynek,-3)[lear_per]~.,size=8,data=nnetdata[lear_per,],linout=TRUE,trace=FALSE)$wts
			}

			val<-foreach(o=1:128,.combine=rbind)%dopar%{
				rmse(predict(
					nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[o,])
				,nnetdata[test_per,])-tail(f_rynek,-3)[test_per])
			}

			siec<-nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=wagi[which.min(val),])
			
			presiecRMSE<-tryCatch({rmse(predict(
					nnet(tail(f_rynek,-3)~.,size=8,data=nnetdata,linout=TRUE,trace=FALSE,Wts=prewagi[,k+1])
				,nnetdata[test_per,])-tail(f_rynek,-3)[test_per])},
					error=function(e){return(Inf)})
			rm(prewagi)
			if(presiecRMSE<min(val)){
				siec<-presiec
			}
			
			WAGI_NNET_MEAN<-cbind(WAGI_NNET_MEAN,siec$wts)	
			
			ff<-forecast(arima(f_rynek,order=c(k %% 3,0,floor(k/3)),method="CSS"),h=12)$mean
			fff<-cbind(c(tail(f_rynek,1),ff[1:11]),c(tail(f_rynek,2),ff[1:10]),c(tail(f_rynek,3),ff[1:9]),neu_gt_fory_mean,c(exo_hype[1:12]))
			colnames(fff)<-colnames(nnetdata)
			NEUNETGTMEAN<-c(NEUNETGTMEAN,(predict(siec,fff))[hor])
			rm(ff,fff,net_wagi,net_blad,nnetdata,siec)
###############################################################
}
stopCluster(cl)
cat("NNETs: CHECKED \n")
write.csv(WAGI_ARMANNET,paste("./TempNNET/ARMANET_",i,".csv",sep=""),row.names=FALSE)
write.csv(WAGI_NNET,paste("./TempNNET/ANN_",i,".csv",sep=""),row.names=FALSE)
write.csv(WAGI_NNET_MEAN,paste("./TempNNET/ANNMEAN_",i,".csv",sep=""),row.names=FALSE)

###############################################################
			VARY<-c(
				forecast(VAR(cbind(f_rynek,f_sent),p=1),h=12)$forecast$f_rynek$mean[hor],
				forecast(VAR(cbind(f_rynek,f_sent),p=2),h=12)$forecast$f_rynek$mean[hor],
				forecast(VAR(cbind(f_rynek,f_sent),p=3),h=12)$forecast$f_rynek$mean[hor]
				)
			cat("VAR: CHECKED \n")
###############################################################
		#ARSENTIS<-rep(0,9*4)
		#ARWUKIS<-rep(0,9*4)
		#GUGLE<-rep(0,9*4)
		#LSTARY<-rep(0,4*4)
		#LSTARYWWUK<-rep(0,4*4)
		#ARNEUNET<-rep(0,9*4)
		NEUNET<-rep(0,9*4)
		#NEUNETGT<-rep(0,9*4)
		#NEUNETGTMEAN<-rep(0,9*4)
		#VARY<-rep(0,3*4)

		
		write.csv(rbind(c(tail(f_dane_var[,2],1),0,
			ARIMAS,
			ARSENTIS,
			ARWUKIS,
			LSTARY,
			LSTARYWWUK,	
			GUGLE,
			ARNEUNET,
			NEUNET,
			NEUNETGT,
			NEUNETGTMEAN,
			VARY
		)),paste("./Temp/Forecast",i,".csv",sep=""))	
		
		cat(paste("SAVED:",i,"\n"))
		rm(ARIMAS,ARSENTIS,ARWUKIS,LSTARY,GUGLE,ARNEUNET,NEUNET,NEUNETGT,VARY)
	
		#write.table(colnames(neu_gt),paste("./TempTXT/Kwerendy_",i,".txt",sep=""))
	}

	for (i in okres){
		OUT_PROGNOZA<-rbind(OUT_PROGNOZA,c(read.csv(paste("./Temp/Forecast",i,".csv",sep=""))))
	}
	return(list(prognoza=OUT_PROGNOZA[,-1],kwerendy=KWERY,pers=PERSY,senty=MEGASENT))
}

plot_senty<-function(f_dane,dl){

	u_oS<-matrix(0,4,4)
	u_oS[1,c(2)]<-NA
	u_oS[2,c(1)]<-NA
	u_oS[3,c(1,4)]<-NA
	u_oS[4,c(2,3)]<-NA
	diag(u_oS)<-1
	
	B_oS<-matrix(0,4,4)
	diag(B_oS)<-NA
	
	TOPLOT<-array(NA,dim=c(57,118))

	for (i in 3:57){
		plsnt<-rep(NA,117)
		f_dane_var<-ts(f_dane[(i-1):(i+dl-1),],start=c(1,1),freq=12)
		f_rynek<-ts(f_dane_var[-(1:p_lag),2],start=c(1,1),freq=12)
		f_wwuk<-ts(f_dane_var[-(1:p_lag),3],start=c(1,1),freq=12)
		f_sent<-ts(c(get_sents(f_dane_var[,-c(3,6,7)])$sent),start=c(1,1),freq=12)
		TOPLOT[i,i:(i+58)]<-f_sent
	}

	TOPLOT<-ts(t(TOPLOT),start=c(2007,5),freq=12)

	plot(-TOPLOT[,3],type="l",ylim=c(-3,1.5),ylab="Sentiment",xlab="Year")
	for(i in 4:57){
		lines(-TOPLOT[,i])
	}
}

my_dm.test<-function(x,y,hor){
	zwroc<-NULL
		zwroc<-c(zwroc,dm.test(na.omit(x[x[,i]<Inf,i]),na.omit(y[x[,i]<Inf,i]),h=hor[i],alternative="less")$p.value)
	return(zwroc)
}

analiza.spektralna<-function(f_sent){
	dev.new()
	x.spec <- spectrum(f_sent,span=12,log="no",plot=FALSE)
	spx <- x.spec$freq
	spx<-spx/mean(diff(spx))
	spy <- 2*x.spec$spec
	plot(spy~spx,xlab="frequency",ylab="density",main="smoothed spectral density",type="l")
	abline(v=which.max(spy),col="blue",lty=4)
	mtext(
		which.max(spy),
		3,
		adj=(which.max(spy)+3)/max(spx),col="blue")
}

fory.Google<-function(f_dane_gg){

		f_dane.gt_fory<-NULL
		for(fi in 1:(ncol(f_dane_gg)/7)){
			f_dane.gt_fory<-cbind(f_dane.gt_fory,forecast(arima(f_dane_gg[,fi],order=ord_for_gg[[fi]],method="CSS"),h=12)$mean)
		}
		f_dane.gt_fory<-cbind(
			f_dane.gt_fory,
			rbind(
				f_dane_gg[nrow(f_dane_gg),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:11,]
			),
			rbind(
				f_dane_gg[(nrow(f_dane_gg)-1):(nrow(f_dane_gg)),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:10,]
			),
			rbind(
				f_dane_gg[(nrow(f_dane_gg)-2):(nrow(f_dane_gg)),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:9,]
			),
			rbind(
				f_dane_gg[(nrow(f_dane_gg)-3):(nrow(f_dane_gg)),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:8,]
			),
			rbind(
				f_dane_gg[(nrow(f_dane_gg)-4):(nrow(f_dane_gg)),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:7,]
			),
			rbind(
				f_dane_gg[(nrow(f_dane_gg)-5):(nrow(f_dane_gg)),1:(ncol(f_dane_gg)/7)],f_dane.gt_fory[1:6,]
			)
		)
	return(f_dane.gt_fory)
}

Type.Combine<-function(f_IS,f_type){
	return(
		table(merge(substr(f_IS,1,nchar(f_IS)-2),f_type,by.x=1,by.y=1)[,2])
	)
}
cor.dist<-function(f_data){

	RESULT<-matrix(0,nrow=ncol(f_data),ncol=ncol(f_data))
	sumn<-apply(f_data,1,sum)
	summ<-apply(f_data,2,sum)
	sumt<-sum(f_data)

	for(n in 1:ncol(f_data)){
		for(m in 1:n){
			RESULT[n,m]<-(-1)*(((cor(f_data[,n],f_data[,m],use="pairwise.complete.obs")+1)/2)-1)
		}
	}
	RESULT<-RESULT+t(RESULT)
	diag(RESULT)<-0
	return(RESULT)
}
sq.chi.dist<-function(f_data,type="chi"){
	if(!(type %in% c('chi','phi'))){return("Not recognized type")}

	RESULT<-matrix(0,nrow=nrow(f_data),ncol=ncol(f_data))
	sumn<-apply(f_data,1,sum)
	summ<-apply(f_data,2,sum)
	sumt<-sum(f_data)

	for(n in 1:nrow(f_data)){
		for(m in 1:n){
			ff_data<-f_data[,c(n,m)]

			sumn<-apply(ff_data,1,sum)
			summ<-apply(ff_data,2,sum)
			sumt<-sum(ff_data)

			xexpn<-sumn*summ[1]/sumt
			xexpm<-sumn*summ[2]/sumt
		
			partn<-(f_data[,n]-xexpn)^2/xexpn
			partm<-(f_data[,m]-xexpm)^2/xexpm
		
			partn[xexpn==0]<-0
			partm[xexpm==0]<-0

			RESULT[n,m]<-sqrt(sum(partn+partm))
			if(type=="phi"){RESULT[n,m]<-RESULT[n,m]/sqrt(sumt)}
		}
	}
	diag(RESULT)<-0
	RESULT<-RESULT+t(RESULT)
	return(RESULT)
}

cluster.dist<-function(f_data){

	RESULT<-matrix(0,nrow=ncol(f_data),ncol=ncol(f_data))
	sumn<-apply(f_data,1,sum)
	summ<-apply(f_data,2,sum)
	sumt<-sum(f_data)

	for(n in 1:ncol(f_data)){
		for(m in 1:n){
			RESULT[n,m]<-(1-(sum(f_data[,n]==f_data[,m])/nrow(f_data)))
		}
	}
	RESULT<-RESULT+t(RESULT)
	diag(RESULT)<-0
	return(RESULT)
}

irf_closed<-function(model_svar,impulse,response,closed_channels,Kmax=24){
	
	if(class(model_svar)!="svarest"){
		print("model in not from class 'svarest'")
		break
	}
	
	model <- model_svar$var
	Nimpulse <- colnames(model$y)[impulse]
	Nresponse <- colnames(model$y)[response]

	IRF <- IRF.closed <-rep(NaN,Kmax)

	F <- matrix(unlist(Acoef(model)),nrow=(model$K))
	F <- rbind(F
          	,cbind(diag(model$p-1) %x% diag((model$K)),matrix(0,nrow=(model$K)*(model$p-1),ncol=(model$K)))
  	)


	F.closed <- F
	for (cc in closed_channels){
		F.closed[response,seq(0,(model$K)*(model$p-1),by=model$K)+cc] <- 0
	}

	A <- model_svar$A
	A.closed <- model_svar$A
	A.closed[response,closed_channels] <- 0

	temp <- diag(1,(model$K))
	temp <- solve(A)%*%model_svar$B
	temp <- rbind(cbind(temp,matrix(0,nrow=(model$K),ncol=(model$K)*(model$p-1)))
		,matrix(0,nrow=(model$K)*(model$p-1),ncol=(model$K)*(model$p))
	)

	temp.closed <- temp

	for(k in 1:Kmax){
 		IRF[k] <- temp[response,impulse]
		IRF.closed[k] <- temp.closed[response,impulse]
		temp <- F %*% temp
 		temp.closed <- F.closed%*%temp.closed
	}

	return(list(
		model=irf(model,impulse=Nimpulse,response=Nresponse,n.ahead=Kmax)$irf[[1]]
		,calc=IRF
		,Psi=Psi(model,nstep=(Kmax-1))[response,impulse,]
		,closed=IRF.closed
	))
}

	



