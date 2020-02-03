gridSearch<-function(var_dflt,sent,type,mean.sent=mean(sent),sd.sent=sd(sent)){

	if(type !="mean" & type !="sd"){
		cat("Wrong type! Please input either 'mean' or 'sd'")
		break
	}

	if(type=="sd"){
		grid<-seq(0.005,3,by=0.005)
	}

	if(type=="mean"){
		grid<-qnorm(seq(0.001,0.999,0.001),mean.sent,sd.sent)
	}


	sum.SSE<-NULL
	sum.AIC<-NULL
	sum.ROOTS<-NULL

	pb <- txtProgressBar(min = 0, max = length(grid), style = 3)

	for(Fct in grid){

		setTxtProgressBar(pb, which(Fct==grid)-1)

		var1<-var2<-var_dflt

		#prob<-abs(pnorm(sent,bestMeanFct,sd(qwe1)*sdFct)-0.5)*2
	
		if(type=="sd"){
			prob<-pnorm(sent,mean.sent,sd(qwe1)*Fct)
		}
	
		if(type=="mean"){
			prob<-pnorm(sent,Fct,sd.sent)
		}

		dane_monet1<-var_dflt$y*prob
		dane_monet2<-var_dflt$y*(1-prob)

		colnames(dane_monet1)<-c("GDP","CPI","WIBOR","REER")
		colnames(dane_monet2)<-c("GDP","CPI","WIBOR","REER")

		dane_all<-dane_monet
		for (i in 1:var_dflt$p){
  			dane_all<-cbind(dane_all,lag(dane_monet1,-i))
		}
		for (i in 1:var_dflt$p){
 			dane_all<-cbind(dane_all,lag(dane_monet2,-i))
		}

		r1<-lm(dane_all[,1]~dane_all[,-c(1:4)]-1)
		r2<-lm(dane_all[,2]~dane_all[,-c(1:4)]-1)
		r3<-lm(dane_all[,3]~dane_all[,-c(1:4)]-1)
		r4<-lm(dane_all[,4]~dane_all[,-c(1:4)]-1)

		comat1<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c(1:(4*var_dflt$p))]
		comat2<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c((4*var_dflt$p+1):length(r1$coef))]

		for (i in 1:4){
  			nazwy<-names(var1$varresult[[i]]$coefficients)
  			var1$varresult[[i]]$coefficients<-c(comat1[i,],0)  
  			names(var1$varresult[[i]]$coefficients)<-nazwy
  
  			nazwy<-names(var2$varresult[[i]]$coefficients)
  			var2$varresult[[i]]$coefficients<-c(comat2[i,],0)
  			names(var2$varresult[[i]]$coefficients)<-nazwy
  		}

		### DODAJ CONSTANT
		dane<-var_dflt$datamat[,-c(1:4,ncol(var_dflt$datamat))]

		#przeliczamy reszty
 		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-as.matrix(dane)%*%as.matrix(head(var1$varresult[[i]]$coefficients,-1))
    			var1$varresult[[i]]$residuals<-var_dflt$datamat[,i]-var1$varresult[[i]]$fitted.values

    			var2$varresult[[i]]$fitted.values<-as.matrix(dane)%*%as.matrix(head(var2$varresult[[i]]$coefficients,-1))
    			var2$varresult[[i]]$residuals<-var_dflt$datamat[,i]-var2$varresult[[i]]$fitted.values
  		}

		#Obliczamy stale
		c1<-residuals(var1)
		const1<-apply(c1,2,mean)
		#const1

		c2<-residuals(var2)
		const2<-apply(c2,2,mean)
		#const2


		### podmieniamy sta³e
		for (i in 1:4){
  			var1$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const1[i]
  			var2$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const2[i]
		}

  		#przeliczamy reszty uwzgledniajac stala
  		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-var1$varresult[[i]]$fitted.values+tail(var1$varresult[[i]]$coefficients,1)
    			var1$varresult[[i]]$residuals<-var1$varresult[[i]]$residuals-tail(var1$varresult[[i]]$coefficients,1)
    			
			var2$varresult[[i]]$fitted.values<-var2$varresult[[i]]$fitted.values+tail(var2$varresult[[i]]$coefficients,1)
    			var2$varresult[[i]]$residuals<-var2$varresult[[i]]$residuals-tail(var2$varresult[[i]]$coefficients,1)
  		}

		checkpoint.reszta<-list(apply(residuals(var1),2,mean),apply(residuals(var2),2,mean))

		VarFitted<-prob*ts(fitted(var1),freq=12,end=end(var_dflt$y)) +
			(1-prob)*ts(fitted(var2),freq=12,end=end(var_dflt$y))

		colnames(VarFitted)<-colnames(var_dflt$y)

		resTempVar<-(VarFitted-var_dflt$y)

		resTempVar[,1]<-resTempVar[,1]/apply(dane_monet,2,sd)[1]
		resTempVar[,2]<-resTempVar[,2]/apply(dane_monet,2,sd)[2]
		resTempVar[,3]<-resTempVar[,3]/apply(dane_monet,2,sd)[3]
		resTempVar[,4]<-resTempVar[,4]/apply(dane_monet,2,sd)[4]

		sum.SSE<-c(sum.SSE,sqrt(sum(resTempVar^2)))
		sum.ROOTS<-rbind(sum.ROOTS,c(max(roots(var1)),max(roots(var2))))
		sum.AIC<-rbind(sum.AIC,c(AIC(var1),AIC(var2)))

	}

	return(list(grid=grid
		,sse=sum.SSE
		,roots=sum.ROOTS
		,aic=sum.AIC
		))
}

gridSearch2D<-function(var_dflt,sent,type="all",mean.sent=mean(sent),sd.sent=sd(sent)){

	if(type !="mean" & type !="sd"){
		cat("Wrong type! Please input either 'mean' or 'sd'")
		break
	}

	if(type=="sd"){
		grid<-seq(0.005,3,by=0.005)
	}

	if(type=="mean"){
		grid<-qnorm(seq(0.001,0.999,0.001),mean.sent,sd.sent)
	}

	if(type=="all"){
		gridM<-qnorm(seq(0.001,0.999,0.001),mean.sent,sd.sent)
		gridS<-seq(0.005,3,by=0.005)
	}


	mat.SSE<-NULL
	mat.AIC<-NULL
	sum.ROOTS<-NULL
	pb <- txtProgressBar(min = 0, max = length(gridM)*length(gridS), style = 3)

	for(FctM in gridM){
		sum.SSE<-NULL
		sum.AIC<-NULL
		sum.ROOTS<-NULL

		for(FctS in gridS){

		setTxtProgressBar(pb, (which(FctM==gridM)-1*length(gridS)+which(FctS==gridS))

		var1<-var2<-var_dflt

		#prob<-abs(pnorm(sent,bestMeanFct,sd(qwe1)*sdFct)-0.5)*2
	
		if(type=="sd"){
			prob<-pnorm(sent,mean.sent,sd(qwe1)*Fct)
		}
	
		if(type=="mean"){
			prob<-pnorm(sent,Fct,sd.sent)
		}

		if(type=="mean"){
			prob<-pnorm(sent,FctM,FctS)
		}
		dane_monet1<-var_dflt$y*prob
		dane_monet2<-var_dflt$y*(1-prob)

		colnames(dane_monet1)<-c("GDP","CPI","WIBOR","REER")
		colnames(dane_monet2)<-c("GDP","CPI","WIBOR","REER")

		dane_all<-dane_monet
		for (i in 1:var_dflt$p){
  			dane_all<-cbind(dane_all,lag(dane_monet1,-i))
		}
		for (i in 1:var_dflt$p){
 			dane_all<-cbind(dane_all,lag(dane_monet2,-i))
		}

		r1<-lm(dane_all[,1]~dane_all[,-c(1:4)]-1)
		r2<-lm(dane_all[,2]~dane_all[,-c(1:4)]-1)
		r3<-lm(dane_all[,3]~dane_all[,-c(1:4)]-1)
		r4<-lm(dane_all[,4]~dane_all[,-c(1:4)]-1)

		comat1<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c(1:(4*var_dflt$p))]
		comat2<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c((4*var_dflt$p+1):length(r1$coef))]

		for (i in 1:4){
  			nazwy<-names(var1$varresult[[i]]$coefficients)
  			var1$varresult[[i]]$coefficients<-c(comat1[i,],0)  
  			names(var1$varresult[[i]]$coefficients)<-nazwy
  
  			nazwy<-names(var2$varresult[[i]]$coefficients)
  			var2$varresult[[i]]$coefficients<-c(comat2[i,],0)
  			names(var2$varresult[[i]]$coefficients)<-nazwy
  		}

		### DODAJ CONSTANT
		dane<-var_dflt$datamat[,-c(1:4,ncol(var_dflt$datamat))]

		#przeliczamy reszty
 		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-as.matrix(dane)%*%as.matrix(head(var1$varresult[[i]]$coefficients,-1))
    			var1$varresult[[i]]$residuals<-var_dflt$datamat[,i]-var1$varresult[[i]]$fitted.values

    			var2$varresult[[i]]$fitted.values<-as.matrix(dane)%*%as.matrix(head(var2$varresult[[i]]$coefficients,-1))
    			var2$varresult[[i]]$residuals<-var_dflt$datamat[,i]-var2$varresult[[i]]$fitted.values
  		}

		#Obliczamy stale
		c1<-residuals(var1)
		const1<-apply(c1,2,mean)
		#const1

		c2<-residuals(var2)
		const2<-apply(c2,2,mean)
		#const2


		### podmieniamy sta³e
		for (i in 1:4){
  			var1$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const1[i]
  			var2$varresult[[i]]$coefficients[4*var_dflt$p+1]<-const2[i]
		}

  		#przeliczamy reszty uwzgledniajac stala
  		for (i in 1:4){
    			var1$varresult[[i]]$fitted.values<-var1$varresult[[i]]$fitted.values+tail(var1$varresult[[i]]$coefficients,1)
    			var1$varresult[[i]]$residuals<-var1$varresult[[i]]$residuals-tail(var1$varresult[[i]]$coefficients,1)
    			
			var2$varresult[[i]]$fitted.values<-var2$varresult[[i]]$fitted.values+tail(var2$varresult[[i]]$coefficients,1)
    			var2$varresult[[i]]$residuals<-var2$varresult[[i]]$residuals-tail(var2$varresult[[i]]$coefficients,1)
  		}

		checkpoint.reszta<-list(apply(residuals(var1),2,mean),apply(residuals(var2),2,mean))

		VarFitted<-prob*ts(fitted(var1),freq=12,end=end(var_dflt$y)) +
			(1-prob)*ts(fitted(var2),freq=12,end=end(var_dflt$y))

		colnames(VarFitted)<-colnames(var_dflt$y)

		resTempVar<-(VarFitted-var_dflt$y)

		resTempVar[,1]<-resTempVar[,1]/apply(dane_monet,2,sd)[1]
		resTempVar[,2]<-resTempVar[,2]/apply(dane_monet,2,sd)[2]
		resTempVar[,3]<-resTempVar[,3]/apply(dane_monet,2,sd)[3]
		resTempVar[,4]<-resTempVar[,4]/apply(dane_monet,2,sd)[4]

		sum.SSE<-c(sum.SSE,sqrt(sum(resTempVar^2)))
		sum.ROOTS<-rbind(sum.ROOTS,c(max(roots(var1)),max(roots(var2))))
		sum.AIC<-rbind(sum.AIC,c(AIC(var1),AIC(var2)))

	}

mat.SSE<-cbind(mat.SSE,sum.SSE)
mat.ROOTS<-cbind(mat.ROOTS,sum.ROOTS)
mat.AIC<-cbind(mat.AIC,sum.AIC)
}

	return(list(gridM=gridM
		,gridS=gridS
		,sse=mat.SSE
		,roots=mat.ROOTS
		,aic=mat.AIC
		))
}