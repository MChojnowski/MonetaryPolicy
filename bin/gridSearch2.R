
gridSearch2 <- function(data,ph,pl,sent,mean.sent=mean(sent),sd.sent=sd(sent),min.observ=60,grid.length=50,prob.lim=c(0,1),stable_param=1){

  #Models  
  var_dflt <- vars::VAR(data,p=max(ph,pl))
  var_dflt_h <- vars::VAR(data,p=ph)
  var_dflt_l <- vars::VAR(data,p=pl)
  
  # Parameters
  gridM <- qnorm(seq(0.01,0.99,length.out=grid.length),mean.sent,sd.sent)
  gridS <- seq(0.05,3,length.out=grid.length)

  
  if(min.observ<=1){
    min.observ <- min.observ*nrow(var_dflt_h$y)
  }
  
  # Objects

  
  mat.SSE <- NULL
  mat.AIC <- NULL
  mat.ROOTS <- NULL
  mat.IRF <- list(NULL,NULL)
  
  pb <- txtProgressBar(min = 0, max = length(gridM)*length(gridS), style = 3)
  
  for(FctM in gridM){
    
    sum.SSE <- NULL
    sum.AIC <- NULL
    sum.ROOTS <- NULL
    sum.IRF1 <- NULL
    sum.IRF2 <- NULL
    
    for(FctS in gridS){
      
      setTxtProgressBar(pb, ((which(FctM==gridM)-1)*length(gridS)+which(FctS==gridS)))
      
      prob <- pnorm(sent,FctM,FctS)
      
      if(length(prob[round(prob,2)>prob.lim[1] & round(prob,2)<prob.lim[2]])>min.observ){

        stsvar_model <- STSVAR(data,ph,pl,as.vector(prob))

        var_dflt_h <- stsvar_model$varH
        var_dflt_l <- stsvar_model$varL
        
        checkpoint.reszta<-list(apply(residuals(var_dflt_h),2,mean),apply(residuals(var_dflt_l),2,mean))
        
        if(var_dflt_h$K == var_dflt_l$K){
          VarFitted<-prob*ts(fitted(var_dflt_h),freq=12,end=end(var_dflt_h$y)) +
            (1-prob)*ts(fitted(var_dflt_l),freq=12,end=end(var_dflt_l$y))
        } else {
          VarFitted<-prob*ts(cbind(fitted(var_dflt_h),as.matrix(array(0,dim=c(nrow(fitted(var_dflt_h)),sum(!colnames(dane_all_org) %in% colnames(fitted(var_dflt_h))))))),freq=12,end=end(var_dflt_h$y)) +
            (1-prob)*ts(cbind(fitted(var_dflt_l),as.matrix(array(0,dim=c(nrow(fitted(var_dflt_h)),sum(!colnames(dane_all_org) %in% colnames(fitted(var_dflt_l))))))),freq=12,end=end(var_dflt_l$y))
        }
        
        colnames(VarFitted)<-colnames(var_dflt$y)
        
        resTempVar<-(VarFitted-var_dflt$y)
        
        resTempVar[,1]<-resTempVar[,1]/(apply(dane_monet,2,sd)[1])
        resTempVar[,2]<-resTempVar[,2]/(apply(dane_monet,2,sd)[2])
        resTempVar[,3]<-resTempVar[,3]/(apply(dane_monet,2,sd)[3])
        resTempVar[,4]<-resTempVar[,4]/(apply(dane_monet,2,sd)[4])
        
        sum.SSE<-c(sum.SSE,sqrt(sum(resTempVar^2)))
        sum.ROOTS<-c(sum.ROOTS,max(roots(var_dflt_h),roots(var_dflt_l)))
        sum.AIC<-rbind(sum.AIC,c(AIC(var_dflt_h),AIC(var_dflt_l)))
        sum.IRF1<-c(sum.IRF1,-1)
        sum.IRF2<-c(sum.IRF2,-1)
        
      }else{
        
        sum.SSE<-c(sum.SSE,Inf)
        sum.ROOTS<-c(sum.ROOTS,Inf)
        sum.AIC<-rbind(sum.AIC,c(Inf,Inf))
        sum.IRF1<-rbind(sum.IRF1,c(Inf,Inf))
        sum.IRF2<-rbind(sum.IRF2,c(Inf,Inf))
      }
      
      
    }
    
    mat.SSE<-cbind(mat.SSE,sum.SSE)
    mat.ROOTS<-cbind(mat.ROOTS,sum.ROOTS)
    mat.AIC<-cbind(mat.AIC,apply(sum.AIC,1,sum))
    mat.IRF[[1]]<-cbind(mat.IRF[[1]],sum.IRF1)
    mat.IRF[[2]]<-cbind(mat.IRF[[2]],sum.IRF1)
    
  }
  
  cord_sse <- mat.SSE
  cord_sse[mat.ROOTS >= stable_param] <- Inf
  cord_sse <- which(cord_sse == min(cord_sse), arr.ind = TRUE)
  
  return(list(gridM=gridM
              ,gridS=gridS
              ,sse=mat.SSE
              ,roots=mat.ROOTS
              ,aic=mat.AIC
              ,irf=mat.IRF
              ,sent=sent
              ,optparam=c(gridS[cord_sse[1]],gridM[cord_sse[2]])
              ,stable_param=stable_param
  ))
}
