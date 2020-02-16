##### plot_closed_irf ##### 
#
plot_closed_irf<-function(model,title=""){
  
  plot.ts(model$model
          ,ylim=c(min(unlist(model)),max(unlist(model)))
          ,lwd=2
          ,main=title)
  
  lines(model$calc,col="red",lwd=2,lty=2)
  
  lines(model$closed,col="blue",lwd=2)
  
  lines(model$Psi,col="orange",lty=4)
  
  legend("topright"
         ,legend=c("vars::irf","calc IRF","IRF closed channel","Psi")
         ,lty=c(1,1,1,4)
         ,pch=c(NA,NA,NA,NA)
         ,col=c("black","red","blue","orange")
         ,bg="white"
         ,horiz=FALSE
         ,xpd=TRUE
         ,cex = 0.7)
}
