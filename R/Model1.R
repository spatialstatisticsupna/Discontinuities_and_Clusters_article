model.selection <- function(cluster.prior, Y.real, E.real, C, max.cluster=NULL, final.cluster=NULL,
                            plot.dic=TRUE, lincomb=FALSE, extraconstr=FALSE, strategy="simplified.laplace"){
                            
     #### Fit a sequence of models with upto "max.cluster" and choose the best
     #### by minimising the Deviance Information Criterion.
     
time <- system.time({

     ## Define our uniform prior distributions ##
     sdunif="expression:
      logdens=-log_precision/2;
      return(logdens)"
  
     lunif = "expression:
      a = 1;
      b = 1;
      beta = exp(theta)/(1+exp(theta));
      logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
      log_jacobian = log(beta*(1-beta));
      return(logdens+log_jacobian)"
  
     R.Leroux <- diag(dim(C)[1])-C
  
     
     if(is.null(final.cluster)){
       
       ## Store the DIC values
       n <- length(Y.real)
       if(is.null(max.cluster)) {max.cluster <- n}
       dic.list <- rep(NA, max.cluster)
       waic.list <- rep(NA, max.cluster)

     
       ## Fit a model with a single cluster
       i <- 1
       
       data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n)

       formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                             hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
     
       model <- inla(formula, family="poisson", data=data.temp, E=E.real,
                    control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                    control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE),
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.inla=list(strategy = strategy, npoints = 21))

       dic.list[i] <- model$dic$dic
       waic.list[i] <- model$waic$waic
       
       if(plot.dic==TRUE) {
         win.graph()
         y.lim <- dic.list[i]*c(1/1.01,1.005)
         plot(1:max.cluster, dic.list, type="l", xlab="Number of clusters", ylab="DIC",
              main=paste("Model",i,"of",max.cluster), ylim=y.lim)
         points(1,dic.list[i],pch=19)
         mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2)
       }
     
       ## Fit separate models with between 2 and the max number of clusters.
       for(i in 2:max.cluster)
       {
         j <- n+1-i
         factor.clust <- cluster.prior$cluster.store[j,]
         data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n, factor.clust=factor.clust)
         formula <- Y.real ~ factor(factor.clust) +
                             f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                               hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
         model <- inla(formula, family="poisson", data=data.temp, E=E.real,
                       control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                       control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE),
                       control.predictor=list(compute=TRUE, cdf=c(log(1))),
                         control.inla=list(strategy = strategy, npoints = 21))
         dic.list[i] <- model$dic$dic
         waic.list[i] <- model$waic$waic
          
         if(plot.dic==TRUE) {
           if(dic.list[i]<y.lim[1]) {
             y.lim[1] <- y.lim[1]-abs(y.lim[1]-dic.list[i])
           }
           if(dic.list[i]>y.lim[2]) {
             y.lim[2] <- y.lim[2]+abs(y.lim[2]-dic.list[i])
           }
           plot(1:max.cluster, dic.list, type="l", xlab="Number of clusters", ylab="DIC",
                main=paste("Model",i,"of",max.cluster), ylim=y.lim)
           points(i,dic.list[i],pch=19)
           mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2)
         }
    }
  }
  
     cat("\n\n ************* RUNNING FINAL MODEL ************* \n\n")
     ###Fit the final model
     if(is.null(final.cluster)){
       best.model <- which(dic.list==min(dic.list))
     }else{
       best.model <- final.cluster
       plot.dic <- FALSE
       dic.list <- NULL
       waic.list <- NULL
     }
  
     j <- n+1-best.model
     factor.clust <- cluster.prior$cluster.store[j,]
     data.temp <- data.frame(Y.real=Y.real, E.real=E.real, region=1:n, factor.clust=factor.clust)
     
     if(lincomb) {
       lc1 = inla.make.lincomb(Predictor = rep(1/(n),n))
       names(lc1) <- "intercept"
       
       M1 <- matrix(-1/n,n,n)
       diag(M1) <- (n-1)/n
       lc2 = inla.make.lincombs(Predictor = M1)
       names(lc2) <- paste("spatial.",as.character(seq(1:n)),sep="")
       
       all.lc <- c(lc1,lc2)
     }else{
       all.lc <- NULL
     }
     
     if(best.model==1){
          formula <- Y.real ~ f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
          model.final = inla(formula, family="poisson", data=data.temp, E=E.real,
                             control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                             control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE),
                             control.predictor=list(compute=TRUE, cdf=c(log(1))),
                             lincomb=all.lc,
                             control.inla=list(strategy = strategy, npoints = 21))
     }else{

       if(!extraconstr){
          formula <- Y.real ~ factor(factor.clust) + f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                                                           hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
       }else{
         k <- length(unique(factor.clust))
         X <- matrix(0,n,k)
         
         aux <- as.numeric(factor(factor.clust,labels=1:k))
         for(i in 1:n){
           X[i,aux[i]] <- 1
         }
         
        if(sum(table(factor.clust)==1)>0){
          		Area.constr <- matrix(t(X)[table(factor.clust)!=1,], dim(X)[2]-sum(table(factor.clust)==1), n)
	    	}else{
          		Area.constr <- t(X)
	    	}
            
         formula <- Y.real ~ factor(factor.clust) +
                             f(region, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                               extraconstr=list(A=Area.constr, e=rep(0,dim(Area.constr)[1])),
                               hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
         }

	   ## Compute posterior distributions for the fixed effects ##
	   beta <- names(table(factor.clust))[-1]
	   for(i in beta){
		beta.lc <- list(list(list(list(weight=1)),list(list(weight=1))))
		names(beta.lc) <- paste("beta",i,sep="")
		names(beta.lc[[1]][[1]]) <- "(Intercept)"
		names(beta.lc[[1]][[2]]) <- paste("factor(factor.clust)",i,sep="")
		
		all.lc <- c(all.lc,beta.lc)
	   }

         model.final <- inla(formula, family="poisson", data=data.temp, E=E.real,
                             control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                             control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE),
                             control.predictor=list(compute=TRUE, cdf=c(log(1))),
                             lincomb=all.lc,
                             control.inla=list(strategy = strategy, npoints = 21))
         }
     
     if(plot.dic==TRUE) {
       plot(1:max.cluster, dic.list, type="l", xlab="Number of clusters", ylab="DIC",
            main=paste("Model",i,"of",max.cluster), ylim=y.lim)
       points(which.min(dic.list),min(dic.list), pch=19, col="red")
       lines(rep(which.min(dic.list),2),c(0,min(dic.list)), lty=2, col="red")
       axis(1, at=which.min(dic.list), labels=which.min(dic.list), col.axis="red")
       mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2, col="red")
     }
})

     results <- list(model.final=model.final, factor.clust=factor.clust, dic.list=dic.list, waic.list=waic.list, cpu.time=time[3])
     return(results)

}