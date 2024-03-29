source("Cluster_neighbours.R")
source("kronecker_nullspace.R")

model.selection <- function(cluster.prior, Y.real, E.real, Qs, Qt, Carto=NULL,
				    max.cluster=NULL, final.cluster=NULL, plot.dic=TRUE,
				    strategy="gaussian", lincomb=FALSE){
                            
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
      logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*beta+(b-1)*(1-beta);
      log_jacobian = log(beta*(1-beta));
      return(logdens+log_jacobian)"
    
    R.area.Leroux <- diag(dim(Qs)[1])-Qs
  
    
    if(is.null(final.cluster)){
  
      ## Store the DIC values
      n <- dim(Qs)[1]
      t <- dim(Qt)[1]
      if(is.null(max.cluster)) {max.cluster <- n-1}
      dic.list <- rep(NA, max.cluster)
      waic.list <- rep(NA, max.cluster) 
      
        
      ## Fit a model with a single cluster
      i <- 1
      
      data.temp <- data.frame(Y.real=Y.real, E.real=E.real, ID.area=rep(1:n,t),
                              ID.year=rep(1:t,each=n), ID.area.year=seq(1,n*t))
  
      R <- kronecker(Qt,Qs)
      r.def <- n+t-1
      A1 <- kronecker(matrix(1,1,t),diag(n))
      A2 <- kronecker(diag(t),matrix(1,1,n))
      A_delta <- rbind(A1[-1,],A2[-1,])
  
      formula <- Y.real ~ f(ID.area, model="generic1", Cmatrix=R.area.Leroux, constr=TRUE,
                            hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
                      	  f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
                          f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE,
                            extraconstr=list(A=A_delta, e=rep(0,dim(A_delta)[1])),
                            hyper=list(prec=list(prior=sdunif)))
  
      model <- inla(formula, family="poisson", data=data.temp, E=E.real,
                    control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy))
  
      dic.list[i] <- model$dic$dic
      waic.list[i] <- model$waic$waic
      
      if(plot.dic==TRUE) {
        win.graph()
        y.lim <- dic.list[i]*c(1/1.01,1.005)
        plot(1:max.cluster, dic.list, type="l", xlab="Number of spatial clusters", ylab="DIC",
             main=paste("Model",i,"of",max.cluster), ylim=y.lim)
        points(1,dic.list[i],pch=19)
        mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2)
      }
     

      ## Fit separate models with between 2 and the max number of clusters.
      for(i in 2:max.cluster){
      
	      j <- n-i+1
	      factor.clust <- as.numeric(as.factor(cluster.prior$cluster.store[j,]))
	      m <- length(unique(factor.clust))
      
	      data.temp <- data.frame(Y.real=Y.real, E.real=E.real, ID.area=rep(1:n,t),
	                              ID.clust=rep(factor.clust,t), ID.year=rep(1:t,each=n),
	                              ID.clust.year=factor.clust+rep(0:(t-1)*m,each=n))
      
	      ## Cluster configuration map ##
	      cl.nb <- Cluster.neighbours(Carto,factor.clust)
	      R.area.Leroux <- diag(n)-cl.nb$C
      
	      Carto$factor.clust <- factor.clust
        cluster.map <- aggregate(Carto[,"geometry"], by=list(factor.clust=Carto$factor.clust), head)
	      nb2INLA(poly2nb(cluster.map), file="cluster_nb.inla")
                
	      g <- inla.read.graph("cluster_nb.inla")
	      R.clust = matrix(0, g$n, g$n)
	      for (ii in 1:g$n){
	      	R.clust[ii,ii]=g$nnbs[[ii]]
                  R.clust[ii,g$nbs[[ii]]]=-1
	      }
	      R.clust.Leroux <- diag(dim(R.clust)[1])-R.clust
      
	      ## INLA formula ##
	      if(length(cl.nb$single.cluster)>0){
                	Area.constr <- matrix(cl.nb$C.constr[-cl.nb$single.cluster,], dim(cl.nb$C.constr)[1]-length(cl.nb$single.cluster), n)
	      }else{
                	Area.constr <- cl.nb$C.constr
	      }
	      
	      R <- kronecker(Qt,R.clust)
	      n.cluster <- dim(R.clust)[1]
	      r.def <- n.cluster+t-1
	      A1 <- kronecker(matrix(1,1,t),diag(n.cluster))
	      A2 <- kronecker(diag(t),matrix(1,1,n.cluster))
	      A_delta <- rbind(A1[-1,],A2[-1,])
      
	      formula <- Y.real ~ f(ID.area, model="generic1", Cmatrix=R.area.Leroux, constr=ifelse(length(cl.nb$single.cluster)==0,FALSE,TRUE), 
	                            extraconstr=list(A=Area.constr, e=rep(0,dim(Area.constr)[1])),
	                            hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
	      	   		            f(ID.clust, model="generic1", Cmatrix=R.clust.Leroux, constr=TRUE, 
	      	   		              hyper=list(beta=list(prior=lunif, initial=0))) + 
	      			              f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
	      			              f(ID.clust.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE,
	      			                extraconstr=list(A=A_delta, e=rep(0,dim(A_delta)[1])),
	      			                hyper=list(prec=list(prior=sdunif)))
      
	      model <- inla(formula, family="poisson", data=data.temp, E=E.real,
	      		  control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
	      		  control.predictor=list(compute=TRUE, cdf=c(log(1))),
	      		  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
	      		  control.inla=list(strategy=strategy))
      
	      dic.list[i] <- model$dic$dic
	      waic.list[i] <- model$waic$waic
                
	      if(plot.dic==TRUE) {
	      	if(dic.list[i]<y.lim[1]) {
                    y.lim[1] <- y.lim[1]-abs(y.lim[1]-dic.list[i])
                  }
	      	if(dic.list[i]>y.lim[2]) {
	      	  y.lim[2] <- y.lim[2]+abs(y.lim[2]-dic.list[i])
	      	}
	      	plot(1:max.cluster, dic.list, type="l", xlab="Number of spatial clusters", ylab="DIC",
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
	factor.clust <- as.numeric(as.factor(cluster.prior$cluster.store[j,]))
	m <- length(unique(factor.clust))

	data.temp <- data.frame(Y.real=Y.real, E.real=E.real, ID.area=rep(1:n,t),
	                        ID.clust=rep(factor.clust,t), ID.year=rep(1:t,each=n),
	                        ID.clust.year=factor.clust+rep(0:(t-1)*m,each=n),
	                        ID.area.year=seq(1:(n*t)))

	if(lincomb) {
		lc1 = inla.make.lincomb(Predictor = rep(1/(n*t),n*t))
		names(lc1) <- "intercept"

		M1 <- matrix(-1/(n*t),n,n*t)
		for (i in 1:n) {
			M1[i,seq(0,n*(t-1),by=n)+i] <- 1/(n*t)*(n-1)
		}
		lc2 = inla.make.lincombs(Predictor = M1)
		names(lc2) <- paste("spatial.",as.character(seq(1:n)),sep="")

		M2 <- matrix(-1/(n*t),t,n*t)
		for (i in 1:t) {
			M2[i,seq(1,n)+n*(i-1)] <- 1/(n*t)*(t-1)
		}
		lc3 = inla.make.lincombs(Predictor = M2)
		names(lc3) <- paste("temporal.",as.character(seq(1:t)),sep="")

		M3 <- matrix(1/(n*t),n*t,n*t)
		k=1
		for (j in 1:t) {
			for (i in 1:n) {
				M3[k,seq(0,n*(t-1),by=n)+i] <- 1/(n*t)*(1-n)
				M3[k,seq(1,n)+n*(j-1)] <- 1/(n*t)*(1-t)
				M3[k,k] <- (n*t-n-t+1)/(n*t)
				k=k+1
			}
		}
		lc4 = inla.make.lincombs(Predictor = M3)
		names(lc4) <- paste("spatio.temporal.",as.character(seq(1:(n*t))),sep="")

		all.lc <- c(lc1,lc2,lc3,lc4)
  }else{
     all.lc <- NULL
  }
     
  if(best.model==1){
    
    R <- kronecker(Qt,Qs)
    r.def <- n+t-1
    A1 <- kronecker(matrix(1,1,t),diag(n))
    A2 <- kronecker(diag(t),matrix(1,1,n))
    A_delta <- rbind(A1[-1,],A2[-1,])

    formula <- Y.real ~ f(ID.area, model="generic1", Cmatrix=R.area.Leroux, constr=TRUE,
                          hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
                     	 	f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE,
                          extraconstr=list(A=A_delta, e=rep(0,dim(A_delta)[1])),
                          hyper=list(prec=list(prior=sdunif)))

    model.final <- inla(formula, family="poisson", data=data.temp, E=E.real,
                        control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
                        control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   	  	control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   	  	control.inla=list(strategy=strategy))
  }else{

		## Cluster configuration map ##
		cl.nb <- Cluster.neighbours(Carto,factor.clust)
		R.area.Leroux <- diag(n)-cl.nb$C

		Carto$factor.clust <- factor.clust
		cluster.map <- aggregate(Carto[,"geometry"], by=list(factor.clust=Carto$factor.clust), head)
		nb2INLA(poly2nb(cluster.map), file="cluster_nb.inla")
          
		g <- inla.read.graph("cluster_nb.inla")
		R.clust = matrix(0, g$n, g$n)
		for (ii in 1:g$n){
			R.clust[ii,ii]=g$nnbs[[ii]]
      	      R.clust[ii,g$nbs[[ii]]]=-1
		}
		R.clust.Leroux <- diag(dim(R.clust)[1])-R.clust

		## INLA formula ##
		if(length(cl.nb$single.cluster)>0){
      	    	Area.constr <- matrix(cl.nb$C.constr[-cl.nb$single.cluster,], dim(cl.nb$C.constr)[1]-length(cl.nb$single.cluster), n)
		}else{
      	    	Area.constr <- cl.nb$C.constr
		}

		R <- kronecker(Qt,R.clust)
		n.cluster <- dim(R.clust)[1]
		r.def <- n.cluster+t-1
		A1 <- kronecker(matrix(1,1,t),diag(n.cluster))
		A2 <- kronecker(diag(t),matrix(1,1,n.cluster))
		A_delta <- rbind(A1[-1,],A2[-1,])

		formula <- Y.real ~ f(ID.area, model="generic1", Cmatrix=R.area.Leroux, constr=ifelse(length(cl.nb$single.cluster)==0,FALSE,TRUE), 
		                      extraconstr=list(A=Area.constr, e=rep(0,dim(Area.constr)[1])),
		                      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
		   		  	          f(ID.clust, model="generic1", Cmatrix=R.clust.Leroux, constr=TRUE, 
		   		  	            hyper=list(beta=list(prior=lunif, initial=0))) + 
				  	            f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
				  	            f(ID.clust.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE,
				  	              extraconstr=list(A=A_delta, e=rep(0,dim(A_delta)[1])),
				  	              hyper=list(prec=list(prior=sdunif)))

		model.final <- inla(formula, family="poisson", data=data.temp, E=E.real,
		                    control.fixed=list(mean=0, mean.intercept=0, prec=0.1, prec.intercept=0.001),
		                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   		  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   		  control.inla=list(strategy=strategy),
		                    lincomb=all.lc)
  }

  if(plot.dic==TRUE) {
    plot(1:max.cluster, dic.list, type="l", xlab="Number of spatial clusters", ylab="DIC",
         main=paste("Model",best.model,"of",max.cluster), ylim=y.lim)
    points(which.min(dic.list),min(dic.list), pch=19, col="red")
    lines(rep(which.min(dic.list),2),c(0,min(dic.list)), lty=2, col="red")
    axis(1, at=which.min(dic.list), labels=which.min(dic.list), col.axis="red")
    mtext(paste("Minimun DIC value:",round(min(dic.list,na.rm=T),2)), side=3, line=-2, col="red")
  }
})

   results <- list(model.final=model.final, factor.clust=factor.clust, dic.list=dic.list, waic.list=waic.list, cpu.time=time[3])
   return(results)
}