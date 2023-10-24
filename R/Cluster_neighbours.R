Cluster.neighbours <- function(Carto, factor.clust, plot=FALSE){
  ## Function to create the neighbourhood matrix of areas within each cluster ##

  nn <- dim(Carto)[1]

  k <- length(unique(factor.clust))
  X <- matrix(0,nn,k)

  aux <- as.numeric(factor(factor.clust,labels=1:k))
  for(i in 1:nn){
    X[i,aux[i]] <- 1
  }

  single.zones <- apply(as.matrix(X[,apply(X,2,sum)==1]),2,function(x) which(x==1))
  single.cluster <- which(apply(X,2,sum)==1)

  C <- matrix(0,nn,nn)

  carto.clust <- vector("list",k)
  for (i in 1:k){
    lista <- 1:nn
    lista[which(factor.clust!=sort(unique(factor.clust))[i])] <- NA
  
    Carto$lista <- lista
    carto.aux <- aggregate(Carto[,"geometry"], by=list(lista=Carto$lista), head)
    carto.clust[[i]] <- carto.aux
  
    if(plot){
      windows()
      plot(carto.aux, main=paste("Cluster",i))
    }
  
    if(nrow(carto.aux)>1) {
      x.nb <- poly2nb(carto.aux, snap=0.01)
      nb2INLA(file=paste("Cluster_",i,".inla",sep=""), x.nb)
      if(plot) {plot(x.nb, coordinates(carto.aux), pch=19, add=TRUE)}
  
      l <- seq(1:nn)[!is.na(lista)]
      g <- inla.read.graph(paste(paste("Cluster_",i,sep=""),".inla",sep=""))
      for (j in 1:g$n){
        C[l[j],l[g$nbs[[j]]]]=-1
      }
    }else{
      if(plot) {points(coordinates(carto.aux)[1], coordinates(carto.aux)[2], pch=19)}
    }
  }

  ## Add the number of neighbour to the diagonal ##
  diag(C) <- -apply(C, 1, sum)

  ## Add a "1" in the diagonal entries corresponding to single zones ##
  for(i in 1:length(single.zones)){
    C[single.zones[i],single.zones[i]] <- 1
  }

  return(list(C=C,C.constr=t(X),single.cluster=single.cluster,single.zones=single.zones))
}