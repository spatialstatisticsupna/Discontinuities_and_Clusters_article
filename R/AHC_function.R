########Clustering Function################
clustering.function <- function(data, W){

  n.prior <- ncol(data)

  ###Creating the storage matrices###
  cluster.store <- matrix(rep(0,nrow(data)*nrow(data)), ncol=nrow(data))
  cluster.store[1,] <- 1:nrow(data)
  cluster.size <- rep(1,nrow(data))
  dissim.size <- rep(0,nrow(data)-1)
  merged <- matrix(rep(0,2*(nrow(data)-1)),ncol=2)

  merge.ind <- -cluster.store[1,]
  initial.W <- W


  ############Centroid Linkage###############
  update.data <- data



  ###Loop##

  for (i in 1:(nrow(data)-1)){

    ###Count cluster size###
    for (j in 1:nrow(data)){
      cluster.size[j] <- sum(cluster.store[i,]==j)
    }


    ###Induce cluster structure in dissimilarity matrix###
    clusts <- cluster.store[i,][cluster.size>1]


    if(length(clusts)>0){
      for (j in 1:length(clusts)){
        num.row <- sum(cluster.store[i,]==clusts[j])
        update.data[cluster.store[i,]==clusts[j],] <- matrix(rep(apply(as.matrix(data[cluster.store[i,]==clusts[j],]),2,mean),num.row),nrow=num.row, byrow=T)
      }
    }

    
    ## Compute the distance matrix ##
    sim.mat <- as.matrix(dist(update.data, diag=TRUE, upper=TRUE))

    sim.mat[lower.tri(sim.mat)] <- max(sim.mat)+1
    sim.mat[W==0] <- max(sim.mat)
    diag(sim.mat) <- max(sim.mat)

    
    ###Ensure points in same cluster aren't compared###
    clust.mat <- as.matrix(dist(cbind(cluster.store[i,], cluster.store[i,]), method="maximum", diag=TRUE, upper=TRUE))
    sim.mat[clust.mat==0] <- max(sim.mat)

    ###Find most similar###
    similar <- which(sim.mat==min(sim.mat), arr.ind=TRUE)
    choice <- sample(x=1:nrow(similar), size=1)
    join.1 <- cluster.store[i,similar[choice,1]]
    join.2 <- cluster.store[i,similar[choice,2]]


    ###Store results###
    cluster.store[(i+1),] <- cluster.store[i,]
    cluster.store[(i+1),][cluster.store[(i+1),]==join.1] <- min(join.1,join.2)
    cluster.store[(i+1),][cluster.store[(i+1),]==join.2] <- min(join.1,join.2)

    dissim.size[i] <- min(sim.mat)
    merged[i,] <- c(merge.ind[similar[choice,1]], merge.ind[similar[choice,2]])
    old.merge.ind <- merged[i,]
    merge.ind[merge.ind==old.merge.ind[1]] <- i
    merge.ind[merge.ind==old.merge.ind[2]] <- i



    ###Update the W matrix ###
    new.W.row <- apply((W[cluster.store[(i+1),]==min(join.1,join.2),]),2,sum)
    num.replaced <- nrow(W[cluster.store[(i+1),]==min(join.1,join.2),])
    W[cluster.store[(i+1),]==min(join.1,join.2),] <- matrix(rep(new.W.row,num.replaced),nrow=num.replaced, byrow=TRUE)
    W[,cluster.store[(i+1),]==min(join.1,join.2)] <- matrix(rep(new.W.row,num.replaced),ncol=num.replaced, byrow=FALSE)

    W[W>0] <- 1
  }



  W.list <- vector("list", nrow(data))
  for (i in 1:nrow(data)){
    W.list[[i]] <- matrix(data = NA, nrow = nrow(data), ncol = nrow(data),
                          byrow= FALSE, dimnames = NULL)
  }


  for(i in 1:nrow(data)){
    clust.mat <- as.matrix(dist(cbind(cluster.store[i,], cluster.store[i,]), method="maximum", diag=TRUE, upper=TRUE))
    W.list[[i]] <- initial.W
    W.list[[i]][clust.mat>0] <- 0
  }

  result <- list(W.list=W.list, cluster.store=cluster.store, merge.ind=merge.ind, dissim.size=dissim.size, merged=merged)
  return(result)
}