######################################################################################
##  The function "kronecker.null.space" returns the Kronecker product of 2 matrix   ##
##  and the eigenvectors which span the null space of this Kronecker product		##
######################################################################################
kronecker.null.space <- function(Matrix1, Matrix2){

M1.sparse <- as(Matrix1, "dgTMatrix")
M2.sparse <- as(Matrix2, "dgTMatrix")

n1 <- dim(M1.sparse)[1]
n2 <- dim(M2.sparse)[1]


## Null-space of first matrix ##
eigen.M1 <- eigen(M1.sparse)
if (sum(eigen.M1$values<1e-8)>0) {
	null.space.M1 <- as.matrix(eigen.M1$vectors[,eigen.M1$values<1e-8])
	mult.M1 <- dim(null.space.M1)[2]
} else {
	null.space.M1 <- NULL
	mult.M1 <- 0
}

## Null-space of second matrix ##
eigen.M2 <- eigen(M2.sparse)
if (sum(eigen.M2$values<1e-8)>0) {
	null.space.M2 <- as.matrix(eigen.M2$vectors[,eigen.M2$values<1e-8])
	mult.M2 <- dim(null.space.M2)[2]
} else {
	null.space.M2 <- NULL
	mult.M2 <- 0
}


## Null-space of Kronecker product ##
if (mult.M1>0) {
	kronecker.null.space.1 <-  as(null.space.M1[,1] %x% eigen.M2$vectors, "dgTMatrix")
	if (mult.M1>1) {
		for (i in 2:mult.M1) {
			kronecker.null.space.1 <- cBind(kronecker.null.space.1, null.space.M1[,i] %x% eigen.M2$vectors)
		}
	}

	if (mult.M2>0) {
		kronecker.null.space.2 <-  as(eigen.M1$vectors[,-((n1-mult.M1+1):n1)] %x% null.space.M2[,1], "dgTMatrix")
		if (mult.M2>1) {
			for (i in 2:mult.M2) {
				kronecker.null.space.2 <- cBind(kronecker.null.space.2, eigen.M1$vectors[,-((n1-mult.M1+1):n1)] %x% null.space.M2[,i])
			}
		}
	} else {
		kronecker.null.space.2 <- NULL
	}

} else {
	kronecker.null.space.1 <- NULL

	if (mult.M2>0) {
		kronecker.null.space.2 <-  as(eigen.M1$vectors %x% null.space.M2[,1], "dgTMatrix")
		if (mult.M2>1) {
			for (i in 2:mult.M2) {
				kronecker.null.space.2 <- cBind(kronecker.null.space.2, eigen.M1$vectors %x% null.space.M2[,i])
			}
		}
	} else {
		kronecker.null.space.2 <- NULL
	}
}

if(!is.null(kronecker.null.space.1) & !is.null(kronecker.null.space.2)) {
	kronecker.null.space <- cBind(kronecker.null.space.1,kronecker.null.space.2)
} else if(!is.null(kronecker.null.space.1)) {
	kronecker.null.space <- kronecker.null.space.1
} else if(!is.null(kronecker.null.space.2)) {
	kronecker.null.space <- kronecker.null.space.2
} else {
	kronecker.null.space <- NULL
}

## Return Kroncker product matrix and the corresponding Null-space ##
result <- vector("list",2)
	
result[[1]] <- M1.sparse %x% M2.sparse
if(!is.null(kronecker.null.space))	result[[2]] <- t(kronecker.null.space)
names(result) <- c("kronecker","null.space")

result
}



##########################################################################
##  The function "check.kronecker.null.space" checks that the results   ##
##  given by "kronecker.null.space" function are correct  			##
##########################################################################
check.kronecker.null.space <- function(result) {

## Check if the input object is the a list with the structure of the "kronecker.null.space" function result ##
if(class(result[[1]])!="dgTMatrix")
	stop("The first element of the input object is not of class dgTMatrix")

if(is.null(result[[2]]))
	stop("This Kronecker product matrix null-space has dimension 0")

for (i in 1:dim(result[[2]])[1]) {
	aux <- as.matrix(result[[1]])%*%t(result[[2]])[,i]
	if(length(aux[aux<1e-8])!=dim(result[[2]])[2])
		stop(paste("The ",i,".th vector is not a eigenvector of the Kronekcer product matrix null-space",sep=""))
}

print("All the vectors are eigenvectors of the Kronecker product matrix null-space")
}


###################################################################################################
##  The function "triple.kronecker.null.space" returns the Kronecker product of 3 matrix   	 ##
##  and the eigenvectors which span the null space of this triple Kronecker product		 	 ##
###################################################################################################
triple.kronecker.null.space <- function(Matrix1, Matrix2, Matrix3){

## Call "kronecker.null.space" function ##
aux <- kronecker.null.space(Matrix1, Matrix2)

M1.sparse <- aux[[1]]
M2.sparse <- as(Matrix3, "dgTMatrix")

n1 <- dim(M1.sparse)[1]
n2 <- dim(M2.sparse)[1]

## Null-space of first matrix ##
eigen.M1 <- eigen(M1.sparse)
if(!is.null(aux[[2]])){
	null.space.M1 <- t(aux[[2]])
	mult.M1 <- dim(null.space.M1)[2]
} else {
	null.space.M1 <- NULL
	mult.M1 <- 0
}

## Null-space of second matrix ##
eigen.M2 <- eigen(M2.sparse)
if (sum(eigen.M2$values<1e-8)>0) {
	null.space.M2 <- as.matrix(eigen.M2$vectors[,eigen.M2$values<1e-8])
	mult.M2 <- dim(null.space.M2)[2]
} else {
	null.space.M2 <- NULL
	mult.M2 <- 0
}


## Null-space of Kronecker product ##
if (mult.M1>0) {
	kronecker.null.space.1 <-  as(null.space.M1[,1] %x% eigen.M2$vectors, "dgTMatrix")
	if (mult.M1>1) {
		for (i in 2:mult.M1) {
			kronecker.null.space.1 <- cBind(kronecker.null.space.1, null.space.M1[,i] %x% eigen.M2$vectors)
		}
	}

	if (mult.M2>0) {
		kronecker.null.space.2 <-  as(eigen.M1$vectors[,-((n1-mult.M1+1):n1)] %x% null.space.M2[,1], "dgTMatrix")
		if (mult.M2>1) {
			for (i in 2:mult.M2) {
				kronecker.null.space.2 <- cBind(kronecker.null.space.2, eigen.M1$vectors[,-((n1-mult.M1+1):n1)] %x% null.space.M2[,i])
			}
		}
	} else {
		kronecker.null.space.2 <- NULL
	}

} else {
	kronecker.null.space.1 <- NULL

	if (mult.M2>0) {
		kronecker.null.space.2 <-  as(eigen.M1$vectors %x% null.space.M2[,1], "dgTMatrix")
		if (mult.M2>1) {
			for (i in 2:mult.M2) {
				kronecker.null.space.2 <- cBind(kronecker.null.space.2, eigen.M1$vectors %x% null.space.M2[,i])
			}
		}
	} else {
		kronecker.null.space.2 <- NULL
	}
}

if(!is.null(kronecker.null.space.1) & !is.null(kronecker.null.space.2)) {
	kronecker.null.space <- cBind(kronecker.null.space.1,kronecker.null.space.2)
} else if(!is.null(kronecker.null.space.1)) {
	kronecker.null.space <- kronecker.null.space.1
} else if(!is.null(kronecker.null.space.2)) {
	kronecker.null.space <- kronecker.null.space.2
} else {
	kronecker.null.space <- NULL
}

## Return Kroncker product matrix and the corresponding Null-space ##
result <- vector("list",2)
	
result[[1]] <- M1.sparse %x% M2.sparse
if(!is.null(kronecker.null.space))	result[[2]] <- t(kronecker.null.space)
names(result) <- c("kronecker","null.space")

result
}





########################################################
########################################################
## 				EXAMPLES:				##
##									##
## example1 <- kronecker.null.space(Qs,Qt) 		##
## check.kronecker.null.space(example1)			##
## R <- example1[[1]]						##
## A_delta <- as.matrix(example1[[2]])			##
##									##
########################################################
##									##
## example2 <- triple.kronecker.null.space(Qt,Is,It)	##
## R <- example2[[1]]									##
## A_delta <- as.matrix(example2[[2]])								##
## 									##
########################################################

