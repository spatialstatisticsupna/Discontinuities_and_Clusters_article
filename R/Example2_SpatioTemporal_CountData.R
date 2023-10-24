##########################################################
####	INLA code to fit the models described in:       ####
####	"A two-stage approach to estimate spatial and   ####
####   spatio-temporal disease risks in the presence  ####
####   of local discontinuities and clusters"         ####
####	(Adin et al., 2019)                             ####
##########################################################
rm(list=ls())

# install.packages("spdep", dependencies=TRUE)
# install.packages("sp", dependencies=TRUE)
# install.packages("Hmisc", dependencies=TRUE)
# install.packages("RColorBrewer", dependencies=TRUE)
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)


#####################################################
#### Two-stage SPATIO-TEMPORAL cluster modelling ####
#####################################################
library(Hmisc)
library(INLA)
library(RColorBrewer)
library(sf)
library(spdep)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


##############################################
##  1) Load the original data:              ##
##	- Brain cancer incidence data in the    ##
##    regions of Navarre and Basque Country ##
##    during the period 2000-2008           ##
##############################################
load("../data/BrainCancer_MUN.Rdata")

str(Data.INCI)
str(Data.MORT)


## Map of the standardized incidence ratio (SIR) for the period 2000-2008 ##
n <- length(unique(Data.INCI$region)) # n <- 27
t <- length(unique(Data.INCI$year))   # t <- 9
t.from <- min(unique(Data.INCI$year)) # t.from <- 2000
t.to <- max(unique(Data.INCI$year))   # t.to <- 2008

SMR.matrix <- matrix(Data.INCI$SMR,n,t,byrow=F)
SMR.data.frame <- data.frame(Carto.COM$region,SMR.matrix)
colnames(SMR.data.frame) <- c("region",paste("SMR",seq(t.from,t.to),sep=""))

carto.SMR <- merge(Carto.COM, SMR.data.frame)

paleta <- brewer.pal(8,"YlOrRd")
values <- c(0,0.74,0.82,0.90,1,1.11,1.22,1.35,4.5)

Map.SMR <- tm_shape(carto.SMR) +
  tm_polygons(col=paste("SMR",seq(t.from,t.to),sep=""),
              palette=paleta, legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title="Brain cancer incidence data", main.title.position=0.2, panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=paste("Year",seq(t.from,t.to)),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0))
print(Map.SMR)


############################################################################
## 2) Fit spatio-temporal models with INLA                                ## 
## + Model with no clusters:                                              ##
##   - Type IV interaction LCAR model                                     ##
## + Models with purely spatial clusters                                  ##
##   - TL-Model 1a (fixed + random effects / area-level interaction)      ##
##   - TL-Model 1b (fixed + random effects / cluster-level interaction)   ##
##   - TL-Model 2a (two-level random effects / area-level interaction)    ##
##   - TL-Model 2b (two-level random effects / cluster-level interaction) ##
## + Models with spatio-tempora clusters                                  ##
##   - Option I (fixed + random effects)                                  ##
##   - Option II, III and IV (two-level random effects)                   ##
############################################################################

Rs <- as(diag(apply(W,2,sum))-W,"Matrix")  ## Spatial neighborhood matrix

D1 <- diff(diag(t),differences=1) ## RW1 structure matrix
Rt <- as(t(D1)%*%D1,"Matrix")

## Hyperprior distributions using the "expression()" function
## - Unif(0,Inf) for standard deviations
## - Unif(0,1) for the spatial smoothing parameter

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


####################################
## Type IV interaction LCAR model ##
####################################
R.Leroux <- Diagonal(n)-Rs

R <- kronecker(Rt,Rs) ## Interaction structure matrix

## Identifiability constraints (Goicoa et al., 2017)
r.def <- n+t-1
A1 <- kronecker(matrix(1,1,t),diag(n))
A2 <- kronecker(diag(t),matrix(1,1,n))
A.constr <- rbind(A1[-1,],A2[-1,])

Data.INLA <- data.frame(O=Data.INCI$obs, E=Data.INCI$exp,
                        ID.area=rep(1:n,t), ID.year=rep(1:t,each=n),
                        ID.area.year=seq(1,n*t))

formula <- O ~ f(ID.area, model="generic1", Cmatrix = R.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
               f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
               f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE,
                 extraconstr=list(A=A.constr, e=rep(0,n+t-2)),
                 hyper=list(prec=list(prior=sdunif)))

LCAR.model <- inla(formula, family="poisson", data=Data.INLA, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy="simplified.laplace"))

summary(LCAR.model)


###############################
## TL-Model 1a / TL-Model 1b ##
###############################
source("AHC_function.R") ## Agglomerative Hierarchical Clustering algorithm (Anderson et al, 2014)

## Stage 1: applied the AHC algorithm to mortality data (aggregated by regions) ##
SMR.prior <- aggregate(Data.MORT[,c("obs","exp")], by=list(region=Data.MORT$region), sum)
SMR.prior$SMR <- SMR.prior$obs/SMR.prior$exp

cluster.conf <- clustering.function(data=matrix(log(SMR.prior$SMR),n,1), W)


## Stage 2: Select the best model according the Deviance Information Criteria (DIC) ##
source("TLmodel1a.R")
# source("TLmodel1b.R")

TLModel1a <- model.selection(cluster.prior=cluster.conf, Qs=Rs, Qt=Rt, Carto=Carto.COM,
                             Y.real=Data.INCI$obs, E.real=Data.INCI$exp,
                             lincomb=TRUE,  plot.dic=TRUE, strategy="simplified.laplace")

str(TLModel1a,1)
summary(TLModel1a$model.final)


## Significant clusters
## We must compute the posterior probability distributions
##  Pr(beta_1*>0)=Pr(alpha>alpha*)
##	Pr(beta_k*>0)=Pr(beta+alpha>alpha*)
##
## where
##  alpha*  = mean(log(r_i))
##  beta_1* = alpha-alpha*
##  beta_k* = beta_k+beta_1*, for k=2,...,K
alpha.star <- TLModel1a$model.final$summary.lincomb.derived[1,"mean"]
beta.prob <- data.frame(ID=unique(TLModel1a$factor.clust),
                        prob=c(1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.fixed$'(Intercept)'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta2'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta3'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta4'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta5'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta6'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta7'),
                               1-inla.pmarginal(alpha.star,TLModel1a$model.final$marginals.lincomb.derived$'beta8')))


Carto.COM$ID <- TLModel1a$factor.clust
cluster.map <- aggregate(Carto.COM[,"geometry"], by=list(ID=Carto.COM$ID), head)
cluster.map <- merge(cluster.map,beta.prob)

paleta <- brewer.pal(9,"RdYlGn")[c(9,8,7,5,3,2,1)]
values <- c(0,0.05,0.1,0.2,0.8,0.9,0.95,1)

Map.cluster <- tm_shape(cluster.map) +
  tm_polygons(col="prob", palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title=expression(Pr(beta[j]>0 ~"|"~ O)),
            main.title.position=0.3, main.title.size=1,
            legend.outside=T, legend.outside.position="right", legend.frame=F)
print(Map.cluster)


###############################
## TL-Model 2a / TL-Model 2b ##
###############################
source("AHC_function.R") ## Agglomerative Hierarchical Clustering algorithm (Anderson et al, 2014)

## Stage 1: applied the AHC algorithm to mortality data (aggregated by regions) ##
SMR.prior <- aggregate(Data.MORT[,c("obs","exp")], by=list(region=Data.MORT$region), sum)
SMR.prior$SMR <- SMR.prior$obs/SMR.prior$exp

cluster.conf <- clustering.function(data=matrix(log(SMR.prior$SMR),n,1), W)


## Stage 2: Select the best model according the Deviance Information Criteria (DIC) ##
source("TLmodel2a.R")
#source("TLmodel2b.R")

TLModel2a <- model.selection(cluster.prior=cluster.conf, Qs=Rs, Qt=Rt, Carto=Carto.COM,
                             Y.real=Data.INCI$obs, E.real=Data.INCI$exp,
                             lincomb=TRUE, plot.dic=TRUE, strategy="simplified.laplace")

str(TLModel2a,1)
summary(TLModel2a$model.final)


## Significant clusters: 
cl.prob <- data.frame(ID=rownames(TLModel2a$model.final$summary.random$ID.clust),
                      prob=unlist(lapply(TLModel2a$model.final$marginals.random$ID.clust, function(x) 1-inla.pmarginal(0,x))))

Carto.COM$ID <- TLModel2a$factor.clust
cluster.map <- aggregate(Carto.COM[,"geometry"], by=list(ID=Carto.COM$ID), head)
cluster.map$prob <- unlist(lapply(TLModel2a$model.final$marginals.random$ID.clust, function(x) 1-inla.pmarginal(0,x)))

paleta <- brewer.pal(9,"RdYlGn")[c(9,8,7,5,3,2,1)]
values <- c(0,0.05,0.1,0.2,0.8,0.9,0.95,1)

Map.cluster <- tm_shape(cluster.map) +
  tm_polygons(col="prob", palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title=expression(Pr(psi[j(i)]>0 ~"|"~ O)),
            main.title.position=0.3, main.title.size=1,
            legend.outside=T, legend.outside.position="right", legend.frame=F)
print(Map.cluster)


###########################################
## ST clusters: Option I / II / III / IV ##
###########################################
source("AHC_clusteringST_function.R") ## Spatio-temporal AHC algorithm 

## Stage 1: applied the AHC algorithm to mortality data (in the same regions and time period) ##
SMR.prior <- cbind(Data.MORT$SMR)

cluster.conf <- clusteringST.function(data=log(SMR.prior+0.00001), W)


## Stage 2: Select the best model according the Deviance Information Criteria (DIC) ##
# source("OptionI.R")
source("OptionII.R")
# source("OptionIII.R")
# source("OptionIV.R")


STcluster.model <- model.selection(cluster.prior=cluster.conf, Qs=Rs, Qt=Rt, Carto=Carto.COM,
                                   Y.real=Data.INCI$obs, E.real=Data.INCI$exp, max.cluster=60,
                                   lincomb=TRUE, plot.dic=TRUE, strategy="simplified.laplace")

str(STcluster.model,1)
summary(STcluster.model$model.final)


## Significant clusters 
########################
delta <- unlist(lapply(STcluster.model$model.final$marginals.random$ID.clust, function(x) inla.emarginal(exp,x)))
delta.prob <- unlist(lapply(STcluster.model$model.final$marginals.random$ID.clust, function(x){1-inla.pmarginal(0, x)}))

ID <- which(delta.prob>0.95)
# ID <- which(delta.prob<0.05)

ClusterH.loc <- STcluster.model$factor.clust
ClusterH.loc[!(ClusterH.loc %in% ID)] <- NA
ClusterH.loc <- as.data.frame(matrix(ClusterH.loc,n,t,byrow=FALSE))

for(i in 1:ncol(ClusterH.loc)){
  ClusterH.loc[,i] <- factor(ClusterH.loc[,i], levels=ID, labels=paste("Cluster",length(ID):1))
}
names(ClusterH.loc) <- paste("High",t.from:t.to,sep="")
carto.clust <- cbind(Carto.COM,ClusterH.loc)

paletaH <- brewer.pal(length(ID)+3,"Reds")[c(2,4)]
paletaL <- brewer.pal(length(ID)+2,"Greens")[-c(1,4)]

Map.cluster <- tm_shape(carto.clust) +
  tm_polygons(col=paste("High",t.from:t.to,sep=""), title="",
              labels=as.character(round(sort(delta[ID])[length(ID):1],2)),
              colorNA="white", showNA=FALSE, legend.show=T, legend.reverse=T,
              palette=paletaH[order(delta[ID])],
              style="fixed", interval.closure="left") + 
  tm_layout(main.title="High risk clusters", main.title.position=0.2, panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=paste("Year",seq(t.from,t.to)),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0)) + 
  tm_facets(nrow=3, ncol=3)
print(Map.cluster)