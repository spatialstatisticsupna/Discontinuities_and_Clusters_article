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


#############################################
#### Two-stage SPATIAL cluster modelling ####
#############################################
library(Hmisc)
library(INLA)
library(RColorBrewer)
library(sf)
library(spdep)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


##############################################
##  1) Load the original data:              ##
##	- Spanish stomach cancer mortality data ##
##############################################
load("../data/StomachCancer_ESP.Rdata")
str(Data)

## Map of the standardized mortality ratio (SMR) for the year 2013 ##
Data.2013 <- Data[Data$year==2013,]
carto.SMR <- merge(Carto.ESP, Data.2013)

paleta <- brewer.pal(6,"YlOrRd")
values <- c(0.55,0.67,0.82,1,1.22,1.49,1.82)

Map.SMR <- tm_shape(carto.SMR) +
  tm_polygons(col="SMR", palette=paleta, legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title="Stomach cancer mortality data in Spanish provinces during the year 2013",
            main.title.position="left", main.title.size=0.8,
            legend.outside=T, legend.outside.position="right", legend.frame=F)
print(Map.SMR)


##########################################
## 2) Fit spatial models with INLA      ##
## - LCAR model (no clusters)           ##
## - Model 1 (fixed + random effects)   ##
## - Model 2 (two-level random effects) ##
##########################################
n <- length(unique(Data$prov)) ## n=47 provinces

Rs <- diag(apply(W,2,sum))-W  ## Spatial neighborhood matrix


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


################
## LCAR model ##
################
R.Leroux <- diag(n)-Rs
Data.INLA <- data.frame(O=Data.2013$obs, E=Data.2013$exp, ID.area=1:n)

formula <- O ~ f(ID.area, model="generic1", Cmatrix=R.Leroux, constr=TRUE,
                 hyper=list(prec=list(prior=sdunif), beta=list(prior=lunif,initial=0)))

LCAR.model <- inla(formula, family="poisson", data=Data.INLA, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy="laplace"))

summary(LCAR.model)


#############
## Model 1 ##
#############
source("AHC_function.R") ## Agglomerative Hierarchical Clustering algorithm (Anderson et al, 2014)

## Stage 1: applied the AHC algorithm for data corresponding to years 2010-2012 ##
SMR.prior <- cbind(log(Data[Data$year==2012,"SMR"]),
                   log(Data[Data$year==2011,"SMR"]),
                   log(Data[Data$year==2010,"SMR"]))

cluster.conf <- clustering.function(SMR.prior, W)


## Stage 2: Select the best model according the Deviance Information Criteria (DIC) ##
source("Model1.R")

Model1 <- model.selection(cluster.prior=cluster.conf, C=Rs, lincomb=TRUE,
                          plot.dic=TRUE, Y.real=Data.2013$obs, E.real=Data.2013$exp)

str(Model1,1)
summary(Model1$model.final)


## Significant clusters: 
## We must compute the posterior probability distributions
##  Pr(beta_1*>0)=Pr(alpha>alpha*)
##	Pr(beta_k*>0)=Pr(beta+alpha>alpha*)
##
## where
##  alpha*  = mean(log(r_i))
##  beta_1* = alpha-alpha*
##  beta_k* = beta_k+beta_1*, for k=2,...,K

alpha.star <- Model1$model.final$summary.lincomb.derived[1,"mean"]
beta1.star <- Model1$model.final$summary.fixed['(Intercept)',"mean"]+alpha.star
beta.prob <- data.frame(ID=unique(Model1$factor.clust),
                        mean=c(beta1.star,Model1$model.final$summary.fixed[-1,"mean"]+beta1.star),
                        prob=c(1-inla.pmarginal(alpha.star,Model1$model.final$marginals.fixed$'(Intercept)'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta2'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta3'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta5'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta6'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta7'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta8'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta14'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta15'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta25'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta27'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta28'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta39'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta40'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta41'),
                               1-inla.pmarginal(alpha.star,Model1$model.final$marginals.lincomb.derived$'beta46')))

Carto.ESP$ID <- Model1$factor.clust
cluster.map <- aggregate(Carto.ESP[,"geometry"], by=list(ID=Carto.ESP$ID), head)
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


#############
## Model 2 ##
#############
source("AHC_function.R") ## Agglomerative Hierarchical Clustering algorithm (Anderson et al, 2014)

## Stage 1: applied the AHC algorithm for data corresponding to years 2010-2012 ##
SMR.prior <- cbind(log(Data[Data$year==2012,"SMR"]),
                   log(Data[Data$year==2011,"SMR"]),
                   log(Data[Data$year==2010,"SMR"]))

cluster.conf <- clustering.function(SMR.prior, W)


## Stage 2: Select the best model according the Deviance Information Criteria (DIC) ##
source("Model2.R")

Model2 <- model.selection(cluster.prior=cluster.conf, C=Rs, lincomb=TRUE, Carto=Carto.ESP,
                          plot.dic=TRUE, Y.real=Data.2013$obs, E.real=Data.2013$exp)

str(Model2,1)
summary(Model2$model.final)


## Significant clusters: 
cl.prob <- data.frame(ID=rownames(Model2$model.final$summary.random$clust),
                      prob=unlist(lapply(Model2$model.final$marginals.random$clust, function(x) 1-inla.pmarginal(0,x))))

Carto.ESP$ID <- Model2$factor.clust
cluster.map <- aggregate(Carto.ESP[,"geometry"], by=list(ID=Carto.ESP$ID), head)
cluster.map$prob <- unlist(lapply(Model2$model.final$marginals.random$clust, function(x) 1-inla.pmarginal(0,x)))

paleta <- brewer.pal(9,"RdYlGn")[c(9,8,7,5,3,2,1)]
values <- c(0,0.05,0.1,0.2,0.8,0.9,0.95,1)

Map.cluster <- tm_shape(cluster.map) +
  tm_polygons(col="prob", palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title=expression(Pr(psi[j(i)]>0 ~"|"~ O)),
            main.title.position=0.3, main.title.size=1,
            legend.outside=T, legend.outside.position="right", legend.frame=F)
print(Map.cluster)


###################################
## 3) Model comparison (Table 3) ##
###################################
data.frame(mean.deviance=c(LCAR.model$dic$mean.deviance,Model1$model.final$dic$mean.deviance,Model2$model.final$dic$mean.deviance),
           p.eff=c(LCAR.model$dic$p.eff,Model1$model.final$dic$p.eff,Model2$model.final$dic$p.eff),
           DIC=c(LCAR.model$dic$dic,Model1$model.final$dic$dic,Model2$model.final$dic$dic),
           WAIC=c(LCAR.model$waic$waic,Model1$model.final$waic$waic,Model2$model.final$waic$waic),
           LS=c(-sum(log(LCAR.model$cpo$cpo)),-sum(log(Model1$model.final$cpo$cpo)),-sum(log(Model2$model.final$cpo$cpo))),
           row.names=c("LCAR model", "Model 1", "Model 2"),
           signif.clusters=c(NA,sum(beta.prob$prob>0.95 | beta.prob$prob<0.05),sum(cl.prob$prob>0.95 | cl.prob$prob<0.05)))
