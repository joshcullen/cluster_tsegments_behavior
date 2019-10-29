# rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(1)

#get functions
sourceCpp('aux1.cpp')
source('gibbs functions.R')
source('gibbs sampler.R')

#get data
dat=read.csv('fake data.csv',as.is=T)


#run gibbs sampler
res=cluster.tsegm.behavior(dat=dat,a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,gamma1=gamma1,
                           nclustmax=nclustmax,ngibbs=ngibbs,nburn=nburn)

