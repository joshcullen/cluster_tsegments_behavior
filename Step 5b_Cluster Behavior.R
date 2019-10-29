set.seed(1)

library('MCMCpack')
library('Rcpp')

#get functions
sourceCpp('aux1.cpp')
source('gibbs functions.R')
source('gibbs sampler2.R')
source('helper functions.R')

#get data
dat<- read.csv('Snail Kite Gridded Data_behav.csv', header = T, sep = ',')
obs<- get.summary.stats_behav(dat)
obs.list<- df.to.list(obs)


#################################
#### Run Gibbs Sampler by ID ####
#################################

#basic settings
ngibbs=1000
nclustmax=20



## ID 1 ##

dat1.res=cluster.tsegm.behavior(dat=obs.list$`1`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat1.res$loglikel,type='l')
plot(dat1.res$phi[ngibbs,],type='h')


MAP1<- which(dat1.res$loglikel==max(dat1.res$loglikel))  # iteration 119 of MAP
tbsp.clust<- dat1.res$z[MAP1,]
table(tbsp.clust)  # 11 clusters




## ID 12 ##

dat12.res=cluster.tsegm.behavior(dat=obs.list$`12`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat12.res$loglikel,type='l')
plot(dat12.res$phi[ngibbs,],type='h')


MAP12<- which(dat12.res$loglikel==max(dat12.res$loglikel))  # iteration 73 of MAP
tbsp.clust<- dat12.res$z[MAP12,]
table(tbsp.clust)  # 8 clusters




## ID 19 ##

dat19.res=cluster.tsegm.behavior(dat=obs.list$`19`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                 gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat19.res$loglikel,type='l')
plot(dat19.res$phi[ngibbs,],type='h')


MAP19<- which(dat19.res$loglikel==max(dat19.res$loglikel))  # iteration 32 of MAP
tbsp.clust<- dat19.res$z[MAP19,]
table(tbsp.clust)  # 3 clusters




## ID 27 ##

dat27.res=cluster.tsegm.behavior(dat=obs.list$`27`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                 gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat27.res$loglikel,type='l')
plot(dat27.res$phi[ngibbs,],type='h')


MAP27<- which(dat27.res$loglikel==max(dat27.res$loglikel))  # iteration 1 & 300 of MAP
tbsp.clust<- dat27.res$z[MAP27[2],]
table(tbsp.clust)  # 2 clusters



#### All MAP values found before end of burn-in ####
