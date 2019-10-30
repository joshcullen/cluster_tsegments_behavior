rm(list=ls(all=TRUE))
set.seed(3)

#basic settings; l1 = # bins for TA; l2 = # bins for SL
nclust=8
l1=8; l2=6
n.tsegm=200
n=floor(runif(n.tsegm,min=20,max=200))

#parameters
theta1=matrix(runif(nclust*l1),nclust,l1)
theta1.true=theta1=theta1/apply(theta1,1,sum); apply(theta1,1,sum)
theta2=matrix(runif(nclust*l2),nclust,l2)
theta2.true=theta2=theta2/apply(theta2,1,sum); apply(theta2,1,sum)
theta3.true=theta3=runif(nclust)
z.true=z=sample(nclust,size=n.tsegm,replace=T)

#generate data
y1=matrix(NA,n.tsegm,l1)
y2=matrix(NA,n.tsegm,l2)
y3=rep(NA,n.tsegm)

for (i in 1:n.tsegm){
  y1[i,]=rmultinom(1,size=n[i],prob=theta1[z[i],])
  y2[i,]=rmultinom(1,size=n[i],prob=theta2[z[i],])
  y3[i]=rbinom(1,size=n[i],prob=theta3[z[i]])
}
image(y1[z==1,]/n[z==1])
image(y2[z==1,]/n[z==1])
hist(y3[z==3]/n[z==3])

#combinate datasets
colnames(y1)=paste0('y1.',1:l1)
colnames(y2)=paste0('y2.',1:l2)
fim=cbind(y1,y2,y3)

setwd('U:\\GIT_models\\cluster_tsegments_behavior')
write.csv(fim,'fake data.csv',row.names=F)

