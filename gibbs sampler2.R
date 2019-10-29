cluster.tsegm.behavior=function(dat,a.theta3,b.theta3,psi,gamma1,nclustmax,ngibbs,nburn){
  ind=grep('TA',colnames(dat)); l1=length(ind); y1=data.matrix(dat[,ind])
  ind=grep('SL',colnames(dat)); l2=length(ind); y2=data.matrix(dat[,ind])
  n.tsegm=nrow(dat)
  y3=matrix(dat[,'TAA'],n.tsegm,1)
  n=apply(y2,1,sum)
  
  #general settings
  lo=0.000000000000001
  hi=1-lo
  
  #starting values
  z=sample(1:nclustmax,size=n.tsegm,replace=T)
  theta1=matrix(1/l1,nclustmax,l1)
  theta2=matrix(1/l2,nclustmax,l2)
  theta3=rep(0.5,nclustmax)
  phi=rep(1/nclustmax,nclustmax)
  
  #store results
  store.phi=matrix(NA,ngibbs,nclustmax)
  store.z=matrix(NA,ngibbs,n.tsegm)
  store.theta1=matrix(NA,ngibbs,nclustmax*l1)
  store.theta2=matrix(NA,ngibbs,nclustmax*l2)
  store.theta3=matrix(NA,ngibbs,nclustmax)
  store.loglikel=matrix(NA,ngibbs,1)
  
  #gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    
    #occasionally re-order this
    if (i<nburn & i%%50==0){
      ind=order(phi,decreasing=T)
      theta1=theta1[ind,]
      theta2=theta2[ind,]
      theta3=theta3[ind]
      phi=phi[ind]
      znew=z
      for (j in 1:nclustmax){
        znew[z==ind[j]]=j
      }
      z=znew
    }
    
    #draw samples from FCD's
    v=sample.v(z=z,nclustmax=nclustmax,gamma1=gamma1)
    phi=GetPhi(vec=c(v,1),nclustmax=nclustmax)
    
    theta1=sample.theta1(y1=y1,nclustmax=nclustmax,l1=l1,z=z,psi=psi)
    theta1[theta1<lo]=lo #to avoid numerical issues
    # theta1=theta1.true
    
    theta2=sample.theta2(y2=y2,nclustmax=nclustmax,l2=l2,z=z,psi=psi)
    theta2[theta2<lo]=lo #to avoid numerical issues
    # theta2=theta2.true
    
    theta3=sample.theta3(y3=y3,n=n,nclustmax=nclustmax,z=z,a.theta3=a.theta3,b.theta3=b.theta3)
    theta3[theta3<lo]=lo #to avoid numerical issues
    theta3[theta3>hi]=hi #to avoid numerical issues
    # theta3=theta3.true
    
    z=sample.z(y1=y1,y2=y2,y3=y3,n=n,
               ltheta1=log(theta1),ltheta2=log(theta2),ltheta3=log(theta3),l.1minus.theta3=log(1-theta3),
               lphi=log(phi),z=z,
               n.tsegm=n.tsegm,nclustmax=nclustmax,l1=l1,l2=l2,
               a.theta3=a.theta3,b.theta3=b.theta3,psi=psi)
    # z=z.true
    
    #get logl
    logl=sum(y1*log(theta1)[z,])+
      sum(y2*log(theta2)[z,])+
      sum(y3*log(theta3)[z]+(n-y3)*log(1-theta3[z]))
    lprior=sum(dbeta(v,1,gamma1,log=T))+sum((psi-1)*log(theta1))+sum((psi-1)*log(theta2))+
      sum((a.theta3-1)*log(theta3)+(b.theta3-1)*log(1-theta3))
    
    #store results
    store.loglikel[i]=logl+lprior
    store.theta1[i,]=theta1
    store.theta2[i,]=theta2
    store.theta3[i,]=theta3
    store.phi[i,]=phi
    store.z[i,]=z
  }
  list(loglikel=store.loglikel,theta1=store.theta1,theta2=store.theta2,theta3=store.theta3,
       phi=store.phi,z=store.z)
}
