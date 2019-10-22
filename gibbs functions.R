sample.theta1=function(y1,nclustmax,l1,z,psi){
  theta1=matrix(NA,nclustmax,l1)
  for (i in 1:nclustmax){
    cond=z==i
    soma=sum(cond)
    if (soma==0) tmp=rep(0,l1)
    if (soma==1) tmp=y1[cond,]
    if (soma >1) tmp=colSums(y1[cond,])
    theta1[i,]=rdirichlet(1,tmp+psi)
  }
  theta1
}
#----------------------------------------------
sample.theta2=function(y2,nclustmax,l2,z,psi){
  theta2=matrix(NA,nclustmax,l2)
  for (i in 1:nclustmax){
    cond=z==i
    soma=sum(cond)
    if (soma==0) tmp=rep(0,l2)
    if (soma==1) tmp=y2[cond,]
    if (soma >1) tmp=colSums(y2[cond,])
    theta2[i,]=rdirichlet(1,tmp+psi)
  }
  theta2
}
#----------------------------------------------
sample.theta3=function(y3,n,nclustmax,z,a.theta3,b.theta3){
  n1=n0=rep(NA,nclustmax)
  for (i in 1:nclustmax){
    cond=z==i
    soma=sum(cond)
    if (soma==0) {n1[i]=0;             n0[i]=0}
    if (soma==1) {n1[i]=y3[cond];      n0[i]=n[cond]-n1[i]}
    if (soma >1) {n1[i]=sum(y3[cond]); n0[i]=sum(n[cond])-n1[i]}
  }
  rbeta(nclustmax,n1+a.theta3,n0+b.theta3)
}
#----------------------------------------------
sample.v=function(z,nclustmax,gamma1){
  #get n
  n=rep(0,nclustmax)
  tmp=table(z)
  n[as.numeric(names(tmp))]=tmp

  #get ngreater
  seq1=nclustmax:1
  tmp=cumsum(n[seq1])[seq1]
  ngreater=tmp[-1]
  
  #get v's
  rbeta(nclustmax-1,n[-nclustmax]+1,ngreater+gamma1)
}
#----------------------------------------------
sample.z=function(y1,y2,y3,
                  ltheta1,ltheta2,ltheta3,l.1minus.theta3,
                  lphi,z,n.tsegm,nclustmax,l1,l2,n,a.theta3,b.theta3,psi){

  #for clusters that exist
  lprob=matrix(NA,n.tsegm,nclustmax)
  for (i in 1:nclustmax){
    tmp1=y1*matrix(ltheta1[i,],n.tsegm,l1,byrow=T)
    tmp2=y2*matrix(ltheta2[i,],n.tsegm,l2,byrow=T)
    tmp3=y3*ltheta3[i]+(n-y3)*l.1minus.theta3[i]
    lprob[,i]=rowSums(tmp1)+rowSums(tmp2)+tmp3+lphi[i]
  }

  #sample each z separately
  p2=lgamma(l1*psi)-l1*lgamma(psi)  
  p4=lgamma(l2*psi)-l2*lgamma(psi)
  p6=lgamma(a.theta3+b.theta3)-lgamma(a.theta3)-lgamma(b.theta3)
  
  for (i in 1:n.tsegm){
    max.z=max(z)
    lprob1=lprob[i,1:max.z]
    
    #for clusters that do not exist
    if (max.z!=nclustmax){
      new.z=max.z+1
      p1=lphi[new.z]
      p3=sum(lgamma(y1[i,]+psi))-lgamma(n[i]+l1*psi)
      p5=sum(lgamma(y2[i,]+psi))-lgamma(n[i]+l2*psi)
      p7=lgamma(y3[i]+a.theta3)+lgamma(n[i]-y3[i]+b.theta3)-lgamma(a.theta3+n[i]+b.theta3)
      lprob1=c(lprob1,NA)
      lprob1[new.z]=p1+p2+p3+p4+p5+p6+p7
    }
    tmp=lprob1-max(lprob1)
    tmp=exp(tmp)
    prob=tmp/sum(tmp)
    tmp=rmultinom(1,size=1,prob=prob)
    z[i]=which(tmp==1)
  }
    
  z
}
