plot(res$loglikel,type='l')

plot(res$phi[ngibbs,],type='h')

fim=data.frame(zestim=res$z[ngibbs,],ztrue=z.true)
table(fim)

<<<<<<< HEAD
seq1=c(6,5,2,7,3,8,4,1)
tmp=matrix(res$theta1[ngibbs,],nclustmax,l1)
theta1.estim=tmp[seq1,]
rango=range(c(theta1.estim,theta1.true))
plot(theta1.estim[1:length(seq1),],theta1.true,xlim=rango,ylim=rango)
lines(rango,rango)
=======
seq1=c(4,5,6,2,10,3,8,9,1,7)
# tmp=matrix(res$theta1[ngibbs,],nclustmax,nloc)
# theta.estim=tmp[seq1,]
# rango=range(c(theta.estim,theta.true))
# plot(theta.estim[1:length(seq1),],theta.true,xlim=rango,ylim=rango)
# lines(rango,rango)
>>>>>>> 63fa85df4d961e25479a8c2fd6a84faf38fb6e37

tmp=matrix(res$theta2[ngibbs,],nclustmax,l2)
theta2.estim=tmp[seq1,]
rango=range(c(theta2.estim,theta2.true))
plot(theta2.estim[1:length(seq1),],theta2.true,xlim=rango,ylim=rango)
lines(rango,rango)
