plot(res$loglikel,type='l')

plot(res$phi[ngibbs,],type='h')

fim=data.frame(zestim=res$z[ngibbs,],ztrue=z.true)
table(fim)

seq1=c(4,5,6,2,10,3,8,9,1,7)
# tmp=matrix(res$theta1[ngibbs,],nclustmax,nloc)
# theta.estim=tmp[seq1,]
# rango=range(c(theta.estim,theta.true))
# plot(theta.estim[1:length(seq1),],theta.true,xlim=rango,ylim=rango)
# lines(rango,rango)

