set.seed(1)

library('Rcpp')
library('MCMCpack')
library(dplyr)
library(purrr)
library(tidyr) #for gather function
library(ggnewscale) #for multiple fill scales in ggplot2
library(pals) # for more color palettes
library(progress) #for progress bar

#get functions
sourceCpp('aux1.cpp')
source('gibbs functions.R')
source('gibbs sampler2.R')
source('helper functions.R')

#get data
dat<- read.csv('Snail Kite Gridded Data_behav.csv', header = T, sep = ',')
dat.list<- df.to.list(dat)
obs<- get.summary.stats_behav(dat)
obs.list<- df.to.list(obs)


#################################
#### Run Gibbs Sampler by ID ####
#################################

#basic settings
ngibbs=1000
nclustmax=6



## ID 1 ##

#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

dat1.res=cluster.tsegm.behavior(dat=obs.list$`1`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat1.res$loglikel,type='l')
plot(dat1.res$phi[ngibbs,],type='h')


MAP1<- which(dat1.res$loglikel==max(dat1.res$loglikel))  # iteration 581 of MAP
tbsp.clust1<- dat1.res$z[MAP1,]
table(tbsp.clust1)

time.seg<- 1:nrow(obs.list$`1`)
tbsp.clust1<- cbind(tbsp.clust1,time.seg) %>% data.frame()
dat.list$`1`<- left_join(dat.list$`1`, tbsp.clust1, by="time.seg")  #### ADD TIME.SEG TO DAT


## Plot heatmap of clusters

#format data
colnames(obs.list$`1`)[-1]=1:ncol(obs.list$`1`[,-1]) %>% as.character()
nobs=nrow(obs.list$`1`)
nloc=ncol(obs.list$`1`[,-1])
obs1.long<- obs.list$`1`[,-1] %>% data.frame() %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs1.long$key<- as.factor(obs1.long$key)
levels(obs1.long$key)<- 1:nloc
obs1.long$key<- as.numeric(obs1.long$key)

tbsp.clust1[,1]<- tbsp.clust1[,1] %>% as.numeric()
tbsp.clust1[,2]<- tbsp.clust1[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust1$tbsp.clust1)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust1")
rect.lims.new<- data.frame(rect.lims.new)

#plot
ggplot() +
  geom_tile(data=obs1.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax, ymin = max(obs1.long$key) + 0.5,
                                    ymax = max(obs1.long$key) + 0.75, fill = tbsp.clust1),
            color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  labs(x = "Time Segment", y = "Activity Center") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))




## ID 12 ##

pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

dat12.res=cluster.tsegm.behavior(dat=obs.list$`12`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(dat12.res$loglikel,type='l')
plot(dat12.res$phi[ngibbs,],type='h')


# MAP12<- which(dat12.res$loglikel==max(dat12.res$loglikel))  # iteration 128 of MAP
MAP12<- dat12.res$loglikel %>% order(decreasing = T) %>% subset(. > 500) %>% first() #iteration 591
tbsp.clust12<- dat12.res$z[MAP12,]
table(tbsp.clust12)  




## ID 19 ##

pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

dat19.res=cluster.tsegm.behavior(dat=obs.list$`19`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                 gamma1=gamma1,nclustmax=3,ngibbs=ngibbs)

plot(dat19.res$loglikel,type='l')
plot(dat19.res$phi[ngibbs,],type='h')


MAP19<- which(dat19.res$loglikel==max(dat19.res$loglikel))  # iteration 838 of MAP
tbsp.clust19<- dat19.res$z[MAP19,]
table(tbsp.clust19)  # 3 clusters




## ID 27 ##

pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

dat27.res=cluster.tsegm.behavior(dat=obs.list$`27`[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                 gamma1=gamma1,nclustmax=2,ngibbs=ngibbs)

plot(dat27.res$loglikel,type='l')
plot(dat27.res$phi[ngibbs,],type='h')


# MAP27<- which(dat27.res$loglikel==max(dat27.res$loglikel))  # iteration 1 & 300 of MAP
# tbsp.clust27<- dat27.res$z[MAP27[2],]
# table(tbsp.clust27)  # 2 clusters

### APPEARS TO ONLY BE ASSOCIATED WITH A SINGLE CLUSTER



#### All MAP values found before end of burn-in ####
