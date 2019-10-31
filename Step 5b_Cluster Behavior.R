set.seed(1)

library('Rcpp')
library('MCMCpack')
library(dplyr)
library(purrr)
library(tidyr) #for gather function
library(ggplot2)
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


#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

#basic settings
ngibbs=1000
nclustmax=5



#progress bar
pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

# Run Gibbs sampler
res=cluster.tsegm.behavior(dat=obs[,-1],a.theta3=a.theta3,b.theta3=b.theta3,psi=psi,
                                gamma1=gamma1,nclustmax=nclustmax,ngibbs=ngibbs)

plot(res$loglikel,type='l')
plot(res$phi[ngibbs,],type='h')

#Find MAP value for assignment of clusters
MAP<- which(res$loglikel==max(res$loglikel))  # iteration 748
tbsp.clust<- res$z[MAP,]
table(tbsp.clust)

#Create DF of clusters by behavioral time segment
behav.seg<- 1:nrow(obs)
tbsp.clust<- cbind(tbsp.clust,behav.seg) %>% data.frame()





## Plot heatmap of clusters

#format data
nobs=nrow(obs)
nloc=ncol(obs[,-1])
obs.long<- obs[,-1] %>% data.frame() %>% gather(key, value) %>%
           mutate(time=rep(tbsp.clust$behav.seg, times=nloc))
obs.long$key<- as.factor(obs.long$key)

tbsp.clust[,1]<- tbsp.clust[,1] %>% as.numeric()
tbsp.clust[,2]<- tbsp.clust[,2] %>% as.numeric()


#generate boxes denoting clusters
rect.lims<- rle(tbsp.clust$tbsp.clust)
rect.lims$lengths<- cumsum(rect.lims$lengths)+0.5
rect.lims$lengths<- c(0.5, rect.lims$lengths)

rect.lims.new<- matrix(0, length(rect.lims$values), 3)
for (i  in 2:length(rect.lims$lengths)) {
  rect.lims.new[i-1,]<- c(rect.lims$lengths[i-1], rect.lims$lengths[i], rect.lims$values[i-1])
}
colnames(rect.lims.new)<- c("xmin","xmax","tbsp.clust")
rect.lims.new<- data.frame(rect.lims.new)


#generate boxes denoting IDs
id.rect<- rle(obs$id)
id.rect$lengths<- cumsum(id.rect$lengths) + 0.5
id.rect$lengths<- c(0.5, id.rect$lengths)

id.rect.new<- matrix(0, length(id.rect$values), 3)
for (i  in 2:length(id.rect$lengths)) {
  id.rect.new[i-1,]<- c(id.rect$lengths[i-1], id.rect$lengths[i], id.rect$values[i-1])
}
colnames(id.rect.new)<- c("xmin","xmax","id")
id.rect.new<- data.frame(id.rect.new)
id.rect.new$id<- as.factor(id.rect.new$id)


#plot
ggplot() +
  geom_tile(data=obs.long, aes(x=time, y=key, fill=log10(value+1))) +
  scale_fill_viridis_c("log10(N+1)") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  geom_vline(data = rect.lims.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims.new, aes(xmin = xmin, xmax = xmax,
                                    ymin = max(as.integer(obs.long$key)) + 0.5,
                                    ymax = max(as.integer(obs.long$key)) + 1,
                                    fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Time Cluster", colours = ocean.amp(6)) +
  new_scale_fill() +
  geom_rect(data = id.rect.new, aes(xmin = xmin, xmax = xmax,
                                    ymin = max(as.integer(obs.long$key)) + 1,
                                    ymax = max(as.integer(obs.long$key)) + 1.5,
                                    fill = id), color = NA, size = 1.5) +
  scale_fill_viridis_d("ID", option = "plasma") +
  labs(x = "Time Segment", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))




###################################
#### Export Updated Data Frame ####
###################################
reps<- obs %>% group_by(id) %>% tally()
for (i in 1:nrow(reps)) {
  tmp<- seq(1, reps$n[i])
  tbsp.clust$behav.seg[which(obs$id == reps$id[i])]<- tmp
}
tbsp.clust<- cbind(tbsp.clust, id = obs$id) %>% data.frame()
 

#Add cluster assignments to original data
for (i in 1:length(dat.list)) {
  tmp<- tbsp.clust %>% filter(id == unique(dat.list[[i]]$id)) %>% dplyr::select(-id)
  dat.list[[i]]<- left_join(dat.list[[i]], tmp, by="behav.seg")
}

dat2<- map_dfr(dat.list, `[`)
# write.csv(dat2, "Snail Kite Gridded Data_behavClust.csv", row.names = F)
