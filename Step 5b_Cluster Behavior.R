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
MAP<- which(res$loglikel==max(res$loglikel))  # iteration 711
tbsp.clust<- res$z[MAP,]
table(tbsp.clust)

#Create DF of clusters by behavioral time segment
behav.seg<- 1:nrow(obs)
tbsp.clust<- cbind(tbsp.clust,behav.seg) %>% data.frame()




#####################
#### Plot output ####
#####################


## All IDs Clustered by Time Segment
#format data
nobs=nrow(obs)
nloc=ncol(obs[,-1])
obs.long<- obs[,-1] %>% data.frame() %>% gather(key, value) %>%
           mutate(behav.seg=rep(tbsp.clust$behav.seg, times=nloc))
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
  geom_tile(data=obs.long, aes(x=behav.seg, y=key, fill=log10(value+1))) +
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






## Separate by Behavior to Compare Distribs (Using Proportions) Over Time

#convert from raw obs to proportion
y1<- obs[,grep('y1',colnames(obs))]
y1.prop<- apply(y1, 1, function(x) {x/sum(x)}) %>% t()

y2<- obs[,grep('y2',colnames(obs))]
y2.prop<- apply(y2, 1, function(x) {x/sum(x)}) %>% t()

y3<- obs[,grep('y3',colnames(obs))]
n<- apply(y2,1,sum)
y3.prop<- y3/n

obs.prop<- cbind(id = obs$id, y1.prop, y2.prop, y3 = y3.prop) %>% data.frame()


#format data
nobs=nrow(obs.prop)
nloc=ncol(obs.prop[,-1])
obs.prop.long<- obs.prop[,-1] %>% data.frame() %>% gather(key, value) %>%
  mutate(behav.seg=rep(1:nobs, times=nloc))  #time is now obs #, not segment #
obs.prop.long$key<- as.factor(obs.prop.long$key)


#separate behavior DFs
b<- left_join(obs.prop.long, tbsp.clust, by = "behav.seg")
b1<- b %>% filter(tbsp.clust == 1)
b2<- b %>% filter(tbsp.clust == 2)
b3<- b %>% filter(tbsp.clust == 3)
b4<- b %>% filter(tbsp.clust == 4)
b5<- b %>% filter(tbsp.clust == 5)

#plot b1 (Directed; large SL, TA near 0, high TAA)
ggplot() +
  geom_tile(data=b1, aes(x=as.factor(behav.seg), y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_hline(aes(yintercept = c(8.5, 14.5)), color = "white", size = 0.35) +
  labs(x = "Time Segment", y = "Behavior Profile", title = "Behavior 1") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


#plot b2 (Exploratory; intermediate SL, varied TA, intermediate TAA)
ggplot() +
  geom_tile(data=b2, aes(x=as.factor(behav.seg), y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_hline(aes(yintercept = c(8.5, 14.5)), color = "white", size = 0.35) +
  labs(x = "Time Segment", y = "Behavior Profile", title = "Behavior 2") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


#plot b3 (Encamped, very low SL, TA near pi, intermediate TAA)
ggplot() +
  geom_tile(data=b3, aes(x=as.factor(behav.seg), y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_hline(aes(yintercept = c(8.5, 14.5)), color = "white", size = 0.35) +
  labs(x = "Time Segment", y = "Behavior Profile", title = "Behavior 3") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


#plot b4 (ARS, mainly low SL, varied TA, relatively high TAA)
ggplot() +
  geom_tile(data=b4, aes(x=as.factor(behav.seg), y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_hline(aes(yintercept = c(8.5, 14.5)), color = "white", size = 0.35) +
  labs(x = "Time Segment", y = "Behavior Profile", title = "Behavior 4") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))

#plot b5 (ARS2?, mainly low SL, TA closer to pi, relatively high TAA)
ggplot() +
  geom_tile(data=b5, aes(x=as.factor(behav.seg), y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_hline(aes(yintercept = c(8.5, 14.5)), color = "white", size = 0.35) +
  labs(x = "Time Segment", y = "Behavior Profile", title = "Behavior 5") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## Separate by ID Showing Time Series

#create DF with matching proportional distrib by time segment for all observations
obs.prop.time<- obs.prop %>% slice(rep(1:n(), times = n))


#format data
obs.prop.list<- df.to.list(obs.prop.time) %>% lapply(function(x) x[!(names(x) %in% c("id"))])
obs.prop.list.long<- lapply(obs.prop.list, function(x) {x %>% gather(key, value) %>%
    mutate(time=rep(1:nrow(x), times=ncol(x)))})
obs.prop.list.long<- lapply(obs.prop.list.long, function(x) {mutate_at(x, "key", as.factor)})

time<- lapply(obs.prop.list, function(x) data.frame(time = 1:nrow(x))) %>% map_dfr(`[`)
tbsp.clust.time<- tbsp.clust %>% slice(rep(1:n(), times = n)) %>%
  mutate(time = time[,1]) %>% mutate(id = dat$id)



# ID 1

#generate boxes denoting clusters
rect.lims1<- rle(tbsp.clust.time[tbsp.clust.time$id == 1, "tbsp.clust"])
rect.lims1$lengths<- cumsum(rect.lims1$lengths)+0.5
rect.lims1$lengths<- c(0.5, rect.lims1$lengths)

rect.lims1.new<- matrix(0, length(rect.lims1$values), 3)
for (i  in 2:length(rect.lims1$lengths)) {
  rect.lims1.new[i-1,]<- c(rect.lims1$lengths[i-1], rect.lims1$lengths[i], rect.lims1$values[i-1])
}
colnames(rect.lims1.new)<- c("xmin","xmax","tbsp.clust")
rect.lims1.new<- data.frame(rect.lims1.new)


ggplot() +
  geom_tile(data=obs.prop.list.long[[1]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims1.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims1.new, aes(xmin = xmin, xmax = xmax,
                                    ymin = max(as.integer(obs.prop.list.long[[1]]$key)) + 0.5,
                                    ymax = max(as.integer(obs.prop.list.long[[1]]$key)) + 1,
                                    fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 12

#generate boxes denoting clusters
rect.lims12<- rle(tbsp.clust.time[tbsp.clust.time$id == 12, "tbsp.clust"])
rect.lims12$lengths<- cumsum(rect.lims12$lengths)+0.5
rect.lims12$lengths<- c(0.5, rect.lims12$lengths)

rect.lims12.new<- matrix(0, length(rect.lims12$values), 3)
for (i  in 2:length(rect.lims12$lengths)) {
  rect.lims12.new[i-1,]<- c(rect.lims12$lengths[i-1], rect.lims12$lengths[i],
                            rect.lims12$values[i-1])
}
colnames(rect.lims12.new)<- c("xmin","xmax","tbsp.clust")
rect.lims12.new<- data.frame(rect.lims12.new)


ggplot() +
  geom_tile(data=obs.prop.list.long[[2]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims12.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims12.new, aes(xmin = xmin, xmax = xmax,
                                     ymin = max(as.integer(obs.prop.list.long[[2]]$key)) + 0.5,
                                     ymax = max(as.integer(obs.prop.list.long[[2]]$key)) + 1,
                                     fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 19

#generate boxes denoting clusters
rect.lims19<- rle(tbsp.clust.time[tbsp.clust.time$id == 19, "tbsp.clust"])
rect.lims19$lengths<- cumsum(rect.lims19$lengths)+0.5
rect.lims19$lengths<- c(0.5, rect.lims19$lengths)

rect.lims19.new<- matrix(0, length(rect.lims19$values), 3)
for (i  in 2:length(rect.lims19$lengths)) {
  rect.lims19.new[i-1,]<- c(rect.lims19$lengths[i-1], rect.lims19$lengths[i],
                            rect.lims19$values[i-1])
}
colnames(rect.lims19.new)<- c("xmin","xmax","tbsp.clust")
rect.lims19.new<- data.frame(rect.lims19.new)


ggplot() +
  geom_tile(data=obs.prop.list.long[[3]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims19.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims19.new, aes(xmin = xmin, xmax = xmax,
                                      ymin = max(as.integer(obs.prop.list.long[[3]]$key)) + 0.5,
                                      ymax = max(as.integer(obs.prop.list.long[[3]]$key)) + 1,
                                      fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 27

#generate boxes denoting clusters
rect.lims27<- rle(tbsp.clust.time[tbsp.clust.time$id == 27, "tbsp.clust"])
rect.lims27$lengths<- cumsum(rect.lims27$lengths)+0.5
rect.lims27$lengths<- c(0.5, rect.lims27$lengths)

rect.lims27.new<- matrix(0, length(rect.lims27$values), 3)
for (i  in 2:length(rect.lims27$lengths)) {
  rect.lims27.new[i-1,]<- c(rect.lims27$lengths[i-1], rect.lims27$lengths[i],
                            rect.lims27$values[i-1])
}
colnames(rect.lims27.new)<- c("xmin","xmax","tbsp.clust")
rect.lims27.new<- data.frame(rect.lims27.new)


ggplot() +
  geom_tile(data=obs.prop.list.long[[4]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims27.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims27.new, aes(xmin = xmin, xmax = xmax,
                                      ymin = max(as.integer(obs.prop.list.long[[4]]$key)) + 0.5,
                                      ymax = max(as.integer(obs.prop.list.long[[4]]$key)) + 1,
                                      fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





## Separate by ID Showing Time Series of Thetas (not raw proportions)

#take mean of posterior theta estimates (after nburn since log likelihood converged)
theta1.post<- res$theta1[501:ngibbs,] %>% apply(2, mean)
theta2.post<- res$theta2[501:ngibbs,] %>% apply(2, mean)
theta3.post<- res$theta3[501:ngibbs,] %>% apply(2, mean)

#create DF with matching thetas per behavior for all observations
theta1<- matrix(theta1.post, nclustmax, ncol(y1), byrow = F)
theta2<- matrix(theta2.post, nclustmax, ncol(y2), byrow = F)
theta3<- matrix(theta3.post, nclustmax, 1, byrow = F)
colnames(theta1)=paste0('theta1.',1:ncol(y1))
colnames(theta2)=paste0('theta2.',1:ncol(y2))
colnames(theta3)='theta3'

theta<- cbind(tbsp.clust = 1:5, theta1, theta2, theta3) %>% data.frame()
obs.theta.time<- left_join(tbsp.clust.time, theta, by = "tbsp.clust")


#format data
obs.theta.list<- df.to.list(obs.theta.time) %>% lapply(function(x) x[!(names(x) %in% c("tbsp.clust","behav.seg","time","id"))])
obs.theta.list.long<- lapply(obs.theta.list, function(x) {x %>% gather(key, value) %>%
    mutate(time=rep(1:nrow(x), times=ncol(x)))})
obs.theta.list.long<- lapply(obs.theta.list.long, function(x) {mutate_at(x, "key", as.factor)})



# ID 1

ggplot() +
  geom_tile(data=obs.theta.list.long[[1]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims1.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims1.new, aes(xmin = xmin, xmax = xmax,
                                     ymin = max(as.integer(obs.theta.list.long[[1]]$key)) + 0.5,
                                     ymax = max(as.integer(obs.theta.list.long[[1]]$key)) + 1,
                                     fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile", title = "ID 1") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 12

ggplot() +
  geom_tile(data=obs.theta.list.long[[2]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims12.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims12.new, aes(xmin = xmin, xmax = xmax,
                                      ymin = max(as.integer(obs.theta.list.long[[2]]$key)) + 0.5,
                                      ymax = max(as.integer(obs.theta.list.long[[2]]$key)) + 1,
                                      fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile", title = "ID 12") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 19

ggplot() +
  geom_tile(data=obs.theta.list.long[[3]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims19.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims19.new, aes(xmin = xmin, xmax = xmax,
                                      ymin = max(as.integer(obs.theta.list.long[[3]]$key)) + 0.5,
                                      ymax = max(as.integer(obs.theta.list.long[[3]]$key)) + 1,
                                      fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile", title = "ID 19") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


# ID 27

ggplot() +
  geom_tile(data=obs.theta.list.long[[4]], aes(x=time, y=key, fill=value)) +
  scale_fill_viridis_c("Proportion") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  new_scale_fill() +
  # geom_vline(data = rect.lims27.new, aes(xintercept = xmin), color = "white", size = 0.35) +
  geom_rect(data=rect.lims27.new, aes(xmin = xmin, xmax = xmax,
                                      ymin = max(as.integer(obs.theta.list.long[[4]]$key)) + 0.5,
                                      ymax = max(as.integer(obs.theta.list.long[[4]]$key)) + 1,
                                      fill = tbsp.clust), color = NA, size = 1.5) +
  scale_fill_gradientn("Behavior", colours = ocean.amp(6)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile", title = "ID 27") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))






## Show Proportion of Behaviors Associated with Each Time Segment (from posterior)



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
