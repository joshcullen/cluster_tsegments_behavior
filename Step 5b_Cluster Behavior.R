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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(reshape2)
library(stringr)
library(lubridate)

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
MAP<- which(res$loglikel==max(res$loglikel))  # iteration 711 (756 for 4 clust; 693 for 3 clust)
tbsp.clust<- res$z[MAP,]
table(tbsp.clust)

#Create DF of clusters by behavioral time segment
behav.seg<- 1:nrow(obs)
tbsp.clust<- cbind(tbsp.clust,behav.seg) %>% data.frame()




#####################
#### Plot output ####
#####################


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


#plot b3 (Resting/Encamped, very low SL, TA near pi, intermediate TAA)
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

#plot b5 (Resting2?, mainly low SL, TA closer to pi, relatively high TAA)
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


# Example of proportions using ID 1

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
rect.lims1.new$tbsp.clust<- factor(rect.lims1.new$tbsp.clust,
                                   labels = c("Transit","Exploratory","Resting","ARS","Resting2"))
rect.lims1.new$tbsp.clust<- factor(rect.lims1.new$tbsp.clust,
                                   levels(rect.lims1.new$tbsp.clust)[c(1:2,4:5,3)], ordered = T)


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
  scale_fill_manual("Behavior", values = ocean.amp(5)) +
  geom_hline(aes(yintercept = c(8.5, 14.5, 15.5)), color = "white", size = 0.5) +
  labs(x = "Time", y = "Behavior Profile") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))





# All IDs

tbsp.clust.time$tbsp.clust<- factor(tbsp.clust.time$tbsp.clust,
                                    labels = c("Transit","Exploratory","Resting","ARS","Resting2"))
tbsp.clust.time$tbsp.clust<- factor(tbsp.clust.time$tbsp.clust,
                                    levels(tbsp.clust.time$tbsp.clust)[c(1:2,4:5,3)], ordered = T)
tbsp.clust.time$id<- factor(tbsp.clust.time$id)
tbsp.clust.time$date<- dat$ESTtime %>% as_datetime()

#Window of peak breeding (March 1 - June 30)
breed<- data.frame(xmin = as_datetime(c("2016-03-01 00:00:00","2017-03-01 00:00:00",
                                        "2018-03-01 00:00:00","2019-03-01 00:00:00")),
                   xmax = as_datetime(c("2016-06-30 23:59:59","2017-06-30 23:59:59",
                                        "2018-06-30 23:59:59","2019-06-30 23:59:59")),
                   ymin = -Inf, ymax = Inf)

#Aligned by first observation
ggplot() +
  geom_tile(data=tbsp.clust.time, aes(x=time, y=id, fill=tbsp.clust)) +
  scale_fill_viridis_d("Behavior", direction = -1) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Time", y = "ID") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))


#Aligned by date
ggplot() +
  geom_rect(data = breed, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5) +
  geom_tile(data=tbsp.clust.time, aes(x=date, y=id, fill=tbsp.clust)) +
  scale_fill_viridis_d("Behavior", direction = -1) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_datetime(expand = c(0,0)) +
  labs(x = "Time", y = "ID") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))




## Separate Histograms of Thetas (not raw proportions) by Behavior

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
theta$tbsp.clust<- factor(theta$tbsp.clust,
                          labels = c("Transit","Exploratory","Resting","ARS","Resting2"))
theta$tbsp.clust<- factor(theta$tbsp.clust, levels(theta$tbsp.clust)[c(1:2,4:5,3)], ordered = T)
theta<- melt(theta, "tbsp.clust")
theta<- theta %>% mutate(param = ifelse(str_detect(theta$variable, "theta1"), "TA",
                                    ifelse(str_detect(theta$variable, "theta2"), "SL", "TAA")))

ggplot(theta, aes(x = variable, y = value, fill = tbsp.clust)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis_d(guide = F, direction = -1) +
  labs(x = "Bin", y = "Proportion") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  facet_grid(param ~ tbsp.clust, scales = "free_y")
  




## Show Proportion of Behaviors Associated with Each Time Segment (from posterior)

#sample z's from posterior (after burn-in since log-likelihood converged)
z.post<- res$z[501:ngibbs,]
z.post.count<- matrix(0, ncol(z.post), nclustmax)
tmp<-  apply(z.post, 2, table)
for (i in 1:nrow(z.post.count)) {
  z.post.count[i, as.integer(names(tmp[[i]]))]<- tmp[[i]]
}

#convert from count to proportion by time seg
z.post.prop<- apply(z.post.count, 1, function(x) x/sum(x)) %>% t()
colnames(z.post.prop)<- 1:5

#convert to long form for plotting
z.post.prop.time<- cbind(id = obs$id, z.post.prop) %>% data.frame() %>% slice(rep(1:n(), times = n))
z.post.list<- df.to.list(z.post.prop.time) %>% lapply(function(x) x[!(names(x) %in% c("id"))])
z.post.list.long<- lapply(z.post.list, function(x) {x %>% gather(key, value) %>%
    mutate(time=rep(1:nrow(x), times=ncol(x)))})
z.post.list.long<- lapply(z.post.list.long, function(x) {mutate_at(x, "key", as.factor)})



# ID 1 as example

ggplot() +
  geom_area(data=z.post.list.long[[1]], aes(x=time, y=value, fill=key)) +
  scale_fill_viridis_d("Behavior", direction = -1) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Time", y = "Proportion of Behavior", title = "ID 1") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16))






## Generate Geographic Maps with Pts Assigned a Behavior (from MAP estimate)

#load map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N


# ID 1

dat.list$`1`$behav<- tbsp.clust.time[tbsp.clust.time$id == 1, "tbsp.clust"]

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-120000), max(dat$utmlong+40000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat.list$`1`, aes(x=utmlong, y=utmlat), color="gray60", size=0.25) +
  geom_point(data = dat.list$`1`, aes(utmlong, utmlat, color=behav), size=1) +
  scale_color_viridis_d("Behavior", direction = -1) +
  labs(x = "Longitude", y = "Latitude", title = "ID 1") +
  theme_bw()


# ID 12

dat.list$`12`$behav<- tbsp.clust.time[tbsp.clust.time$id == 12, "tbsp.clust"]

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-120000), max(dat$utmlong+40000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat.list$`12`, aes(x=utmlong, y=utmlat), color="gray60", size=0.25) +
  geom_point(data = dat.list$`12`, aes(utmlong, utmlat, color=behav), size=1) +
  scale_color_viridis_d("Behavior", direction = -1) +
  labs(x = "Longitude", y = "Latitude", title = "ID 12") +
  theme_bw()


# ID 19

dat.list$`19`$behav<- tbsp.clust.time[tbsp.clust.time$id == 19, "tbsp.clust"]

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-120000), max(dat$utmlong+40000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat.list$`19`, aes(x=utmlong, y=utmlat), color="gray60", size=0.25) +
  geom_point(data = dat.list$`19`, aes(utmlong, utmlat, color=behav), size=1) +
  scale_color_viridis_d("Behavior", direction = -1) +
  labs(x = "Longitude", y = "Latitude", title = "ID 19") +
  theme_bw()


# ID 27

dat.list$`27`$behav<- tbsp.clust.time[tbsp.clust.time$id == 27, "tbsp.clust"]

ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-120000), max(dat$utmlong+40000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_path(data = dat.list$`27`, aes(x=utmlong, y=utmlat), color="gray60", size=0.25) +
  geom_point(data = dat.list$`27`, aes(utmlong, utmlat, color=behav), size=1) +
  scale_color_manual("Behavior", values = viridis(5)[c(3,1)]) +
  labs(x = "Longitude", y = "Latitude", title = "ID 27") +
  theme_bw()





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
