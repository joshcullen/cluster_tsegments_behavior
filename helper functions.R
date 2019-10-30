df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
get.summary.stats_behav=function(dat){  #dat must have time.seg assigned; for all IDs
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  #calculate # of obs in each bin (per move param) by behav.seg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$behav.seg)
    
    
    #SL
    SL<- matrix(0, ntseg, max(dat.ind$SL, na.rm = T))
    colnames(SL)<- paste0("SL.",1:max(dat.ind$SL, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.ind %>% filter(behav.seg == j) %>% dplyr::select(SL) %>% table()
      SL[j,as.numeric(names(tmp))]<- tmp
    }
    
    #TA
    TA<- matrix(0, ntseg, max(dat.list[[i]]$TA, na.rm = T))
    colnames(TA)<- paste0("TA.",1:max(dat.list[[i]]$TA, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.list[[i]] %>% filter(behav.seg == j) %>% dplyr::select(TA) %>% table()
      TA[j,as.numeric(names(tmp))]<- tmp
    }
    
    #TAA
    TAA<- matrix(0, ntseg, 1)
    colnames(TAA)<- 'TAA.1'
    for (j in 1:ntseg){
      tmp<- dat.list[[i]] %>% filter(behav.seg == j) %>% dplyr::select(TAA) %>% table()
      tmp1<- tmp %>% data.frame() %>% filter(names(tmp) == 1) %>% dplyr::select(Freq) %>%
                     as.numeric()
      TAA[j,]<- tmp1
    }
    
    
    id<- rep(unique(dat.ind$id), ntseg)
    behav.res<- cbind(id, TA, SL, TAA) %>% data.frame()
    obs.list[[i]]<- behav.res
  }
  #obs<- do.call(rbind.data.frame, obs.list)
  obs<- map_dfr(obs.list, `[`)
  obs
}