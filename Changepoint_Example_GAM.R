#### uisng the GAM analysis from the Rosenberge et al. 2019 paper to estimate the change in trends
library(tidyverse)
library(ggrepel)
library(ggforce)
library(rjags)
library(jagsUI)
library(mgcv)
source("utility_functions.R")

ind = "index"
lci = "lci"
uci = "uci"
year = "year"
base.yr = 1970
popsource = "Pop.source"

set.seed(2019)




indicesraw = read.csv("Rosenberg et al annual indices of abundance.csv",
                      stringsAsFactors = F)

indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]
# ###################### GAM smoothing of indices

indicesraw <- indicesraw %>% 
  mutate(lind = log(index.raw),
         lsd = (log(uci.raw)-log(lci.raw))/(1.96*2))

meanlsd = mean(indicesraw$lsd,na.rm = T)
indicesraw[which(is.na(indicesraw$lsd)),"lsd"] <- meanlsd


sps = unique(indicesraw$species)

nyears_recent = 10

jj = 0
for(ss in sps[464:length(sps)]){
  jj = jj+1
  
  wss = which(indicesraw$species == ss)
  tmp = indicesraw[wss,]
  
  fyr = unique(tmp$firstyear)
  lyr = unique(tmp$lastyear)
  
  torep = which(indicesraw$species == ss & (indicesraw$year >= fyr & indicesraw$year <= lyr))
  
  if(any(!is.na(tmp$lind))){
    wprec = T
  }else{
    wprec = F
  }
  
  tmpd = tmp[which(!is.na(tmp$lind)),]
  
  
  
  nknots = min(c(11,max(floor(nrow(tmpd)/3),3)))
  if(ss %in% c("Cackling Goose","Greater White-fronted Goose","Trumpeter Swan")){nknots = 4} 
  
  
  
  form = as.formula(paste("lind","~",
                          "s(year,k =",nknots,")"))
  
  
  
  
  ncounts = nrow(tmpd)
  ###### building the GAM basis function
  # gam basis functions created using function jagam from mgcv package          
  
  yminy = min(tmpd$year,na.rm = T)
  ymaxy = max(tmpd$year,na.rm = T)
  yearvec = tmpd$year-(yminy-1)
  yrs = seq(yminy,ymaxy,by = 1)
  nyears = length(yrs)
  ymin = min(yearvec)
  ymax = max(yearvec)
  lindex = tmpd$lind
  preci = 1/(tmpd$lsd^2)  
  
  preddat = data.frame(lind = 1,
                       year = yrs)
  
  
  

    
    gamprep = jagam(formula = form,
                    data = tmpd,
                    file = "tempgam.txt",
                    centred = T)
    
    gamprep.pred = jagam(formula = form,
                         data = preddat,
                         file = "tempgampred.txt",
                         centred = T)
    
    
    
    dat = list(X = gamprep$jags.data$X,
               S1 = gamprep$jags.data$S1,
               ncounts = nrow(tmpd),
               lindex = lindex,
               nknots = nknots,
               preci = preci,
               X.pred = gamprep.pred$jags.data$X,
               nyears = nyears,
               zero = gamprep$jags.data$zero,
               Ys = as.integer(ymin),
               Ye = as.integer(ymax),
               Yb = as.integer(ymax-nyears_recent))
    
    
    mgo = jagsUI(data = dat,
                 model.file = paste0("JAGS_model_GAM.R"),
                 n.chains = 3,
                parameters.to.save = c("ind.pred",
                         "C1",
                         "C2",
                         "T1",
                         "T2",
                         "Tdif",
                         "Tdif_neg"
                         #"rho"
                         #"mu",
                         #"b"
                       ),
                       n.iter = 20000,
                       n.burnin = 10000,
                       n.thin = 10,
                parallel = T)
    
    
    mgosum = data.frame(mgo$summary)
    
    
    mgosum$parameter <- row.names(mgosum)
    
    mgosum$species = ss
    
    if(jj == 1){
      out = mgosum
    }else{
      out <- bind_rows(out,mgosum)
    }
    
    

}

names(out)[which(grepl(names(out),pattern = ".",fixed = T))] <- paste0(gsub(names(out)[which(grepl(names(out),pattern = ".",fixed = T))] ,pattern = ".",replacement = "_", fixed = TRUE))

mu_ <- filter(out,parameter %in% c(paste0("ind.pred[",1:53,"]")))
mu_ <- mutate(mu_,yr = jags_dim(dat = mu_, cl = "parameter",var = "ind.pred"),
              est = "mu")
fyr_sp <- unique(indicesraw[,c("species","firstyear")])
mu_ <- left_join(mu_,fyr_sp,by = "species")
mu_$year <- mu_$yr+(mu_$firstyear-1)






Tdif <- filter(out,parameter %in% c(paste0("Tdif")))
Tdif_neg <- filter(out,parameter %in% c(paste0("Tdif_neg")))

T1 <- filter(out,parameter %in% c(paste0("T1")))
T2 <- filter(out,parameter %in% c(paste0("T2")))


prob_annot <- select(Tdif_neg,mean,species)
names(prob_annot)[1] <- "mean_p"
bdif <- select(Tdif,X2_5_,mean,X97_5_,species)
bdif <- left_join(bdif,prob_annot,by = "species")
bdif$Difference_late_minus_early_trend <- round(bdif$mean,1)
bdif$UCI_90_Difference_late_minus_early_trend <- round(bdif$X97_5_,1)
bdif$LCI_90_Difference_late_minus_early_trend <- round(bdif$X2_5_,1)
bdif$prob_decreasing_trend <- round(bdif$mean_p,2)

T1$early_trend <- round(T1$mean,1)
T1$LCI_early_trend <- round(T1$X2_5_,1)
T1$UCI_early_trend <- round(T1$X97_5_,1)

T2$late_trend <- round(T2$mean,1)
T2$LCI_late_trend <- round(T2$X2_5_,1)
T2$UCI_late_trend <- round(T2$X97_5_,1)



bdif <- left_join(bdif,T1[,c("species","early_trend","LCI_early_trend","UCI_early_trend")])
bdif <- left_join(bdif,T2[,c("species","late_trend","LCI_late_trend","UCI_late_trend")])


write.csv(bdif[,c("species","Difference_late_minus_early_trend","LCI_90_Difference_late_minus_early_trend","UCI_90_Difference_late_minus_early_trend",
                  "prob_decreasing_trend",
                  "early_trend","LCI_early_trend","UCI_early_trend",
                  "late_trend","LCI_late_trend","UCI_late_trend")],"output/Differences_in_Trends_GAM.csv",row.names = F)

i90 <- filter(indicesraw,year == 2000)
bdiflab <- left_join(i90,bdif,by = "species")
bdiflab$lab = paste(bdiflab$Difference_late_minus_early_trend,":",bdiflab$LCI_90_Difference_late_minus_early_trend,"-",bdiflab$UCI_90_Difference_late_minus_early_trend,"p =",bdiflab$prob_decreasing_trend)
bdiflab$index.raw <- bdiflab$index.raw*1.1

npag <- ceiling(nrow(bdif)/9)
pdf("output/changepoint_graphs_gam.pdf",
    height = 8.5,width = 11)
for(ij in 1:npag){
  trajs <- ggplot(data = mu_,aes(x = year,y = mean))+
    geom_pointrange(data = indicesraw,inherit.aes = FALSE,aes(x = year,y = index.raw,ymax = (uci.raw),ymin = (lci.raw)),alpha = 0.2, size = 0.5)+
    geom_ribbon(aes(ymin = X2_5_,ymax = X97_5_,fill = est),alpha = 0.2)+
    geom_line(aes(colour = est))+
    scale_color_viridis_d(aesthetics = c("colour","fill"),direction = 1,end = 0.8,begin = 0.2)+
    geom_text(data = bdiflab,inherit.aes = FALSE,aes(x = year,y = index.raw,label = lab))+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    facet_wrap_paginate(facets = ~species,scales = "free_y",nrow = 3,ncol = 3,page = ij)
  print(trajs)
}
dev.off()

