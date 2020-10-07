

### Alternative analysis to half-life
## fitting abreakpoint model to the annual indices to estimate the difference in trends at a breakpoint ~10-15 years ago.

## fitting the breakpoint model in Stan
library(tidyverse)
library(rstan)
library(shinystan)
library(tidybayes)
library(ggforce)
# Data - annual indices used in Rosenberg et al 2019 ----------------------

spp.indall <- read.csv('Rosenberg et al annual indices of abundance.csv',stringsAsFactors = F) 


# lookup file relating AOU# to species name
spdf=read.fwf(file='species 2011.txt',
              width = c(11,30,31),
              col.names = c("AOU", "English_Common_Name", "Species"),
              colClasses = c("integer","character","character"),
              strip.white = T)


spp.indall$English_Common_Name <- as.character(spp.indall$species)
spp.indall$yr <- spp.indall$year-1969
sp1 = unique(spp.indall$English_Common_Name)
#sp2 = unique(spdf$English_Common_Name)



# 
# newpal <- viridis::viridis_pal(option = "D")
# 
# cl.orig = newpal(6)[6]
# cl.pred = newpal(6)[5]
# cl.qe50 = newpal(6)[4]
# cl.hl30 = newpal(6)[3]
# cl.hl50 = newpal(6)[2]
# cl.hl70 = newpal(6)[1]


sp.states.i <- unique(spp.indall[,c("species","English_Common_Name","firstyear","lastyear")])
# 


 
spp.ind <- spp.indall


spp.ind <- spp.ind %>% 
  mutate(lind = log(index.raw),
         lsd = (log(uci.raw)-log(lci.raw))/(1.96*2))

meanlsd = mean(spp.ind$lsd,na.rm = T)
spp.ind[which(is.na(spp.ind$lsd)),"lsd"] <- meanlsd

sp = "Red Knot"


mod = stan_model("Changepoint3.stan")

jj = 0
for(sp in unique(spp.ind$species)){
  jj = jj+1
  df = filter(spp.ind,species == sp,
              !is.na(lind))
  df <- df[order(df$yr),]
  stan_dat <- list(Y = nrow(df),
                   Yb = 15,#number of years in later-term trend analysis, could be 10, but that's a very short period of time with which to precisely estimate a trend
                   i = df$lind,
                   sdi = df$lsd,
                   sig_mean = 0,
                   sig_sd = 1)
  #t1 = Sys.time()
  fit <- sampling( 
    object = mod,
    data = stan_dat,    
    chains = 3,             
    warmup = 3000,          
    iter = 5000,            
    cores = 3,              
    open_progress = F,
    control = list(adapt_delta = 0.999,
                   max_treedepth = 15),
    verbose = F
    )         
  
  #launch_shinystan(fit)  
  #Sys.time()-t1
    
  # plot(fit,pars = "mu")
  # 
  # plot(fit,pars = "B")
  # plot(fit,pars = c("beta1","beta2"))
  # plot(fit,pars = c("alpha1","alpha2"))
  # pairs(fit,pars = c("beta1","beta2","alpha1","alpha2"))
  # print(fit)
  outt <- as.data.frame(summary(fit,probs = c(0.025,0.05,0.1, 0.25, 0.50, 0.75,0.9,0.95, 0.975))$summary)
  outt$parameter <- row.names(outt)
  
  outt$species = sp
  
  if(jj == 1){
    out = outt
  }else{
    out <- bind_rows(out,outt)
  }
  
}

b_ <- filter(out,parameter %in% c(paste0("B[",1:53,"]")))

mu_ <- filter(out,parameter %in% c(paste0("mu[",1:53,"]")))

beta_early <- filter(out,parameter %in% c(paste0("beta1")))
beta_late <- filter(out,parameter %in% c(paste0("beta2")))

source("utility_functions.R")
b_ <- mutate(b_,yr = jags_dim(dat = b_, cl = "parameter",var = "B"),
             est = "B")
mu_ <- mutate(mu_,yr = jags_dim(dat = mu_, cl = "parameter",var = "mu"),
              est = "mu")

b_mu = bind_rows(b_,mu_)

fyr_sp <- unique(spp.ind[,c("species","firstyear")])
b_mu <- left_join(b_mu,fyr_sp,by = "species")
b_mu$year <- b_mu$yr+(b_mu$firstyear-1)

names(b_mu)[which(grepl(names(b_mu),pattern = "%",fixed = T))] <- paste0("CL",gsub(names(b_mu)[which(grepl(names(b_mu),pattern = "%",fixed = T))] ,pattern = "%",replacement = "", fixed = TRUE))


betadif <- filter(out,parameter %in% c("betadif"))
p_betaneg <- filter(out,parameter %in% c("beta_neg"))
names(p_betaneg)[which(grepl(names(p_betaneg),pattern = "%",fixed = T))] <- paste0("CL",gsub(names(p_betaneg)[which(grepl(names(p_betaneg),pattern = "%",fixed = T))] ,pattern = "%",replacement = "", fixed = TRUE))
names(betadif)[which(grepl(names(betadif),pattern = "%",fixed = T))] <- paste0("CL",gsub(names(betadif)[which(grepl(names(betadif),pattern = "%",fixed = T))] ,pattern = "%",replacement = "", fixed = TRUE))

names(beta_early)[which(grepl(names(beta_early),pattern = "%",fixed = T))] <- paste0("CL",gsub(names(beta_early)[which(grepl(names(beta_early),pattern = "%",fixed = T))] ,pattern = "%",replacement = "", fixed = TRUE))
names(beta_late)[which(grepl(names(beta_late),pattern = "%",fixed = T))] <- paste0("CL",gsub(names(beta_late)[which(grepl(names(beta_late),pattern = "%",fixed = T))] ,pattern = "%",replacement = "", fixed = TRUE))

prob_annot <- select(p_betaneg,mean,species)
bdif <- select(betadif,CL50,CL5,CL95,species)
bdif <- left_join(bdif,prob_annot,by = "species")
bdif$Difference_late_minus_early_trend <- round((exp(bdif$CL50)-1)*100,1)
bdif$UCI_90_Difference_late_minus_early_trend <- round((exp(bdif$CL95)-1)*100,1)
bdif$LCI_90_Difference_late_minus_early_trend <- round((exp(bdif$CL5)-1)*100,1)
bdif$prob_decreasing_trend <- round(bdif$mean,2)

beta_early$early_trend <- round((exp(beta_early$CL50)-1)*100,1)
beta_early$LCI_early_trend <- round((exp(beta_early$CL5)-1)*100,1)
beta_early$UCI_early_trend <- round((exp(beta_early$CL95)-1)*100,1)

beta_late$late_trend <- round((exp(beta_late$CL50)-1)*100,1)
beta_late$LCI_late_trend <- round((exp(beta_late$CL5)-1)*100,1)
beta_late$UCI_late_trend <- round((exp(beta_late$CL95)-1)*100,1)

bdif <- left_join(bdif,beta_early[,c("species","early_trend","LCI_early_trend","UCI_early_trend")])
bdif <- left_join(bdif,beta_late[,c("species","late_trend","LCI_late_trend","UCI_late_trend")])


write.csv(bdif[,c("species","Difference_late_minus_early_trend","LCI_90_Difference_late_minus_early_trend","UCI_90_Difference_late_minus_early_trend",
                  "prob_decreasing_trend",
                  "early_trend","LCI_early_trend","UCI_early_trend",
                  "late_trend","LCI_late_trend","UCI_late_trend")],"output/Differences_in_Trends.csv",row.names = F)

i90 <- filter(spp.ind,year == 1997)
bdif <- left_join(i90,bdif,by = "species")
bdif$lab = paste(round((exp(bdif$CL50)-1)*100,1),":",round((exp(bdif$CL5)-1)*100,1),"-",round((exp(bdif$CL95)-1)*100,1),"p =",round(bdif$mean,2))
bdif$lind <- bdif$lind+0.1
b_mu <- filter(b_mu,est == "B")

npag <- ceiling(nrow(betadif)/9)
pdf("output/changepoint_graphs.pdf",
    height = 8.5,width = 11)
for(ij in 1:npag){
trajs <- ggplot(data = b_mu,aes(x = year,y = mean))+
  geom_pointrange(data = spp.ind,inherit.aes = FALSE,aes(x = year,y = lind,ymax = log(uci.raw),ymin = log(lci.raw)),alpha = 0.2, size = 0.5)+
  geom_ribbon(aes(ymin = CL2.5,ymax = CL97.5,fill = est),alpha = 0.2)+
  geom_line(aes(colour = est))+
  scale_color_viridis_d(aesthetics = c("colour","fill"),direction = -1,end = 0.8,begin = 0.2)+
  geom_text(data = bdif,inherit.aes = FALSE,aes(x = year,y = lind,label = lab))+
  facet_wrap_paginate(facets = ~species,scales = "free_y",nrow = 3,ncol = 3,page = ij)
print(trajs)
}
dev.off()



