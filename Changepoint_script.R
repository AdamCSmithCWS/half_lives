

### Alternative analysis to half-life
## fitting abreakpoint model to the annual indices to estimate the difference in trends at a breakpoint ~10-15 years ago.

## fitting the breakpoint model in Stan
library(tidyverse)
library(rstan)
library(shinystan)
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


stan_dat <- list(Y = nrow(df))
sp = "Red Knot"


mod = stan_model("Changepoint3.stan")

for(sp in unique(spp.ind$species)[1:10]){
  df = filter(spp.ind,species == sp,
              !is.na(lind))
  df <- df[order(df$yr),]
  stan_dat <- list(Y = nrow(df),
                   #zeros = c(0,0,0),
                   i = df$lind,
                   sdi = df$lsd)
  #t1 = Sys.time()
  fit <- sampling( 
    object = mod,
    data = stan_dat,    
    chains = 3,             
    warmup = 2000,          
    iter = 4000,            
    cores = 3,              
    open_progress = F,
    control = list(adapt_delta = 0.999,
                   max_treedepth = 15),
    verbose = F
    )         
  
  launch_shinystan(fit)  
  #Sys.time()-t1
    
  # plot(fit,pars = "mu")
  # 
  # plot(fit,pars = "B")
  # plot(fit,pars = c("beta1","beta2"))
  # plot(fit,pars = c("alpha1","alpha2"))
  # pairs(fit,pars = c("beta1","beta2","alpha1","alpha2"))
  # print(fit)
  
  extract(fit,parm)
  
}








