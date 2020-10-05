

### Alternative analysis to half-life
## fitting abreakpoint model to the annual indices to estimate the difference in trends at a breakpoint ~10-15 years ago.

## fitting the breakpoint model in Stan
library(tidyverse)
library(rstan)
library(shinystan)
library(tidybayes)
library(ggforce)
# Data - annual indices used in Rosenberg et al 2019 ----------------------
sp = "Barn Swallow"

load(file = "Barn_Swallow_Example_Data.RData")

mod = stan_model("Changepoint3.stan")

  stan_dat <- list(Y = nrow(df),
                   Yb = 12,
                   i = df$lind,
                   sdi = df$lsd)

  fit <- sampling( 
    object = mod,
    data = stan_dat,    
    chains = 3,             
    warmup = 3000,          
    iter = 5000,            
    cores = 3,              
    open_progress = F,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 15),
    verbose = F
    )         
  
  launch_shinystan(fit)  
 



