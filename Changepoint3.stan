//
// This Stan program defines a simple model, of a breakpoint
// regression for the continental population trajectory of a bird species
// The breakpoint is fixed at 12-years before the end of the time-series
// 

// The input data are:
// Y = number of years in the time-series (usually ~ 48 years, 1970 - 2017)
// Yb = 12, Fixed breakpoint at Yb years before end of time-series, usually 2005 (2017-Yb)
// vector 'i' of length Y = log-transformed, annual indices of abundance in each year, 
  // estimated with uncertainty from some monitoring program - waterfowl surveys, North American Breeding Bird Survey, Christmas Bird Count, etc.
// vector 'sdi' of length Y = standard deviation of log-transformed annual indices of abundance
  // this 'sdi' value supplies the information on the uncertainty in the species status in a given year
  
data {
  int<lower=0> Y; //number of years
  int<lower=0> Yb; // fixed breakpoint = Y-Yb, 
  vector[Y] i; //log transformed annual indices
  vector[Y] sdi; //log scale sd of the annual indices
  real sig_mean; //prior mean on sigma of regression
  real sig_sd; //prior sd for sigma of regression
  
}


parameters {
     real mu[Y];
  real beta1;
  real beta2;
  real alpha1;
  real<lower=0> sigma;
}

transformed parameters {
   real betadif;
  real beta_neg;
  real B[Y];
      
      for (n in 1:Y) {
          B[n] = n < (Y-Yb) ? alpha1+beta1*(n-(Y-Yb)) : alpha1+(beta2)*(n-(Y-Yb));
      }
      betadif = beta2-beta1; //difference between the slopes
      beta_neg = step(-1*betadif);// = 1 if difference is negative, and = 0 if difference is positive
}



model {
 
  i ~ normal(mu, sdi);
  //the data-estimation model estimating the vector of annual indices mu, 
  //using i (data), and accounting for the uncertainty in i using sdi (data)
  beta1 ~ normal(0,0.1); // early slope 1970 - 2005, 0.1 = ~10%/year which is rather extreme, so prior has a VERY mild regularizing effect
  beta2 ~ normal(0,0.1); // late slope 2005 - 2017
  sigma ~ normal(sig_mean,sig_sd); // variance of the breakpoint regression
  alpha1 ~ normal(0,10); //intercept - index at breakpoint

target += normal_lpdf(B | mu, sigma); //I'll admit, I don't really understand why this works and mu ~ normal(B,sigma) doesn't
}





