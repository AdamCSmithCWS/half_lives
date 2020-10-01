//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> Y; //number of years

  vector[Y] i; //log transformed annual indices
  vector[Y] sdi; //log scale sd of the annual indices
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
     real mu[Y];

  real beta1;
  real beta2;
  real alpha1;
  real alpha2;
  real<lower=0> sigma;
}

transformed parameters {
   real betadif;
  real B[Y];
      
      for (n in 1:Y) {
          B[n] = n < (Y-12) ? alpha1+beta1*n : alpha2+(beta2)*(n-(Y-13));
      }
      betadif = beta2-beta1;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
 
  i ~ normal(mu, sdi);
  //the data-estimation model estimating the vector of annual indices accounting for the uncertainty
  beta1 ~ normal(0,0.1); // early slope
  beta2 ~ normal(0,0.1); // late slope
  sigma ~ normal(0,0.5);
  alpha1 ~ normal(0,10); //intercept
  alpha2 ~ normal(0,10); //intercept


target += normal_lpdf(mu | B, sigma);
}





