functions {
#include /functions/poisson_log_single_hit.stan
#include /functions/log_hit_rate.stan
}


data {
  int<lower = 1> n_total_datapoints;
  int<lower = 0, upper = n_total_datapoints> n_used_datapoints;

  int<lower = 1> n_samples;

  int<lower = 0, upper = 1> well_status[n_total_datapoints];
  vector[n_total_datapoints] dilution;

  int<lower = 1, upper = n_samples> sample_id[n_total_datapoints];

  int<lower = 0, upper = 1> debug;

  // priors
  real titer_prior_mean;
  real<lower = 0> titer_prior_sd;

}

transformed data {}


parameters{

  vector[n_samples] sampled_titer;

}

transformed parameters {
}

model {

  // observation process

  for (i_dat in 1:n_used_datapoints) {
    int i_sample = sample_id[i_dat];
    
    real log_dose = log_hit_rate(sampled_titer[i_sample],
                                 dilution[i_dat],
                                 10);
      
    well_status[i_dat] ~ poisson_log_single_hit(log_dose);
  }

  // priors
  sampled_titer ~ normal(titer_prior_mean,
                         titer_prior_sd);

}


generated quantities {}
