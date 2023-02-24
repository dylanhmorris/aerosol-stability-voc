functions {
#include /functions/poisson_log_single_hit.stan
#include /functions/log_hit_rate.stan
#include /functions/recenter.stan
}


data {
// observed quantities
  int<lower = 1> n_total_datapoints;
  int<lower = 0, upper = n_total_datapoints> n_used_datapoints;
  int<lower = 0, upper = 1> well_status[n_total_datapoints];
  vector[n_total_datapoints] dilution;

  int<lower = 1> n_samples;
  int<lower = 1> n_experiments;

  int<lower = 1, upper = n_samples> sample_id[n_total_datapoints];
  int<lower = 1, upper = n_samples> sample_experiment_id[n_samples];
  vector<lower = 0>[n_samples] sample_times;

  real intercept_prior_mean;
  real<lower = 0> intercept_prior_sd;
  real sd_intercept_prior_mode;
  real<lower = 0> sd_intercept_prior_sd;

  real log_hl_prior_mean;
  real<lower = 0> log_hl_prior_sd;

  // flags
  int<lower = 0, upper = 1> debug;

}

transformed data {}

parameters{

// this model uses non-mechanistic halflives
// with no transient phase
  vector[n_experiments] log_half_life;    

// this model uses data in which each point is
// an observation with a potentially different intercept
  vector[n_experiments] mean_intercept;
  vector[n_samples] error_intercept;
  vector<lower = 0>[n_experiments] sd_intercept;

}

transformed parameters {
// declarations
  vector[n_samples] sampled_titer;
  vector[n_samples] intercept;
  vector[n_experiments] decay_rate;

  // calculations
  decay_rate = log10(2) ./ exp(log_half_life);

  for(i_sample in 1:n_samples){
    int i_exp = sample_experiment_id[i_sample];
    
    intercept[i_sample] = recenter(mean_intercept[i_exp],
                                   sd_intercept[i_exp],
                                   error_intercept[i_sample]);

    sampled_titer[i_sample] = intercept[i_sample] - 
      decay_rate[i_exp] * sample_times[i_sample];
  }
}

model {

  // observation process: poisson single hit

for (i_dat in 1:n_used_datapoints) {

  int i_sample = sample_id[i_dat];

  real log_dose = log_hit_rate(sampled_titer[i_sample],
                               dilution[i_dat],
                               10);
  well_status[i_dat] ~ poisson_log_single_hit(log_dose);
  
 }

  // priors
  log_half_life ~ normal(log_hl_prior_mean,
                         log_hl_prior_sd);
  mean_intercept ~ normal(intercept_prior_mean,
                          intercept_prior_sd);
  

  error_intercept ~ normal(0, 1);
  
  sd_intercept  ~ normal(sd_intercept_prior_mode,
                         sd_intercept_prior_sd);
}

generated quantities {
  vector[n_samples] intercept_pred;
  vector[n_samples] sampled_titers_pred;

  for(i_sample in 1:n_samples){
    int i_exp = sample_experiment_id[i_sample];
    
    intercept_pred[i_sample] =
      normal_rng(mean_intercept[i_exp],
                 sd_intercept[i_exp]);
    
    sampled_titers_pred[i_sample] = intercept_pred[i_sample] - 
      decay_rate[i_exp] * sample_times[i_sample];
  }
}
