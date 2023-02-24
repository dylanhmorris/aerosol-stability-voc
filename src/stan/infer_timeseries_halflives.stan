functions {
#include /functions/poisson_log_single_hit.stan
#include /functions/log_hit_rate.stan
#include /functions/recenter.stan
}

data {
// observed quantities
  int<lower = 1> n_total_datapoints;
  int<lower = 0, upper = n_total_datapoints> n_used_datapoints;

  int<lower = 1> n_samples;

  int<lower = 0, upper = 1> well_status[n_total_datapoints];
  vector[n_total_datapoints] dilution;

  int<lower = 1, upper = n_samples> sample_id[n_total_datapoints];

  int<lower = 1> n_runs;
  int<lower = 1, upper = n_runs> sample_run_id[n_samples];
  vector[n_samples] log_change_virus_genomes;

  int<lower = 1> n_experiments;
  int<lower = 1, upper = n_experiments> experiment_id[n_total_datapoints];
  int<lower = 1, upper = n_samples> sample_experiment_id[n_samples];
  vector<lower = 0>[n_samples] sample_times;
  int<lower = 1, upper = n_experiments> run_experiment_id[n_runs];


  // parameters for priors  
  real log_hl_prior_mean;
  real<lower = 0> log_hl_prior_sd;

  real intercept_prior_mean;
  real<lower = 0> intercept_prior_sd;
  real sd_intercept_prior_mode;
  real<lower = 0> sd_intercept_prior_sd;
  
// flags
  int<lower = 0, upper = 1> debug;
}

transformed data {}

parameters{

  vector[n_experiments] log_half_life;    

  // this model uses data in which each run
  // has a single true intercept
  vector[n_experiments] experiment_mean_intercept;
  vector<lower = 0>[n_experiments] experiment_sd_intercept;
  
  vector[n_runs] run_intercept_errors;
}

transformed parameters {
// declarations
  vector[n_runs] run_intercept;
  vector[n_samples] sampled_titer;
  vector[n_samples] intercept;
  vector[n_experiments] decay_rate;

  decay_rate = log10(2) ./ exp(log_half_life);
  
  // calculations

  for(i_run in 1:n_runs){
    int exp_id = run_experiment_id[i_run];
    run_intercept[i_run] = recenter(experiment_mean_intercept[exp_id],
                                    experiment_sd_intercept[exp_id],
                                    run_intercept_errors[i_run]);
  }


  for(i_sample in 1:n_samples){
    int exp_id = sample_experiment_id[i_sample];
    int run_id = sample_run_id[i_sample];
    real t = sample_times[i_sample];
    real ith_intercept = run_intercept[run_id];
    real ith_predicted_titer =
      ith_intercept
      - decay_rate[exp_id] * t
      + log_change_virus_genomes[i_sample];

    // save values
    intercept[i_sample] = ith_intercept;
    sampled_titer[i_sample] = ith_predicted_titer;
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


  experiment_mean_intercept ~ normal(intercept_prior_mean,
                                     intercept_prior_sd);
  experiment_sd_intercept ~ normal(sd_intercept_prior_mode,
                                   sd_intercept_prior_sd);
  
  run_intercept_errors ~ std_normal();
      
}

// predictive checks
generated quantities {
  vector[n_samples] sampled_titers_pred;
  vector[n_runs] intercept_pred;

  for(r_id in 1:n_runs){
    int i_exp = run_experiment_id[r_id];
    intercept_pred[r_id] = normal_rng(experiment_mean_intercept[i_exp],
                                      experiment_sd_intercept[i_exp]);
  }

  for(s_id in 1:n_samples){
    int exp_id = sample_experiment_id[s_id];
    int run_id = sample_run_id[s_id];
    real t = sample_times[s_id];
    real ith_intercept = intercept_pred[run_id];
    sampled_titers_pred[s_id] =
      ith_intercept
      - decay_rate[exp_id] * t
      + log_change_virus_genomes[s_id];
  }
}
