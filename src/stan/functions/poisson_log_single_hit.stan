real poisson_log_single_hit_lpmf(int well_status,
                                 real log_virions){
    real result = 0;
  
  if (well_status == 0) {
    result = poisson_log_lpmf(0 | log_virions);
  } else if (well_status == 1) {
    result = log_diff_exp(0, poisson_log_lpmf(0 | log_virions));
  } else {
    reject("well_status data must be one or zero; given",
           well_status);
  }

  return result;
}
