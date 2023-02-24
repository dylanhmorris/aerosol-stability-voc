real log_hit_rate(real true_log_titer,
                  real log_dilution,
                  real log_base){

  real log_dilute_dose = true_log_titer + log_dilution;
  real log_dilute_virions = log(log(2)) +
    log_dilute_dose * log(log_base);
  return log_dilute_virions;
}
