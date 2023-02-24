  /* filename: recenter.stan
   * author: Dylan H. Morris (dylanhmorris.com)
   * 
   * description: header file defining function
   * to undo a non-centered parametrization
   */

real recenter(real param_mean,
              real param_sd,
              real param_error){

  return param_mean + param_sd * param_error;
}
