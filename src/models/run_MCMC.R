#setwd("/Users/nathanaelhoze/Dropbox/Scientific_Projects/Projet_Orchestra/CT_Model/src")
#rm(list = ls())
# adjust to your settings
.libPaths("~/Rlib")

rstan_model_exp = stan_model("models/Stan_scripts/Exp_Model.stan") ###   inference of the infection time.
rstan_model_symptoms_gamma = stan_model("models/Stan_scripts/Symptoms_gamma.stan") ### Infer one value for the incubation period 
rstan_model_constant = stan_model("models/Stan_scripts/Constant_Model.stan") 
rstan_model_positivity_3 = stan_model("models/Stan_scripts/Positivity_Model_3.stan") 

run_inference <- function(xdata, init, rstan_model, iter_warmup = 2500){
  
  fit=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init,control = list(adapt_delta= 0.8))
  Chains = rstan::extract(fit)
  
  init_fixed = function(){
    list(logm = mean(Chains$logm),
         log_T_P_mu = mean(Chains$log_T_P_mu),
         log_T_U_mu =  mean(Chains$log_T_U_mu),
         log_V_P_mu = mean(Chains$log_V_P_mu),        
         log_T_P_sd =  mean(Chains$log_T_P_sd),
         log_T_U_sd = mean(Chains$log_T_U_sd),
         log_V_P_sd =  mean(Chains$log_V_P_sd),
         T_P_eta = as.array( colMeans(Chains$T_P_eta)),
         T_U_eta = as.array(colMeans(Chains$T_U_eta)),
         V_P_eta = as.array( colMeans(Chains$V_P_eta)),
         TI_tmp_positive = as.array( colMeans(Chains$TI_tmp_positive)),
         alpha =  mean(Chains$alpha),
         log_Delta_S_mu =  mean(Chains$log_Delta_S_mu),
         Delta_S_sd =  mean(Chains$Delta_S_sd),
         TI =  as.array(colMeans(Chains$TI)),
         sigma_noise = mean(Chains$sigma_noise))
  }
  fit_fixed=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init_fixed, control = list(adapt_delta= 0.8),algorithm = "Fixed_param")
  Chains_fixed= rstan::extract(fit_fixed)
  
  barD =-2*mean(Chains$loglik)
  Dbar =-2*mean(Chains_fixed$loglik)
  pD = barD-Dbar
  DIC = pD +barD
  
  return(list(Chains = Chains, pD = pD, DIC = DIC))
}

run_inference_symptoms <- function(xdata, init, rstan_model, iter_warmup = 2500){
  
  fit=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init,control = list(adapt_delta= 0.8))
  Chains = rstan::extract(fit)
  
  init_fixed = function(){
    list(logm = mean(Chains$logm),
         alpha =  mean(Chains$alpha),
         gamma_infection_shape = mean(Chains$gamma_infection_shape),
         incubation = mean(Chains$incubation))
  } 
  fit_fixed=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init_fixed, control = list(adapt_delta= 0.8),algorithm = "Fixed_param")
  Chains_fixed= rstan::extract(fit_fixed)
  
  barD =-2*mean(Chains$loglik)
  Dbar =-2*mean(Chains_fixed$loglik)
  pD = barD-Dbar
  DIC = pD +barD
  
  return(list(Chains = Chains, pD = pD, DIC = DIC))
}

## HERE HERE HERE

run_inference_positivity <- function(xdata, init, rstan_model, iter_warmup = 2500){
  
  fit=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init,control = list(adapt_delta= 0.8))
  Chains = rstan::extract(fit)
  
  init_fixed = function(){
    list(logm = mean(Chains$logm),
         alpha =  mean(Chains$alpha),
         gamma_infection_shape = mean(Chains$gamma_infection_shape),
         gamma_detection_shape = mean(Chains$gamma_detection_shape),
         gamma_detection_rate = mean(Chains$gamma_detection_rate),
         #log_gamma_infection_shape_mu = mean(Chains$log_gamma_infection_shape_mu),
         #gamma_infection_shape_eta = as.array( colMeans(Chains$gamma_infection_shape_eta)),
         incubation = mean(Chains$incubation))
  } 
  fit_fixed=sampling(rstan_model, data = xdata,  chains =1, cores=1, iter=iter_warmup*2,warmup=iter_warmup, init=init_fixed, control = list(adapt_delta= 0.8),algorithm = "Fixed_param")
  Chains_fixed= rstan::extract(fit_fixed)
  
  barD =-2*mean(Chains$loglik)
  Dbar =-2*mean(Chains_fixed$loglik)
  pD = barD-Dbar
  DIC = pD +barD
  
  return(list(Chains = Chains, pD = pD, DIC = DIC))
}
