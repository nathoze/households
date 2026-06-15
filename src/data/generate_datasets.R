
generate_datasets <-function(n.household, 
                             size.households = NULL,
                             filename,
                             obs.time.step = 3,
                             noise = 0.8,
                             prob_symptomatic = 0.6,
                             T_P_star = 5,
                             T_U_star= 5,
                             V_P_star = 4,
                             Delta_S_star = 3,
                             T_P_eta = 0.4,
                             T_U_eta = 0.4,
                             V_P_eta = 0.2,
                             Delta_S_eta =  0.5,
                             SAR.target = 0.4,
                             SAR.sd = 0, 
                             insee.size = FALSE,
                             epid.model = "exp",
                             #size.dependency="density",
                             size.parameter.alpha = 0.5){
  
  if(insee.size ==TRUE){
    
    sizes = readRDS(file = "../data/processed_data/INSEE_French_HH_size.rds" )
    size.households = sample(sizes,size = n.household)
    nh =size.hh = c()
    for(i in 1: max(size.households)){
      ni =  length(which(size.households ==i )) 
      if(ni>0){
        nh =c(nh, length(which(size.households ==i )))
        size.hh =c(size.hh, i)
      }
    }
  }else{
    nh  = n.household
    size.hh = c(size.households)
  }
  
  parameters = list(number.hh = nh, 
                    size.hh = size.hh,
                    obs.time.step = obs.time.step,
                    noise = noise,
                    prob_symptomatic =prob_symptomatic,
                    T_P_star = T_P_star,
                    T_U_star= T_U_star,
                    V_P_star = V_P_star,
                    Delta_S_star = Delta_S_star,
                    T_P_eta = T_P_eta,
                    T_U_eta = T_U_eta,
                    V_P_eta = V_P_eta,
                    Delta_S_eta =  Delta_S_eta,
                    Delta_S_star =  Delta_S_star,
                    SAR.target = SAR.target,
                    SAR.sd = SAR.sd,
                    epid.model = epid.model,
                    size.parameter.alpha= size.parameter.alpha)
  
  # 
  xdata = generate_observations(number.hh = parameters$number.hh ,
                                size.hh = parameters$size.hh,
                                T_P_star = parameters$T_P_star,
                                T_U_star= parameters$T_U_star,
                                V_P_star = parameters$V_P_star,
                                Delta_S_star = parameters$Delta_S_star,
                                T_P_eta = parameters$T_P_eta,
                                T_U_eta = parameters$T_U_eta,
                                V_P_eta = parameters$V_P_eta,
                                Delta_S_eta =  parameters$Delta_S_eta,
                                obs.time.step = parameters$obs.time.step,
                                noise=parameters$noise,
                                SAR.target = parameters$SAR.target,
                                SAR.sd = parameters$SAR.sd,
                                incubation = TRUE,
                                offset.symptoms = 1,
                                infector.infected = FALSE,
                                epid.model = epid.model,
                                size.parameter.alpha = size.parameter.alpha)
  
  
  
  params = xdata$all.HH
  xdata$sigma_prior[1] = 0.5
  xdata$sigma_prior[2] = 0.1
  xdata$T_P_star = 5
  xdata$T_U_star = 5
  xdata$V_P_star = 2
  xdata$K = 10^(2*xdata$h) # infectious virus
  
  # symptom time
  xdata$t_symptoms = xdata$all.HH$t.symptom
  symptomatic =rep(0, xdata$n_infected_total)
  symptomatic_array = array(0, dim=c(xdata$n_households,xdata$max_infected))
  
  I=0
  for(k in 1:xdata$n_households){
    for(i in 1:xdata$n_infected[k]){
      I=I+1
      symptomatic[I] = runif(n=1)<parameters$prob_symptomatic
      if(i == 1)
        symptomatic[I] = 1
      
      symptomatic_array[k,i] =  symptomatic[I]
      
    }
  }
  xdata$symptomatic = symptomatic
  xdata$symptomatic_array = symptomatic_array
  
  
  K=0
  first.detectable= rep(0, xdata$n_infected_total)
  
  time = seq(0,xdata$max_time,by=xdata$obs_time_step)
  
  
  
  for(i in 1:xdata$n_households){
    for(j in 1:xdata$n_infected[i])  {
      K=K+1
      xdata$LVLObs[i,]
      start= xdata$start[i,j]
      end= xdata$end[i,j]
      first.detectable[K] = time[which(xdata$LVLObs[i,start:end]>detection.limit)[1]]
    }
  }
  xdata$first_detectable_time = first.detectable
  for(k in 1:xdata$n_infected_total){
    if(xdata$symptomatic[k])
      xdata$first_detectable_time[k] = min(first.detectable[k],xdata$t_symptoms[k],na.rm = TRUE)
  }
  
  xdata$all.HH$size.parameter.alpha = size.parameter.alpha
  
  
  # iter_warmup=5000
  # init = function(){
  #   list(logm = runif(n=1, min=-7, max=0),
  #        log_T_P_mu = 0,
  #        log_T_U_mu =  0,
  #        log_V_P_mu = 0,        
  #        log_T_P_sd =  1,
  #        log_T_U_sd = 1,
  #        log_V_P_sd =  1,
  #        T_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
  #        T_U_eta = as.array(rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
  #        V_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
  #        TI =  as.array(runif( n = xdata$n_infected_total,  min = 0, max=5)),
  #        sigma_noise = 0.1)
  # }
  
  
  iter_warmup=2500
  init = function(){
    list(logm = runif(n=1, min=-7, max=0),
         alpha =  runif(n=1, min=0, max=1),
         log_Delta_S_mu =  2,
         Delta_S_sd =  0.5,
         log_T_P_mu = 0,
         log_T_U_mu =  0,
         log_V_P_mu = 0,        
         log_T_P_sd =  1,
         log_T_U_sd = 1,
         log_V_P_sd =  1,
         T_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         T_U_eta = as.array(rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         V_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         #  TI_tmp =  as.array(runif( n = xdata$n_infected_total,  min = -10, max=30)),
         TI_tmp_positive =  as.array(runif( n = xdata$n_infected_total,  min = 0, max=10)),
         sigma_noise = 0.8)
  }
  
  save('xdata',"parameters","init", file = paste0("../data/simulated_data/", filename,"_",epid.model,".rds"))
  
}


generate_datasets_sampling_scheme <-function(n.household, 
                                             size.households = NULL,
                                             filename,
                                             obs.time.step = NULL,
                                             obs.times = c(0,3,7,14,21,28),
                                             noise = 0.8,
                                             prob_symptomatic = 0.6,
                                             T_P_star = 5,
                                             T_U_star= 5,
                                             V_P_star = 4,
                                             Delta_S_star = 3,
                                             T_P_eta = 0.4,
                                             T_U_eta = 0.4,
                                             V_P_eta = 0.2,
                                             Delta_S_eta =  0.5,
                                             SAR.target = 0.4,
                                             SAR.sd = 0, 
                                             insee.size = FALSE,
                                             epid.model = "exp",
                                             #size.dependency="density",
                                             size.parameter.alpha = 0.5){
  
  if(insee.size ==TRUE){
    
    sizes = readRDS(file = "../data/processed_data/INSEE_French_HH_size.rds" )
    size.households = sample(sizes,size = n.household)
    nh =size.hh = c()
    for(i in 1: max(size.households)){
      ni =  length(which(size.households ==i )) 
      if(ni>0){
        nh =c(nh, length(which(size.households ==i )))
        size.hh =c(size.hh, i)
      }
    }
  }else{
    nh  = n.household
    size.hh = c(size.households)
  }
  
  parameters = list(number.hh = nh, 
                    size.hh = size.hh,
                    obs.time.step = obs.time.step,
                    obs.times=obs.times,
                    noise = noise,
                    prob_symptomatic =prob_symptomatic,
                    T_P_star = T_P_star,
                    T_U_star= T_U_star,
                    V_P_star = V_P_star,
                    Delta_S_star = Delta_S_star,
                    T_P_eta = T_P_eta,
                    T_U_eta = T_U_eta,
                    V_P_eta = V_P_eta,
                    Delta_S_eta =  Delta_S_eta,
                    Delta_S_star =  Delta_S_star,
                    SAR.target = SAR.target,
                    SAR.sd = SAR.sd,
                    epid.model = epid.model,
                    size.parameter.alpha= size.parameter.alpha)
  
  # 
  xdata = generate_observations(number.hh = parameters$number.hh ,
                                size.hh = parameters$size.hh,
                                T_P_star = parameters$T_P_star,
                                T_U_star= parameters$T_U_star,
                                V_P_star = parameters$V_P_star,
                                Delta_S_star = parameters$Delta_S_star,
                                T_P_eta = parameters$T_P_eta,
                                T_U_eta = parameters$T_U_eta,
                                V_P_eta = parameters$V_P_eta,
                                Delta_S_eta =  parameters$Delta_S_eta,
                                obs.time.step = parameters$obs.time.step,
                                obs.times = parameters$obs.times,
                                noise=parameters$noise,
                                SAR.target = parameters$SAR.target,
                                SAR.sd = parameters$SAR.sd,
                                incubation = TRUE,
                                offset.symptoms = 1,
                                infector.infected = FALSE,
                                epid.model = epid.model,
                                size.parameter.alpha = size.parameter.alpha)
  
  
  params = xdata$all.HH
  xdata$sigma_prior[1] = 0.5
  xdata$sigma_prior[2] = 0.1
  xdata$T_P_star = 5
  xdata$T_U_star = 5
  xdata$V_P_star = 2
  xdata$K = 10^(2*xdata$h) # infectious virus
  
  # symptom time
  xdata$t_symptoms = xdata$all.HH$t.symptom
  symptomatic =rep(0, xdata$n_infected_total)
  symptomatic_array = array(0, dim=c(xdata$n_households,xdata$max_infected))
  
  I=0
  for(k in 1:xdata$n_households){
    for(i in 1:xdata$n_infected[k]){
      I=I+1
      symptomatic[I] = runif(n=1)<parameters$prob_symptomatic
      if(i == 1)
        symptomatic[I] = 1
      
      symptomatic_array[k,i] =  symptomatic[I]
      
    }
  }
  xdata$symptomatic = symptomatic
  xdata$symptomatic_array = symptomatic_array
  
  K=0
  first.detectable= rep(0, xdata$n_infected_total)
  
  # time = seq(0,xdata$max_time,by=xdata$obs_time_step)
  time = obs.times
  
  for(i in 1:xdata$n_households){
    for(j in 1:xdata$n_infected[i])  {
      K=K+1
      xdata$LVLObs[i,]
      start= xdata$start[i,j]
      end= xdata$end[i,j]
      first.detectable[K] = time[which(xdata$LVLObs[i,start:end]>detection.limit)[1]]
    }
  }
  xdata$first_detectable_time = first.detectable
  for(k in 1:xdata$n_infected_total){
    if(xdata$symptomatic[k])
      xdata$first_detectable_time[k] = min(first.detectable[k],xdata$t_symptoms[k],na.rm = TRUE)
  }
  
  xdata$all.HH$size.parameter.alpha = size.parameter.alpha
  
  iter_warmup=2500
  init = function(){
    list(logm = runif(n=1, min=-7, max=0),
         alpha =  runif(n=1, min=0, max=1),
         log_Delta_S_mu =  2,
         Delta_S_sd =  0.5,
         log_T_P_mu = 0,
         log_T_U_mu =  0,
         log_V_P_mu = 0,        
         log_T_P_sd =  1,
         log_T_U_sd = 1,
         log_V_P_sd =  1,
         T_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         T_U_eta = as.array(rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         V_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         #  TI_tmp =  as.array(runif( n = xdata$n_infected_total,  min = -10, max=30)),
         TI_tmp_positive =  as.array(runif( n = xdata$n_infected_total,  min = 0, max=10)),
         sigma_noise = 0.8)
  }
  
  save('xdata',"parameters","init", file = paste0("../data/simulated_data/", filename,"_",epid.model,".rds"))
  
}


tmp_generate_datasets_sampling_scheme <-function(n.household, 
                                             size.households = NULL,
                                             filename,
                                             obs.time.step = NULL,
                                             obs.times = c(0,3,7,14,21,28),
                                             noise = 0.8,
                                             prob_symptomatic = 0.6,
                                             T_P_star = 5,
                                             T_U_star= 5,
                                             V_P_star = 4,
                                             Delta_S_star = 3,
                                             T_P_eta = 0.4,
                                             T_U_eta = 0.4,
                                             V_P_eta = 0.2,
                                             Delta_S_eta =  0.5,
                                             SAR.target = 0.4,
                                             SAR.sd = 0, 
                                             insee.size = FALSE,
                                             epid.model = "exp",
                                             min.offset.symptoms = 1,
                                             max.offset.symptoms = 1,
                                             #size.dependency="density",
                                             size.parameter.alpha = 0.5){
  
  
  
  
  if(min.offset.symptoms > max.offset.symptoms){
    print("Error : min.offset.symptoms should be lower than max.offset.symptoms")
  }
  
  
  if(insee.size ==TRUE){
    
    sizes = readRDS(file = "../data/processed_data/INSEE_French_HH_size.rds" )
    size.households = sample(sizes,size = n.household)
    nh =size.hh = c()
    for(i in 1: max(size.households)){
      ni =  length(which(size.households ==i )) 
      if(ni>0){
        nh =c(nh, length(which(size.households ==i )))
        size.hh =c(size.hh, i)
      }
    }
  }else{
    nh  = n.household
    size.hh = c(size.households)
  }
  
  parameters = list(number.hh = nh, 
                    size.hh = size.hh,
                    obs.time.step = obs.time.step,
                    obs.times=obs.times,
                    noise = noise,
                    prob_symptomatic =prob_symptomatic,
                    T_P_star = T_P_star,
                    T_U_star= T_U_star,
                    V_P_star = V_P_star,
                    Delta_S_star = Delta_S_star,
                    T_P_eta = T_P_eta,
                    T_U_eta = T_U_eta,
                    V_P_eta = V_P_eta,
                    Delta_S_eta =  Delta_S_eta,
                    Delta_S_star =  Delta_S_star,
                    SAR.target = SAR.target,
                    SAR.sd = SAR.sd,
                    epid.model = epid.model,
                    size.parameter.alpha= size.parameter.alpha)
 
 
  xdata = generate_observations(number.hh = parameters$number.hh ,
                                size.hh = parameters$size.hh,
                                T_P_star = parameters$T_P_star,
                                T_U_star= parameters$T_U_star,
                                V_P_star = parameters$V_P_star,
                                Delta_S_star = parameters$Delta_S_star,
                                T_P_eta = parameters$T_P_eta,
                                T_U_eta = parameters$T_U_eta,
                                V_P_eta = parameters$V_P_eta,
                                Delta_S_eta =  parameters$Delta_S_eta,
                                obs.time.step = parameters$obs.time.step,
                                obs.times = parameters$obs.times,
                                noise=parameters$noise,
                                SAR.target = parameters$SAR.target,
                                SAR.sd = parameters$SAR.sd,
                                incubation = TRUE,
                                min.offset.symptoms = min.offset.symptoms,
                                max.offset.symptoms = max.offset.symptoms,
                                infector.infected = FALSE,
                                epid.model = epid.model,
                                size.parameter.alpha = size.parameter.alpha)
  
  
  params = xdata$all.HH
  xdata$sigma_prior[1] = 0.5
  xdata$sigma_prior[2] = 0.1
  xdata$T_P_star = 5
  xdata$T_U_star = 5
  xdata$V_P_star = 2
  xdata$K = 10^(2*xdata$h) # infectious virus
  
  # symptom time
  xdata$t_symptoms = xdata$all.HH$t.symptom
  symptomatic =rep(0, xdata$n_infected_total)
  symptomatic_array = array(0, dim=c(xdata$n_households,xdata$max_infected))
  
  I=0
  for(k in 1:xdata$n_households){
    for(i in 1:xdata$n_infected[k]){
      I=I+1
      symptomatic[I] = runif(n=1)<parameters$prob_symptomatic
      if(i == 1)
        symptomatic[I] = 1
      
      symptomatic_array[k,i] =  symptomatic[I]
      
    }
  }
  xdata$symptomatic = symptomatic
  xdata$symptomatic_array = symptomatic_array
  
  K=0
  first.detectable= rep(0, xdata$n_infected_total)
  
  # time = seq(0,xdata$max_time,by=xdata$obs_time_step)
  time = obs.times
  
  for(i in 1:xdata$n_households){
    for(j in 1:xdata$n_infected[i])  {
      K=K+1
      xdata$LVLObs[i,]
      start= xdata$start[i,j]
      end= xdata$end[i,j]
      first.detectable[K] = time[which(xdata$LVLObs[i,start:end]>detection.limit)[1]]
    }
  }
  xdata$first_detectable_time = first.detectable
  for(k in 1:xdata$n_infected_total){
    if(xdata$symptomatic[k])
      xdata$first_detectable_time[k] = min(first.detectable[k],xdata$t_symptoms[k],na.rm = TRUE)
  }
  
  xdata$all.HH$size.parameter.alpha = size.parameter.alpha
  
  iter_warmup=2500
  init = function(){
    list(logm = runif(n=1, min=-7, max=0),
         alpha =  runif(n=1, min=0, max=1),
         log_Delta_S_mu =  2,
         Delta_S_sd =  0.5,
         log_T_P_mu = 0,
         log_T_U_mu =  0,
         log_V_P_mu = 0,        
         log_T_P_sd =  1,
         log_T_U_sd = 1,
         log_V_P_sd =  1,
         T_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         T_U_eta = as.array(rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         V_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
         #  TI_tmp =  as.array(runif( n = xdata$n_infected_total,  min = -10, max=30)),
         TI_tmp_positive =  as.array(runif( n = xdata$n_infected_total,  min = 0, max=10)),
         sigma_noise = 0.8)
  }
  
  save('xdata',"parameters","init", file = paste0("../data/simulated_data/", filename))
  
}



