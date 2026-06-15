
bind_VI_household <- function(list.transmission, time){
  df.hh =NULL
  for(i in 1:nrow(list.transmission$HH)){
    # print(list.transmission$Total.VL.discrete.by.individuals)
    df.hh = rbind(df.hh, data.frame(t = time,
                                    #  VI = list.transmission$VI.discrete.by.individuals[[i]] ,
                                    V = list.transmission$VI.discrete.by.individuals[[i]] ,
                                    Total.VL = list.transmission$Total.VL.discrete.by.individuals[[i]] ,
                                    index =i))
  }
  return(df.hh)
} 

generate_observations_one_household_symptoms <- function(hs, 
                                                         T_P_star = 5,
                                                         T_U_star= 5,
                                                         V_P_star = 4,
                                                         Delta_S_star = 3,
                                                         T_P_eta = 0.4,
                                                         T_U_eta = 0.4,
                                                         V_P_eta = 0.2,
                                                         Delta_S_eta =  0.5,
                                                         h = h,
                                                         obs.time.step = NULL,
                                                         obs.times = NULL, # relative to the symptom onset of the index case
                                                         noise,
                                                         observation.threshold  = detection.limit ,
                                                         SAR.target, 
                                                         SAR.sd = 0,
                                                         offset.symptoms = 1, 
                                                         incubation =TRUE, 
                                                         epid.model,
                                                         # size.dependency,
                                                         size.parameter.alpha){
  
  if(is.null(obs.time.step) & is.null(obs.times)){
    print("Error: provide observation times")
  }
  if(!is.null(obs.time.step) & !is.null(obs.times)){
    print("Error: provide observation times with either obs.time.step or obs.times ")
  }
  if(!is.null(obs.time.step)){
    observation.times = seq(0,Tmax,by=obs.time.step)
  }
  if(!is.null(obs.times)){
    observation.times = obs.times
  } 
  
  SAR.target.sd=c()
  for(i in 1:hs){
    SAR.target.sd[i] = max(0,min(rnorm(n = 1, mean = SAR.target, sd = SAR.sd),1))
  }
  
  #print(SAR.target.sd)
  if(model.virus == "triangle" &  epid.model == "exp"){
    Integrale =10^(2*h)/(h*V_P_star*log(10))*(T_P_star+T_U_star)*(10^(h*(V_P_star))-1) #+  10^(2*h)*(T_U_star/V_P_star/(h*log(10)))*(1-10^(-h*V_P_star))
    m =-log(1-SAR.target.sd) /Integrale
  }
  if(model.virus == "triangle" &  epid.model == "constant"){
    m = -log(1-SAR.target.sd) /(T_P_star+T_U_star)
  } 
  
  max.incubation.period  = 30 # To make sure we don't have incubation periods that are too long
  while(max.incubation.period>=30){
    param.list =  generate_parameters_triangle(T_P_star, T_U_star, V_P_star, Delta_S_star,
                                               T_P_eta, T_U_eta, V_P_eta, Delta_S_eta, 
                                               m, 
                                               h = h,
                                               total.samples = hs,
                                               epid.model = epid.model)
    max.incubation.period = max(param.list$incubation)
  }
  
  HH.parameters  = make_household(param.list, household.size =  hs)
  HH.parameters$m = m 
  offset.t=0
  if(incubation  == TRUE ){
    offset.t = 30
  }
  
  #  list.transmission <- generate_transmission_continuous_time(HH.parameters, size.dependency = size.dependency, Tmax= Tmax+offset.t, discrete.dt = discrete.dt, epid.model = epid.model)
  list.transmission <- generate_transmission_continuous_time(HH.parameters, size.parameter.alpha = size.parameter.alpha, Tmax= Tmax+offset.t, discrete.dt = discrete.dt, epid.model = epid.model)
  df.hh =bind_VI_household(list.transmission, time = seq(0,Tmax+offset.t,by=discrete.dt))
  
  infected.individuals = list.transmission$infected.individuals
  symptom.time = list.transmission$infection.times$symptom.time[infected.individuals] 
  discrete.times.long =  seq(0,Tmax+offset.t,by=discrete.dt)
  infection.time = list.transmission$infection.times$infection.time
  infection.events = list.transmission$infection.events
  
  # get a larger time interval to subset an interval of length Tmax (that starts after the first symptoms)
  if(incubation == TRUE){
    t.start = min(symptom.time)+offset.symptoms # We start the study offset.symptoms days after the first symptoms
  }  
  if(incubation == FALSE){
    t.start=0
  }
  index.discrete.start  = which.min(abs(discrete.times-t.start ))
  
  obs = df.hh %>% 
    filter(t %in% discrete.times.long) %>%
    filter(t  >= discrete.times[index.discrete.start]  ) %>%
    mutate(t = t-min(t)  ) %>%
    filter(t<=Tmax) %>% 
    #filter(t<Tmax) %>% 
    mutate(id= (index)) %>% 
    group_by(id) %>% 
    mutate(t = round(t,digits=5)) %>%
    mutate(observed = ifelse( t %in% observation.times,Total.VL,NA)) %>%
    filter(!is.na(observed)) %>% 
    mutate(observed = Total.VL) %>%
    rename(time = t) %>%
    #dplyr::select(c(time, id, VI,Total.VL)) %>%
    dplyr::select(c(time, id, V ,Total.VL)) %>%
    filter(id %in% infected.individuals) %>%
    mutate(index.case  = id==1) %>%
    mutate(household.size  = hs) %>% 
    mutate(t.infection = ifelse( is.na(infection.time[id]), 0, infection.time[id] )) %>%
    mutate(logV = log10(Total.VL)) %>% # VI+VNI
    rowwise() %>% 
    mutate(observation = ifelse(logV>detection.limit+0.0001,logV + noise*rnorm(1), logV)) %>%
    mutate(observation = ifelse(observation >= observation.threshold, observation, observation.threshold ) ) %>%
    mutate(observation = ifelse(t.infection>=0, observation, observation.threshold ) ) %>%
    mutate(censored = ifelse(observation > observation.threshold, 0,1))  
  
  
  LVLObs = obs$observation 
  
  nt = length(LVLObs)
  start = (1:nt)[!duplicated(obs$id)]
  end = c(start[-1] - 1, nt)
  
  n_infected = length(list.transmission$infected.individuals)
  n_subjects = hs 
  n_non_infected = n_subjects - n_infected 
  
  inf.time = infection.time[infected.individuals] 
  inf.time[1] = 0
  inf.time  = inf.time - t.start
  symptom.time = list.transmission$infection.times$symptom.time[infected.individuals] -t.start
  HH.parameters.output = HH.parameters[infected.individuals,]
  HH.parameters.output$t.infection = inf.time
  
  
  time_common =unique(obs$time) 
  #time_common = seq(0,Tmax,by=obs.time.step)
  L = list(HH.parameters = HH.parameters.output,
           nt = nt,
           start = start,
           end = end, 
           n_subjects = n_subjects, 
           n_infected = n_infected, 
           n_non_infected = n_non_infected,
           time = obs$time,
           time_common = time_common,
           LVLObs = LVLObs,
           censored = obs$censored,
           h = HH.parameters.output$h[1],
           m_input = m,
           infection.time = inf.time,
           symptom.time = symptom.time,
           discrete_times = discrete.times,
           n_discrete_times = n.discrete.times,
           offset.symptoms = offset.symptoms,
           infection.events=infection.events,
           epid.model=epid.model)
  
  return(L)
} 

generate_observations<- function(number.hh,
                                 size.hh, 
                                 T_P_star = 5,
                                 T_U_star= 5,
                                 V_P_star = 4,
                                 Delta_S_star = 3,
                                 T_P_eta = 0.4,
                                 T_U_eta = 0.4,
                                 V_P_eta = 0.2,
                                 Delta_S_eta =  0.5,
                                 h = 0.5,
                                 obs.time.step = NULL,
                                 obs.times = NULL,
                                 noise  = 0, 
                                 SAR.target = 0.4,
                                 SAR.sd = 0,
                                 #offset.symptoms = 1, 
                                 min.offset.symptoms = 1,
                                 max.offset.symptoms = 1,
                                 incubation = TRUE, 
                                 infector.infected = FALSE, 
                                 epid.model = "exp",
                                 #                                 size.dependency = "density",
                                 size.parameter.alpha = 0){ 
  
  
  if(is.null(obs.time.step) & is.null(obs.times)){
    print("Error: provide observation times")
  }
  if(!is.null(obs.time.step) & !is.null(obs.times)){
    print("Error: provide observation times with either obs.time.step or obs.times ")
  }
  if(!is.null(obs.time.step)){
    observation.times = seq(0,Tmax,by=obs.time.step)
  }
  if(!is.null(obs.times)){
    observation.times = obs.times
  } 
  
  
  # if incubation == FALSE, start the measures at the first infection 
  # if incubation == TRUE, start the measures after the symptoms
  
  DF = infection.times = VI.input = infected.individuals =individual.viral.loads = individual.viral.loads.shifted = infection.events = NULL
  n_households = sum(number.hh)
  size_households =  n_infected = n_non_infected =  NT = rep(0,n_households)
  print(number.hh)
  max_infected = max(size.hh)
  NT_max = max_infected*length(observation.times)
  
  start = end = array(data = 0, dim = c(n_households, max_infected))
  time = LVLObs =  censored = array(data = 0, dim = c(n_households, NT_max))
  
  L=0
  if(length(number.hh)!= length(size.hh)){
    return(DF)
  }else{
    I=0
    ma = 0
    all.HH = c()
    for(i in 1:length(number.hh)){
      print(i)
      for(j in 1:number.hh[i]){
        I=I+1
        
        offset.symptoms = runif(n = 1, min = min.offset.symptoms, max = max.offset.symptoms)
        if(infector.infected == TRUE){ 
          L$n_infected=0
          while(L$n_infected <2){
            L = generate_observations_one_household_symptoms(hs = size.hh[i], 
                                                             T_P_star,
                                                             T_U_star,
                                                             V_P_star,
                                                             Delta_S_star,
                                                             T_P_eta,
                                                             T_U_eta,
                                                             V_P_eta,
                                                             Delta_S_eta,
                                                             h=h,
                                                             obs.time.step = obs.time.step,
                                                             obs.times = obs.times,
                                                             noise,
                                                             SAR.target = SAR.target, 
                                                             SAR.sd = SAR.sd,
                                                             offset.symptoms = offset.symptoms,
                                                             incubation = incubation, 
                                                             epid.model = epid.model,
                                                             # size.dependency=size.dependency,
                                                             size.parameter.alpha = size.parameter.alpha)
          }
        }
        
        
        continue.loop = TRUE
        if(infector.infected == FALSE & continue.loop == TRUE){
          L = generate_observations_one_household_symptoms(hs = size.hh[i], 
                                                           T_P_star,
                                                           T_U_star,
                                                           V_P_star,
                                                           Delta_S_star,
                                                           T_P_eta,
                                                           T_U_eta,
                                                           V_P_eta,
                                                           Delta_S_eta,
                                                           h=h,
                                                           obs.time.step = obs.time.step,
                                                           obs.times = obs.times,
                                                           noise,
                                                           SAR.target = SAR.target, 
                                                           SAR.sd = SAR.sd,
                                                           offset.symptoms = offset.symptoms,
                                                           incubation = incubation,
                                                           epid.model = epid.model,
                                                           #size.dependency=size.dependency,
                                                           size.parameter.alpha = size.parameter.alpha)
          
          continue.loop = FALSE
          if(max(L$LVLObs) == detection.limit){
            continue.loop = TRUE
          }
        }
        
        size_households[I] = size.hh[i]
        n_infected[I] = L$n_infected
        n_non_infected[I] = L$n_non_infected
        NT[I] = L$nt
        maxNT = max(NT)
        
        start[I,] = complete_vector(L$start,max_infected)
        end[I,] =  complete_vector(L$end,max_infected)
        time[I,]  =  complete_vector(L$time,NT_max)
        LVLObs[I, ] = complete_vector(L$LVLObs,NT_max)
        censored[I, ] = complete_vector(L$censored,NT_max)
        
        # start_inf_times[I,] = seq(1,max_infected*L$n_discrete_times, by =  L$n_discrete_times)
        # end_inf_times[I,] = start_inf_times[I,]+L$n_discrete_time-1
        L$HH.parameters$t.infection=L$infection.time
        L$HH.parameters$t.symptom=L$symptom.time
        L$HH.parameters$household.id = I
        infection.events[[I]] = L$infection.events
        
        all.HH =rbind(all.HH,L$HH.parameters )
      }
    }
  }
  
  n_households = length(size_households)
  n_time = length(unique(time[1,]) )
  n_inf_times_max = max_infected*L$n_discrete_times
  
  row_infected = as.array(all.HH$household.id)
  n_households = length(size_households)
  C = c()
  for( i in 1:n_households){
    C = c(C, seq(1, n_infected[i]))
  }
  col_infected =  as.array(C )
  
  
  return(list(all.HH = all.HH,
              max_infected = max_infected, 
              start = start, 
              end = end,
              NT_max = NT_max,
              size_households = array(size_households),
              n_infected = array(n_infected),
              n_non_infected = array(n_non_infected),
              n_subjects =  sum(size_households), 
              n_households = n_households,
              LVLObs  = LVLObs,
              censored = censored,
              time = time,
              max_time = max(time),
              n_time = n_time,
              obs_time_step = obs.time.step,
              observation_times = obs.times,
              #  T_infection_max = T_infection_max,
              h = L$h,
              m_input= L$m_input,
              infection.events = infection.events,
              n_inf_times_max=n_inf_times_max,
              n_infected_total = sum(n_infected),
              row_infected = row_infected, 
              col_infected = col_infected, 
              noise=noise,
              T_P_star = T_P_star,
              T_U_star= T_U_star,
              V_P_star =V_P_star,
              Delta_S_star = Delta_S_star,
              T_P_eta = T_P_eta,
              T_U_eta = T_U_eta,
              V_P_eta = V_P_eta,
              Delta_S_eta =  Delta_S_eta,
              h = h,
              epid.model= epid.model,
              size.parameter.alpha = size.parameter.alpha))
}
