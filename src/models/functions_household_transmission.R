generate_household <-function(vector.indices){
  
  HH.parameters  = param.list[vector.indices,]
  HH.parameters$status = 'contact'
  HH.parameters$status[1] = 'index' 
  
  HH.parameters$t.infection = -1
  HH.parameters$t.infection[1] = 0
  
  return(HH.parameters)
}
 
generate_transmission_continuous_time <- function(HH,
                                                  #size.dependency,
                                                  size.parameter.alpha,
                                                  Tmax = 70, 
                                                  discrete.dt = 0.05, 
                                                  epid.model ){
  
  HH.size = nrow(HH)
  transmission.factor.size = 1/(HH.size-1)^size.parameter.alpha

  complete.list = seq(1, HH.size)
  
  discrete.times = seq(0,Tmax, by = discrete.dt) 
  n.discrete.times = length(discrete.times)
  
  total.infection.risk = NULL
  VI.by.individuals = NULL
  Total.VL.by.individuals = NULL
  
  index.case = HH %>% filter(status == 'index')
  #epid.model = index.case$epid.model
  #  params
  VP = max(HH$V_P)
  h = max(HH$h) 
  if(epid.model == "exp"){
    DELTA_T = min(discrete.dt, 1 /(max(HH$m)*10^(h*(2+VP)))/nrow(HH))
    DELTA_T = max(DELTA_T, 0.005)
    
  }
  if(epid.model == "constant"){
    DELTA_T = 0.05
  }
  
  # print(DELTA_T)
  
  #  print(DELTA_T)
  #  DELTA_T = max(DELTA_T, 0.005)
  total.time <- seq(0, Tmax, by = DELTA_T)
  n.times = length(total.time)
  
  non.infected = HH %>% filter(t.infection == -1)
  n.non.infected = nrow(non.infected)
  
  titers.index.case = get_viral_load(params = index.case , delay= 0, Tmax = Tmax ,dt = DELTA_T)   
  
  if(epid.model == "exp"){
    titers.index.case = titers.index.case  %>% 
      mutate(prob.transmission = transmission.factor.size*DELTA_T*HH$m[1]*(V )^h )
  }
  if(epid.model == "constant"){
    titers.index.case = titers.index.case  %>% 
      mutate(prob.transmission = ifelse(logV>=infection.threshold, transmission.factor.size*DELTA_T*HH$m[1], 0 ))
  }
  
  discrete.values = get_viral_load(params = index.case , delay= 0, Tmax = Tmax, dt = discrete.dt)   
  total.infection.risk[[1]] =  titers.index.case$prob.transmission
  # VI.by.individuals[[1]] = discrete.values$VI
  VI.by.individuals[[1]] = discrete.values$V 
  Total.VL.by.individuals[[1]] = discrete.values$V   #discrete.values$VNI
  
  for(index.non.infected in 2:HH.size){
    total.infection.risk[[index.non.infected]] = rep(0,n.times)
    VI.by.individuals[[index.non.infected]] = rep(100,n.discrete.times)
    Total.VL.by.individuals[[index.non.infected]] = rep(100,n.discrete.times)
  } 
  
  infected.list = c(1) # the only infected is the index case
  non.infected.list = seq(2,HH.size) 
  current.index.time = 0 # initialize somehting ?
  j=0 
  new.infection = TRUE
  viral.load.infection.events = NULL # contains the viral load of infectious individuals and whether they infected someone
  # PP=0 
  for(t.index in 1:n.times){
    
    current.time = total.time[t.index]
    infected.list.temp = infected.list
    for(infector.index in infected.list){
      p.inf = total.infection.risk[[infector.index]][t.index]# the probability of infection viral load at the time of infection
      if(length(non.infected.list>0)){
        for(potential.infected.index in non.infected.list){# potentially infected during this time
          r=runif(n=1)
          if(r<p.inf){
            # there is an infection
            # infected.index is the index of the new infection
            viral.load.infection.events = rbind(viral.load.infection.events,data.frame(time = total.time[t.index], 
                                                                                       infector = infector.index,
                                                                                       infected = potential.infected.index,  
                                                                                       infection = 1))
            infected.parameters = HH[potential.infected.index,]
            
            delay = min(Tmax-DELTA_T, current.time)
            titers.case = get_viral_load(params = infected.parameters, delay= delay, Tmax = Tmax, dt = DELTA_T)
            # titers.case = get_viral_load(params = infected.parameters, delay= delay, Tmax = Tmax, dt = DELTA_T) %>% 
            #  mutate(prob.transmission = transmission.factor.size*DELTA_T*infected.parameters$m*(VI)^h )
            #   mutate(prob.transmission = transmission.factor.size*DELTA_T*infected.parameters$m*(V)^h )
            
            if(epid.model == "exp"){
              titers.case = titers.case  %>% 
                mutate(prob.transmission = transmission.factor.size*DELTA_T*infected.parameters$m*(V )^h )
            }
            if(epid.model == "constant"){
              titers.case = titers.case  %>% 
                mutate(prob.transmission = ifelse(logV>=infection.threshold, transmission.factor.size*DELTA_T*infected.parameters$m, 0 ))
            }
            
            discrete.values = get_viral_load(params = infected.parameters , delay= delay, Tmax = Tmax, dt = discrete.dt)   
            
            total.infection.risk[[potential.infected.index]] =  titers.case$prob.transmission 
            VI.by.individuals[[potential.infected.index]] = discrete.values$V  
            Total.VL.by.individuals[[potential.infected.index]] = discrete.values$V 
            
            infected.list.temp  = c(infected.list.temp, potential.infected.index)
            # remove from the list of infected
          }
          if(r>=p.inf){
            # there is no infection
            viral.load.infection.events = rbind(viral.load.infection.events, 
                                                data.frame(time = current.time, 
                                                           infector = infector.index,
                                                           infected = -1,  
                                                           infection = 0))
          }
        }  
      } 
      infected.list.temp = sort(infected.list.temp)
      non.infected.list = complete.list[-infected.list.temp] 
    }  
    infected.list = infected.list.temp
  }
  
  n.infection = length(infected.list)
  # get the viral load at the time of infection
  infection.events = viral.load.infection.events%>%filter(infection==1)
  infectious.viral.load.at.infection = infection.events$V  
  total.viral.load.at.infection = infection.events$Total.VL
  
  # get the first infection time within the HH
  tinf =  infection.events$time
  tinf = tinf[which(tinf>0)]
  
  if(length(tinf)>0){  
    min.infection =  min(tinf)
    max.infection =  max(tinf)
  }else{
    min.infection = -1
    max.infection = 0
  }
  
  infection.times = data.frame(infected =1:HH.size, infection.time= rep(-1,HH.size), symptom.time = rep(-1,HH.size))
  infection.times$infection.time[1] = NA
  infection.times$symptom.time[1] = HH$incubation[1]
  
  if(n.infection>1){
    infection.times[infection.events$infected,]$infection.time = infection.events$time
    infection.times[infection.events$infected,]$symptom.time = infection.events$time + HH$incubation[infection.events$infected]
  } 
  
  return(list(HH = HH, 
              HH.SAR =(n.infection-1)/(HH.size-1),
              n.infection = n.infection, # Including the index case
              min.t.infection = min.infection,
              max.t.infection = max.infection,
              infectious.viral.load.at.infection=infectious.viral.load.at.infection,
              total.viral.load.at.infection=total.viral.load.at.infection,
              infection.events  = infection.events,
              all.events =viral.load.infection.events,
              infection.times=infection.times,
              VI.discrete.by.individuals = VI.by.individuals, 
              Total.VL.discrete.by.individuals = Total.VL.by.individuals, 
              infected.individuals = infected.list))
}

generate_transmission_short <- function(HH, discrete.dt = 0.05){
  
  HH.size = nrow(HH)
  
  complete.list = seq(1, HH.size)
  
  discrete.times = seq(0,Tmax, by = discrete.dt) 
  n.discrete.times = length(discrete.times)
  
  index.case = HH %>% filter(status == 'index')
  VP = max(HH$V_P)
  
  h = max(HH$h)
  DELTA_T = min(discrete.dt, 0.1/(max(HH$m)*10^(h*(2+VP)))/nrow(HH))
  
  titers.index.case = get_viral_load(params = index.case , delay= 0, Tmax = Tmax ,dt = DELTA_T)    %>% 
    mutate(prob.transmission =  DELTA_T*HH$m[1]*(VI)^h ) %>%
    mutate(cumulative.prob.transmission = cumsum(prob.transmission))
  
  r =runif(n=1)
  if(r >max(titers.index.case$cumulative.prob.transmission)){
    n.infection=1
    t.infection=NA
  }
  if(r <= max(titers.index.case$cumulative.prob.transmission)){
    n.infection=2
    index.infection=which.min(abs(r-titers.index.case$cumulative.prob.transmission))
    t.infection = as.numeric(titers.index.case$t[index.infection])
  }
  
  return(list( n.infection = n.infection, 
               t.infection = t.infection))
}
