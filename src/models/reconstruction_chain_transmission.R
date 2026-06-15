
## TO DO HERE THE VIRAL LOAD AT INFECTION ? 

VL_infection <- function(xdata){
  vl_at_infection = NULL
  
  for(hh.index in 1:xdata$n_households){
    vl_at_infection[[hh.index]] = numeric(0)
    infection.events = xdata$infection.events[[hh.index]]
    X = xdata$all.HH %>% filter(household.id == hh.index) 
    if(length(infection.events>0)){
      vl = c()
      t.infection.0 = min(X$t.infection)
      
      X =X %>% mutate(times = t.infection-t.infection.0) %>% mutate(x.id = 1)
      w =which(X$times==0 )[1]
      X$times[-w] = X$times[-w] +0.000001 # solves an issue when the first infection happens at time 0
      times= unique(X$times[-1])
      
      for(t in times){
        w = which(abs(t-infection.events$time) <0.0001)
        infected = infection.events$infected[w] 
        X$x.id[which(X$times == t)] = infected
      }
      for(j in 1:dim(infection.events)[1]){
        infector = infection.events$infector[j]
        # print(infector)
        Xinfector = X%>% filter(x.id == infector)
        # 
        if(nrow(Xinfector >1)){
          Xinfector=Xinfector[1,]
        }
        G = get_viral_load(params = Xinfector ,delay = Xinfector$times, Tmax = infection.events$time[j],dt = 0.01)
        if(sum(!is.na(G))){
          Gfinal= G[nrow(G) ,]
          # vl= c(vl, Gfinal$VI+Gfinal$VNI) 
          vl= c(vl, Gfinal$V) 
        }
      }
      vl_at_infection[[hh.index]] = vl
      
    }
  }
  return(vl_at_infection)
}

get_duration_epidemics <- function(xdata){
  durations = c()
  for(hh.index in 1:xdata$n_households){
    X = xdata$all.HH %>% filter(household.id == hh.index) %>% mutate(end.infection = t.infection+T_P+T_U) %>% summarize(d = max(end.infection)-min(t.infection))
    durations = c(durations, X$d)
  }
  return(durations)
}

count_generations <- function(mat){
  n=nrow(mat)
  s=rep(0, n)
  for( i in 1:n){
    infector = which(mat[i,] == 1)
    if(infector == i){
      index.case =i
      s[index.case]=0
    }
  }
  
  for( i in 1:nrow(mat)){
    infector = i
    while(infector != index.case){
      infector = which(mat[infector,] == 1)
      # if(infector == i){
      s[i] = s[i]+1
      # }
    } 
  }
  
  return(s)
}

estimate_vl_at_infection <- function(mat,TP,TU,VP,tinfection){
  
  # For each infectee infected couple, measure the difference in time of infection
  mat = as.matrix(mat)
  vl=c()
  params = data.frame(T_P = TP,
                      T_U = TU,
                      V_P = VP,
                      tinfection = tinfection-min(tinfection))
  
  for(infector in 1:ncol(mat)){
    for(infected in 1:nrow(mat)){
      if(infector!=infected & mat[infected,infector] ==1){
        G = model_gen_triangle(T_P = params[infector,]$T_P, T_U= params[infector,]$T_U, V_P=params[infector,]$V_P, t = params$tinfection[infected]-params$tinfection[infector],tinfection=0, normalize.time = FALSE)
        # vl= c(vl, G$VI )
        vl= c(vl, G$V )
      }
    }
  }
  
  return(vl)
  
  
}

get_generation_interval_one_household <- function(mat,tinfection){
  # For each infector/infected couple, measure the difference in time of infection
  mat = as.matrix(mat)
  n=nrow(mat)
  s=rep(0, n)
  
  generation_interval = c()
  
  for( i in 1:n){
    infector = which(mat[i,] == 1)
    if(infector != i){
      generation_interval = c(generation_interval, tinfection[i]-tinfection[infector])
    }
  }
  return(generation_interval)
  
}


get_generation_interval<- function(xdata){
  
  n_households = xdata$n_households
  generation.times = NULL
  chains.tr = get_transmission_chain(xdata)
  for(hh.index in 1:n_households){
    t.infections = xdata$all.HH %>% filter(household.id == hh.index) %>% select(t.infection)
    mat=chains.tr[[hh.index]]
    generation.times[[hh.index]] = numeric(0)
    est = get_generation_interval_one_household(chains.tr[[hh.index]],tinfection = as.matrix(t.infections) )
    if(!is.null(est)){
      generation.times[[hh.index]] =est
    }
  }
  return(generation.times)
}

get_serial_interval<- function(xdata){
  n_households = xdata$n_households
  serial.intervals = NULL
  chains.tr = get_transmission_chain(xdata)
  for(hh.index in 1:n_households){
    serial.intervals[[hh.index]] = numeric(0)
    t.symptoms = xdata$all.HH %>% filter(household.id == hh.index) %>% select(t.symptom)
    is.symptomatic = xdata$symptomatic_array[hh.index,] 
    is.symptomatic[which(is.symptomatic == 0)] = NA
    if( xdata$n_infected[hh.index]>1){ # if there is an infection
      se = get_generation_interval_one_household(chains.tr[[hh.index]],tinfection = as.matrix(t.symptoms*is.symptomatic) )
      serial.intervals[[hh.index]] = se[which(!is.na(se))]
    }
  }
  return(serial.intervals)
}

get_transmission_chain<- function(xdata){
  
  n_households = xdata$n_households
  input.chain = NULL
  for(hh.index in 1:n_households){
    size_hh = xdata$size_households[hh.index]
    
    mat = matrix(data=0, nrow = size_hh, ncol = size_hh)
    mat[1,1] =1
    infection.events = xdata$infection.events[[hh.index]]
    
    for(infected in 1:size_hh){
      for( infector in 1:size_hh){
        if(sum(infection.events$infected == infected & infection.events$infector == infector) )
          mat[infected,infector] = 1
      }
    }
    # delete matrix column and rows that correspond to non infected 
    x= which(rowSums(mat) >0) # who are the infected persons? 
    mat = mat[x,x]
    input.chain[[hh.index]]= mat
  } 
  return(input.chain)
}


get_R0_chains <-function(mat){
  mat=as.matrix(mat)
  n.infected = colSums(mat) # a vector with the number of transmissions per individuals (even the non infected ones)
  n.infected[1] = n.infected[1] -1 # because we put the matrix entry of the index case as 1 for convenience elsewhere
  
  is.infected = rowSums(mat) >0
  n.infected2 = c()
  for(i in 1:nrow(mat)){
    if(is.infected[i]){
      n.infected2=c(n.infected2, n.infected[i])
    }
  }
  R0 = mean(n.infected2)
  return(R0)
  
}

## encore un bug ici
get_number_infections_index <-function(mat){
  mat=as.matrix(mat)
  
  index = which(diag(mat) ==1)
  n.infected = sum(mat[,index])-1 # a vector with the number of transmissions per individuals (even the non infected ones)
  return(n.infected)
}

reconstruct_transmission_chain <- function(xdata, Chains, hh.index, model = 'exp', chain_posterior = TRUE){ #, size.parameter.alpha= 0){#size.dependency = "density" ){
  # reconstruct_transmission_chain<- function(xdata, Chains, hh.index, model = 'exp', chain_posterior = TRUE){#size.dependency = "density" ){
  
  n.infected = xdata$n_infected[hh.index]
  size_hh = xdata$size_households[hh.index]
  input.chain = get_transmission_chain(xdata)[[hh.index]]
  # get the interval between infections in the input
  input.generation.times = get_generation_interval(xdata )[[hh.index]]
  input.serial.interval = get_serial_interval(xdata )[[hh.index]]
  input.R0 = get_R0_chains(input.chain)
  input.infections.index = get_number_infections_index(input.chain)
  input.duration.epidemics = get_duration_epidemics(xdata)[hh.index]
  input.vl.at.infection = VL_infection(xdata)[hh.index]
  #t.infection.input = xdata$all.HH %>% filter(household.id == hh.index) %>% select(t.infection)
  if(chain_posterior==FALSE){
    M=1
  }else{
    M = length(Chains$m)
  }
  P = difference.edges = gen = generation.times =serial.times = R0 = vl.at.infection = NULL
  difference.edges = R0 = gen = duration.epidemics = number.infections.index = rep(0,M)
  
  is.symptomatic = xdata$symptomatic_array[hh.index,] 
  is.symptomatic[which(is.symptomatic == 0)] = NA
  if(chain_posterior == TRUE){
    duration.epidemics = Chains$T_P[,hh.index,1]+ Chains$T_U[,hh.index,1]
    
    if(n.infected>1){
      h = xdata$all.HH$h[1]
      prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
      P = matrix(data=0, nrow = n.infected, ncol = n.infected)
      generation.times=c()
      serial.times = c()
      vl.at.infection = c()
      # prob.infection[1,1] = 1
      
      for( s in 1:M){
        m=Chains$m[s]
        TP = Chains$T_P[s,hh.index,]
        TU = Chains$T_U[s,hh.index,]
        VP = Chains$V_P[s,hh.index,]
        tinfection = Chains$tinfection[s,hh.index,]
        TS = Chains$TS[s,hh.index,]
        size.parameter.alpha.fitted= Chains$alpha[s]
        
        
        transmission.factor.size.fitted = 1/(size_hh-1)^size.parameter.alpha.fitted
        titers.index.case= NULL
        
        for(j in 1:n.infected){
          titers.index.case[[j]] = get_viral_load(params = data.frame(T_P = TP[j], T_U = TU[j], V_P = VP[j]), delay= 0, Tmax = 30, dt = 0.01) %>% 
            #   mutate(prob.transmission = m*(VI)^h ) %>%
            mutate(t = t + tinfection[j])
          
          if(model == "exp"){
            titers.index.case[[j]] = titers.index.case[[j]] %>% mutate(prob.transmission = transmission.factor.size.fitted*m*(V)^h ) 
          }
          if(model == "log"){
            titers.index.case[[j]] = titers.index.case[[j]] %>% mutate(prob.transmission = ifelse(log10(V)>= infection.threshold, transmission.factor.size.fitted*m*log10(V), 0) )
          } 
          if(model == "constant"){
            # titers.index.case[[j]] = titers.index.case[[j]]%>% mutate(prob.transmission = ifelse(log10(VI)>=2, m, 0) )
          }
        }
        prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
        for(infected in 1:n.infected){
          index.possible.infector=c()
          pt=c()
          i =0
          for(infector in 1:n.infected){
            TI= titers.index.case[[infector]]
            if(infected!=infector & tinfection[infector]<= tinfection[infected]){
              i=i+1
              x = which.min(abs(TI$t - tinfection[infected] ))
              pt[i] = TI$prob.transmission[x]
              index.possible.infector[i] = infector
            }  
          }
          pt = pt/sum(pt)
          if(length(pt)==0){# index case
            prob.infection[infected,infected] = prob.infection[infected,infected] + 1
          }else{
            r=runif(n=1)
            x = which(r<cumsum(pt))[1]
            infector = index.possible.infector[x]
            prob.infection[infected,infector] = prob.infection[infected,infector] + 1
          }
          ## randomly choose the infector of infected 
        }
        # reconstructed.chain[[s]] = prob.infection
        
        P=P + prob.infection 
        difference.edges[ s]= sum(abs(prob.infection-input.chain))
        gen[s] = max(count_generations(prob.infection))
        generation.times = c(generation.times,get_generation_interval_one_household(prob.infection, tinfection))
        number.infections.index[s]=get_number_infections_index(prob.infection)
        R0[s] = get_R0_chains(prob.infection)
        duration.epidemics[s]= max(TP+TU+tinfection-min(tinfection))
        tsymptoms = TS+tinfection
        se = get_generation_interval_one_household(prob.infection,tinfection = as.matrix(tsymptoms*is.symptomatic) )
        if(!is.null(se[which(!is.na(se))])){
          serial.times = c(serial.times,se[which(!is.na(se))])
        }
        vl.at.infection = c(vl.at.infection,estimate_vl_at_infection(prob.infection, TP,TU,VP,tinfection )) 
      }
      P = P/M
    }
  }
  
  return(list(input.chain = input.chain,
              input.generation.times = input.generation.times,
              input.R0 = input.R0,
              input.serial.interval = input.serial.interval,
              P = P,
              difference.edges=difference.edges,
              input.number.generations = max(count_generations(as.matrix(input.chain))),
              number.generations = gen,
              generation.times=generation.times,
              R0 = R0,
              serial.times = serial.times,
              number.infections.index = number.infections.index,
              input.infections.index = input.infections.index,
              input.duration.epidemics = input.duration.epidemics,
              duration.epidemics = duration.epidemics,
              input.vl.at.infection = input.vl.at.infection,
              vl.at.infection=vl.at.infection))
  
}

reconstruct_transmission_chain_optimized_input<- function(xdata, Chains =FALSE, model = 'exp', chain_posterior = FALSE){
  
  n.hh = length(xdata$n_infected)
  n.infected = xdata$n_infected#[hh.index]
  size_hh = xdata$size_households#[hh.index]
  input.chain = get_transmission_chain(xdata)#[[hh.index]]
  # get the interval between infections in the input
  input.generation.times = get_generation_interval(xdata )#[[hh.index]]
  input.serial.interval = get_serial_interval(xdata )#[[hh.index]]
  input.R0=input.infections.index=c()
  for(hh.index in 1:n.hh){
    input.R0[hh.index] = get_R0_chains(input.chain[[hh.index]])
    input.infections.index[hh.index] = get_number_infections_index(input.chain[[hh.index]])
  }
  
  input.duration.epidemics = get_duration_epidemics(xdata)#[hh.index]
  input.vl.at.infection = VL_infection(xdata)#[hh.index]
  
  M=1
  
  P = difference.edges = gen = generation.times =serial.times = R0 = vl.at.infection = NULL
  difference.edges = rep(0,M)
  R0 = rep(0,M)
  gen = rep(0,M)
  duration.epidemics = rep(0,M)
  number.infections.index = rep(0,M)
  
  GEN = c()
  for(hh.index in 1:n.hh){
    GEN[hh.index] = max(count_generations(as.matrix(input.chain[[hh.index]])))
  }
  
  return(list(input.chain = input.chain,
              input.generation.times = input.generation.times,
              input.R0 = input.R0,
              input.serial.interval = input.serial.interval,
              P = P,
              difference.edges=difference.edges,
              input.number.generations = GEN,# max(count_generations(as.matrix(input.chain))),
              number.generations = gen,
              generation.times=generation.times,
              R0 = R0,
              serial.times = serial.times,
              number.infections.index = number.infections.index,
              input.infections.index = input.infections.index,
              input.duration.epidemics = input.duration.epidemics,
              duration.epidemics = duration.epidemics,
              input.vl.at.infection = input.vl.at.infection,
              vl.at.infection=vl.at.infection))
  
}

get_CT<- function(all.CT, hh.index){
  return(list(input.chain = all.CT$input.chain[[hh.index]],
              input.generation.times = all.CT$input.generation.times[[hh.index]],
              input.R0 = all.CT$input.R0[[hh.index]],
              input.serial.interval = all.CT$input.serial.interval[[hh.index]],
              P = all.CT$P,
              difference.edges=all.CT$difference.edges,
              input.number.generations =all.CT$input.number.generations[[hh.index]],
              number.generations =all.CT$number.generations,
              generation.times=all.CT$generation.times,
              R0 =all.CT$R0,
              serial.times =all.CT$serial.times,
              number.infections.index =all.CT$number.infections.index,
              input.infections.index =all.CT$input.infections.index[hh.index],
              input.duration.epidemics =all.CT$input.duration.epidemics[hh.index],
              duration.epidemics =all.CT$duration.epidemics,
              input.vl.at.infection =all.CT$input.vl.at.infection[hh.index],
              vl.at.infection=all.CT$vl.at.infection))
  
}

reconstruct_transmission_chain_symptoms <- function(xdata, Chains, hh.index, chain_posterior = TRUE){ #, size.parameter.alpha= 0){#size.dependency = "density" ){
  prob.pre.symptoms = 0.0000001
  
  n.infected = xdata$n_infected[hh.index]
  size_hh = xdata$size_households[hh.index]
  input.chain = get_transmission_chain(xdata)[[hh.index]]
  # get the interval between infections in the input
  input.generation.times = get_generation_interval(xdata )[[hh.index]]
  input.serial.interval = get_serial_interval(xdata )[[hh.index]]
  input.R0 = get_R0_chains(input.chain)
  input.infections.index = get_number_infections_index(input.chain)
  input.duration.epidemics = get_duration_epidemics(xdata)[hh.index]
  input.vl.at.infection = VL_infection(xdata)[hh.index]
  #t.infection.input = xdata$all.HH %>% filter(household.id == hh.index) %>% select(t.infection)
  if(chain_posterior==FALSE){
    M=1
  }else{
    M = length(Chains$m)
  }
  P = difference.edges = gen = generation.times =serial.times = R0 = vl.at.infection = NULL
  difference.edges = rep(0,M)
  R0 = rep(0,M)
  gen = rep(0,M)
  duration.epidemics = rep(0,M)
  number.infections.index = rep(0,M)
  
  is.symptomatic = xdata$symptomatic_array[hh.index,] 
  is.symptomatic[which(is.symptomatic == 0)] = NA
  
  if(chain_posterior == TRUE){
    
    if(n.infected>1){
      
      h = xdata$all.HH$h[1]
      prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
      P = matrix(data=0, nrow = n.infected, ncol = n.infected)
      generation.times=c()
      serial.times = c()
      vl.at.infection = c()
      for( s in 1:M){
        m=Chains$m[s]
        gamma_infection_shape = Chains$gamma_infection_shape[s]
        gamma_infection_inverse_scale = Chains$gamma_infection_inverse_scale[s]
        tinfection = Chains$tinfection[s,hh.index,]
        TS = Chains$TS[s,hh.index,]
        size.parameter.alpha.fitted= Chains$alpha[s]
        transmission.factor.size.fitted = 1/(size_hh-1)^size.parameter.alpha.fitted
        titers.index.case= NULL
        tlag = tinfection 
        prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
        
        for(infected in 1:n.infected){
          index.possible.infector=c()
          pt=c()
          i =0
          for(infector in 1:n.infected){
            
            if(infected!=infector & tinfection[infector]<= tinfection[infected]){
              i=i+1
              x = tinfection[infected]-tinfection[infector]
              pt[i] = dgamma(x,shape =gamma_infection_shape,scale =1/gamma_infection_inverse_scale )
              index.possible.infector[i] = infector
            }  
          }
          pt = pt/sum(pt)
          if(length(pt)==0){# index case
            prob.infection[infected,infected] = prob.infection[infected,infected] + 1
          }else{
            r=runif(n=1)
            x = which(r<cumsum(pt))[1]
            infector = index.possible.infector[x]
            prob.infection[infected,infector] = prob.infection[infected,infector] + 1
          }
          ## randomly choose the infector of infected 
        }
        # reconstructed.chain[[s]] = prob.infection
        
        P=P + prob.infection 
        difference.edges[ s]= sum(abs(prob.infection-input.chain))
        gen[s] = max(count_generations(prob.infection))
        generation.times = c(generation.times,get_generation_interval_one_household(prob.infection, tinfection))
        number.infections.index[s]=get_number_infections_index(prob.infection)
        R0[s] = get_R0_chains(prob.infection)
        #duration.epidemics[s]= max(TP+TU+tinfection-min(tinfection))
        
        tsymptoms = TS+tinfection
        se = get_generation_interval_one_household(prob.infection,tinfection = as.matrix(tsymptoms*is.symptomatic) )
        if(!is.null(se[which(!is.na(se))])){
          serial.times = c(serial.times,se[which(!is.na(se))])
        }
        
        #vl.at.infection = c(vl.at.infection,estimate_vl_at_infection(prob.infection, TP,TU,VP,tinfection )) 
        
      }
      P = P/M
    }
  }
  
  return(list(input.chain = input.chain,
              input.generation.times = input.generation.times,
              input.R0 = input.R0,
              input.serial.interval = input.serial.interval,
              P = P,
              difference.edges=difference.edges,
              input.number.generations = max(count_generations(as.matrix(input.chain))),
              number.generations = gen,
              generation.times=generation.times,
              R0 = R0,
              serial.times = serial.times,
              number.infections.index = number.infections.index,
              input.infections.index = input.infections.index,
              input.duration.epidemics = input.duration.epidemics,
              duration.epidemics = duration.epidemics,
              input.vl.at.infection = input.vl.at.infection,
              vl.at.infection=vl.at.infection))
  
}



reconstruct_transmission_chain_positivity <- function(xdata, Chains, hh.index, chain_posterior = TRUE){ #, size.parameter.alpha= 0){#size.dependency = "density" ){
  prob.pre.symptoms = 0.0000001
  
  n.infected = xdata$n_infected[hh.index]
  size_hh = xdata$size_households[hh.index]
  input.chain = get_transmission_chain(xdata)[[hh.index]]
  # get the interval between infections in the input
  input.generation.times = get_generation_interval(xdata )[[hh.index]]
  input.serial.interval = get_serial_interval(xdata )[[hh.index]]
  input.R0 = get_R0_chains(input.chain)
  input.infections.index = get_number_infections_index(input.chain)
  input.duration.epidemics = get_duration_epidemics(xdata)[hh.index]
  input.vl.at.infection = VL_infection(xdata)[hh.index]
  #t.infection.input = xdata$all.HH %>% filter(household.id == hh.index) %>% select(t.infection)
  if(chain_posterior==FALSE){
    M=1
  }else{
    M = length(Chains$m)
  }
  P = difference.edges = gen = generation.times =serial.times = R0 = vl.at.infection = NULL
  difference.edges = rep(0,M)
  R0 = rep(0,M)
  gen = rep(0,M)
  duration.epidemics = rep(0,M)
  number.infections.index = rep(0,M)
  
  is.symptomatic = xdata$symptomatic_array[hh.index,] 
  is.symptomatic[which(is.symptomatic == 0)] = NA
  
  if(chain_posterior == TRUE){
    
    if(n.infected>1){
      
      h = xdata$all.HH$h[1]
      prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
      P = matrix(data=0, nrow = n.infected, ncol = n.infected)
      generation.times=c()
      serial.times = c()
      vl.at.infection = c()
      for( s in 1:M){
        m=Chains$m[s]
        gamma_infection_shape = Chains$gamma_shape[s,,]
        gamma_infection_inverse_scale = Chains$gamma_infection_inverse_scale[s]
        tinfection = Chains$tinfection[s,hh.index,]
        TS = Chains$TS[s,hh.index,]
        size.parameter.alpha.fitted= Chains$alpha[s]
        transmission.factor.size.fitted = 1/(size_hh-1)^size.parameter.alpha.fitted
        titers.index.case= NULL
        tlag = tinfection 
        prob.infection = matrix(data=0, nrow = n.infected, ncol = n.infected)
        
        for(infected in 1:n.infected){
          index.possible.infector=c()
          pt=c()
          i =0
          for(infector in 1:n.infected){
            
            if(infected!=infector & tinfection[infector]<= tinfection[infected]){
              i=i+1
              x = tinfection[infected]-tinfection[infector]
              pt[i] = dgamma(x,shape =gamma_infection_shape[infector],scale =1/gamma_infection_inverse_scale )
              index.possible.infector[i] = infector
            }  
          }
          pt = pt/sum(pt)
          if(length(pt)==0){# index case
            prob.infection[infected,infected] = prob.infection[infected,infected] + 1
          }else{
            r=runif(n=1)
            x = which(r<cumsum(pt))[1]
            infector = index.possible.infector[x]
            prob.infection[infected,infector] = prob.infection[infected,infector] + 1
          }
          ## randomly choose the infector of infected 
        }
        # reconstructed.chain[[s]] = prob.infection
        
        P=P + prob.infection 
        difference.edges[ s]= sum(abs(prob.infection-input.chain))
        gen[s] = max(count_generations(prob.infection))
        generation.times = c(generation.times,get_generation_interval_one_household(prob.infection, tinfection))
        number.infections.index[s]=get_number_infections_index(prob.infection)
        R0[s] = get_R0_chains(prob.infection)
        #duration.epidemics[s]= max(TP+TU+tinfection-min(tinfection))
        
        tsymptoms = TS+tinfection
        se = get_generation_interval_one_household(prob.infection,tinfection = as.matrix(tsymptoms*is.symptomatic) )
        if(!is.null(se[which(!is.na(se))])){
          serial.times = c(serial.times,se[which(!is.na(se))])
        }
        
        #vl.at.infection = c(vl.at.infection,estimate_vl_at_infection(prob.infection, TP,TU,VP,tinfection )) 
        
      }
      P = P/M
    }
  }
  
  return(list(input.chain = input.chain,
              input.generation.times = input.generation.times,
              input.R0 = input.R0,
              input.serial.interval = input.serial.interval,
              P = P,
              difference.edges=difference.edges,
              input.number.generations = max(count_generations(as.matrix(input.chain))),
              number.generations = gen,
              generation.times=generation.times,
              R0 = R0,
              serial.times = serial.times,
              number.infections.index = number.infections.index,
              input.infections.index = input.infections.index,
              input.duration.epidemics = input.duration.epidemics,
              duration.epidemics = duration.epidemics,
              input.vl.at.infection = input.vl.at.infection,
              vl.at.infection=vl.at.infection))
  
}


