generate_parameters_triangle <- function(T_P_star, T_U_star, V_P_star, Delta_S_star,
                                         T_P_eta, T_U_eta, V_P_eta, Delta_S_eta,
                                         m, 
                                         h = 0.5,
                                         total.samples = NULL,
                                         epid.model){
  #  source("main.R")
  model.virus ="triangle"
  #detach("package:plyr", unload = TRUE)
  fixed.parameters= TRUE
  
  DF = NULL
  Tmax = 70
  n.samples =1
  T_P = T_U= V_P = incubation = c()
  if(is.null(total.samples)){
    total.samples = length(m)
    
  }
  for(i in 1:total.samples){
    T_P[i] = T_P_star*exp(rnorm(n = n.samples, sd = T_P_eta)) 
    T_U[i] = T_U_star*exp(rnorm(n = n.samples, sd = T_U_eta))  
    V_P[i] = V_P_star*exp(rnorm(n = n.samples, sd = V_P_eta))
    incubation[i] = Delta_S_star*exp(rnorm(n = n.samples, sd = Delta_S_eta))
  }
  
  params = list(T_P = T_P,
                T_U = T_U,
                V_P = V_P,
                incubation = incubation,
                m  = m,
                h = h)
  #m = -log(1-SAR.target.sd) /( 0.1/(h*V_P_star*log(10))*(T_P_star*(10^(h*(V_P_star+4))-1) +T_U_star*(10^(h*V_P_star)-1)  )) 
  
  if(epid.model == "exp"){
   #SAR = 1-exp(-m*10^(2*h)/(h*V_P*log(10))*(T_P+T_U)*(10^(h*V_P)-1))
    #  SAR = 1-exp(-m*0.1/(h*V_P_star*log(10))*(T_P_star*(10^(h*(V_P_star+4))-1) +T_U_star*(10^(h*V_P_star)-1)  ))
    
   # Integrale =10^(-2*h)/(h*V_P_star*log(10))*T_P_star*(10^(h*(V_P_star+4))-1) +  10^(2*h)*(T_U_star/V_P_star/(h*log(10)))*(1-10^(-h*V_P_star))
    Integrale =10^(2*h)/(h*V_P*log(10))*(T_P+T_U)*(10^(h*(V_P))-1) #+  10^(2*h)*(T_U_star/V_P_star/(h*log(10)))*(1-10^(-h*V_P_star))
    SAR=1-exp(-m*Integrale)
   
  }
  if(epid.model == "constant"){
    SAR = 1-exp(-m*(T_P+T_U))
  }
  
 # print(SAR)

  # SAR = df$SAR[1]
  d = data.frame(id = seq(1,total.samples),
                 T_P = T_P,
                 T_U = T_U,
                 V_P = V_P,
                 incubation = incubation,
                 m = m, 
                 h=h,
                 SAR =SAR, 
                 epid.model )
  DF= rbind(DF,d)
  
  return(DF)
}



SAR_vs_m <-function(T_P_star = 5,T_U_star = 5,V_P_star = 2,Delta_S_star = 5, M,h=0.5, epid.model){
  SAR = c()
  for(m in M){
    G = generate_parameters_triangle(T_P_star = T_P_star,T_U_star = T_U_star,V_P_star = V_P_star,Delta_S_star =Delta_S_star,
                                     T_P_eta = 0,T_U_eta = 0,V_P_eta = 0,Delta_S_eta = 0, m=m,h=h, total.samples = 1,epid.model)
    
    SAR = c(SAR,G$SAR)  
  }
  return(SAR)
}

#param.list generated with get_individual_statistics
get_individual_statistics <-function(param.list){
  # durée d'infectivité
  return(list(SAR = param.list$SAR,
              infectivity = param.list$T_P+param.list$T_U,
              AUC = 1/2*(param.list$T_P+param.list$T_U)*param.list$V_P ,
              T_P = param.list$T_P,
              T_U  = param.list$T_U,
              V_P = param.list$V_P ))
  
}

## generate transmission pairs
get_transmission_statistics<-function(T_P_star = 5,T_U_star = 5,V_P_star = 2,Delta_S_star = 5, 
                                      T_P_eta = 0.4,T_U_eta = 0.4,V_P_eta = 0.2,Delta_S_eta = 0.125, m=1.2* 10^-3 , h=0.5,
                                      total.samples=2000,
                                      epid.model){
  model.virus = "triangle"
  n=0
  generation.interval = serial.interval = c()
  
  T_P_infector=T_U_infector=V_P_infector=Delta_S_infector=c()
  presymptomatic.infection = c()
  while(n <total.samples){
    
    param.list = generate_parameters_triangle(T_P_star = T_P_star,
                                              T_U_star = T_U_star,
                                              V_P_star = V_P_star,
                                              Delta_S_star = Delta_S_star, 
                                              T_P_eta = T_P_eta,
                                              T_U_eta = T_U_eta,
                                              V_P_eta = V_P_eta,
                                              Delta_S_eta = Delta_S_eta,
                                              m = m,
                                              h = h,
                                              total.samples = 2,
                                              epid.model)
    
    HH.parameters  = make_household(param.list,household.size =  2)
    
    list.transmission <- generate_transmission_short(HH.parameters)
    is.infected = list.transmission$n.infection==2
    
    if(is.infected){
      # get statistics about the infected/ infected pair
      n=n+1
      # print(n)
      # serial interval : 
      t.symptoms.1 =    param.list$incubation[1]
      t.symptoms.2 =   list.transmission$t.infection + param.list$incubation[2]
      serial.interval[n] = t.symptoms.2-t.symptoms.1
      generation.interval[n] =  list.transmission$t.infection 
      T_P_infector[n] =  param.list$T_P[1]
      T_U_infector[n] =  param.list$T_U[1]
      V_P_infector[n] =  param.list$V_P[1]
      Delta_S_infector[n] =  param.list$incubation[1]
      presymptomatic.infection[n] = list.transmission$t.infection  < t.symptoms.1
      
    }
  }
  return(list(T_P_star = T_P_star,
              T_U_star = T_U_star,
              V_P_star = V_P_star,
              Delta_S_star = Delta_S_star, 
              T_P_eta = T_P_eta,
              T_U_eta = T_U_eta,
              V_P_eta = V_P_eta,
              Delta_S_eta = Delta_S_eta,
              m = m,
              serial.interval = serial.interval,
              generation.interval = generation.interval,
              T_P_infector=T_P_infector,
              T_U_infector=T_U_infector,
              V_P_infector=V_P_infector,
              Delta_S_infector=Delta_S_infector,
              prop.negative.serial.interval = length(which(serial.interval<0))/length(serial.interval),
              presymptomatic.infection = mean(presymptomatic.infection))
  )
  
}


## generate transmission pairs
infectiousness<-function(T_P_star,
                         T_U_star,
                         V_P_star,
                         Delta_S_star, 
                         T_P_eta,
                         T_U_eta,
                         V_P_eta,
                         Delta_S_eta,
                         m=m,
                         h,
                         total.samples=2000,
                         Tmax=20,
                         epid.model){
  model.virus = "triangle"
  
  generation.interval = serial.interval = c()
  
  T_P_infector=T_U_infector=V_P_infector=Delta_S_infector=c()
  presymptomatic.infection = c()
  dt = 0.05
  N=Tmax/dt+1
  prob = rep(0,N)
  for(i in 1:total.samples){
    
    param.list = generate_parameters_triangle(T_P_star = T_P_star,
                                              T_U_star = T_U_star,
                                              V_P_star = V_P_star,
                                              Delta_S_star = Delta_S_star, 
                                              T_P_eta = T_P_eta,
                                              T_U_eta = T_U_eta,
                                              V_P_eta = V_P_eta,
                                              Delta_S_eta = Delta_S_eta,
                                              m = m,
                                              h = h,
                                              total.samples = 1,
                                              epid.model)
    
    index.case = HH %>% filter(status == 'index')
    titers.index.case = get_viral_load(params = param.list , delay= 0, Tmax =  Tmax ,dt = dt)    %>% 
      mutate(prob.transmission =  dt*m*(VI)^h ) %>%
      mutate(cumulative.prob.transmission = cumsum(prob.transmission))
    
    prob = prob + titers.index.case$prob.transmission/total.samples
    
  }
  
  return(data.frame(t = titers.index.case$t, prob.infection  = prob))
  
}







