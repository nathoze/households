# 
# source("MCMC_generate_observed_VI_symptoms.R")

# number.hh = 5
# size.hh = c(5)
# obs.time.step = 3
# m = 1.2* 10^-3
# noise = 0.3
# Household.list = NULL
# 
# 
# A = generate_observations(number.hh = number.hh  ,size.hh = size.hh, obs.time.step = obs.time.step, noise=0, m = m40)
# 


colors = c('#c45441',"#1212ca","#128a18","#11a1a8","#000000")
for(k in 1:A$n_households){
  
  params = A$all.HH
  Xlim = 70
  discrete_times = seq(0,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  
  g=ggplot()
  
  P = params %>% 
    filter(household.id == k) %>%
    mutate(tinfection = t.infection)
  for(j in 1:A$n_infected[k]){
    
    P1=P[j,]
    
    X = model_gen_triangle(T_P = P1$T_P,
                           T_U = P1$T_U, 
                           V_P = P1$V_P,  
                           t = seq(0, max(discrete_times), by=0.1), 
                           tinfection = 0 )
    X = X%>% 
      mutate(t = t+P1$t.infection) %>%
      mutate(logV = log10(VNI))  %>% 
      mutate(logV = ifelse(logV>6,logV,6)) %>%
      #mutate(t=t+ P1$t.infection)  %>% 
      filter(t <=Tmax)
    
    # obs.data = data.frame(time=A$time_common,   
    #                       x = A$LVLObs[L$start[j]:L$end[j]])
    
    obs.data = data.frame(time=A$time[k,A$start[k,j]:A$end[k,j]],   
                          x =  A$LVLObs[k,A$start[k,j]:A$end[k,j]])
    real.curve = data.frame(time= X$t, y = X$logV)
    Xlim=25
    g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j],size=3) + 
      geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
      xlim(c(-10,Xlim))+
      #  ylim(c(0,13))+
      theme_bw() +
      xlab('Days post inclusion')+
      ylab('Log Viral load')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+
      geom_vline(xintercept = 0, color = 'black')+
      ggtitle(paste0("Household ", k))
 
    
  }
  
  for(x in A$time[k,A$start[k,j]:A$end[k,j]]){
    g=g +geom_vline(xintercept = x, color = 'grey', linetype= 'dashed')
  }
  print(g)
}


