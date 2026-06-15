
MCMC_plot_predicted_viral_load_hhindex<- function(Chains, xdata, Xlim = Tmax,hh.index){
  
  G=NULL
  
  
  # load(filename)
  #  Chains= rstan::extract(fit)
  S = length(Chains$logm)
  discrete_times = seq(-10,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  
  colors = c('#c45441',"#1212ca","#128a18","#ffb20d")
  g=ggplot()
  for(j in 1:xdata$n_infected[hh.index]){
    print(j)
    
    LVL = array(data=NA, dim = c(S, n_discrete_times))
    for (i in 1:S){
      mod =   model_gen_triangle(T_P = Chains$T_P[i,hh.index,j], 
                                 T_U = Chains$T_U[i,hh.index,j],
                                 V_P = Chains$V_P[i,hh.index,j],
                                 t =  discrete_times,
                                 tinfection = Chains$tinfection[i,hh.index,j],
                                 normalize.time = FALSE)
      LVL[i,] = mod$logV
      
    }
    
    med = apply(LVL,2,median)
    lower = apply(LVL,2,quantile025)
    upper = apply(LVL,2,quantile975)
    
    
    P1= xdata$all.HH %>% 
      filter(household.id == hh.index) %>% 
      mutate(tinfection = t.infection)
    
    P1=P1[j,]
    
    X = model_gen_triangle(T_P = P1$T_P,
                           T_U = P1$T_U, 
                           V_P = P1$V_P,  
                           t = seq(P1$t.infection, max(discrete_times), by=0.1), 
                           tinfection = P1$t.infection,normalize.time = FALSE)
    
    X = X %>%
      mutate(logV = log10(V))  %>% 
      mutate(logV = ifelse(logV>detection.limit,logV,detection.limit)) %>%
      filter(t <=Tmax)
    
    pred.data =  data.frame(time=discrete_times,
                            mean = med, 
                            quantile025=lower, 
                            quantile975=upper) %>% 
      mutate(mean = ifelse(mean>detection.limit,mean,detection.limit)) %>%
      mutate(quantile025 = ifelse(quantile025>detection.limit,quantile025,detection.limit)) %>%
      mutate(quantile975 = ifelse(quantile975>detection.limit,quantile975,detection.limit))
    
    obs.data = data.frame(time=xdata$time[hh.index,xdata$start[hh.index,j]:xdata$end[hh.index,j]],   
                          x = xdata$LVLObs[hh.index,xdata$start[hh.index,j]:xdata$end[hh.index,j]])
    real.curve = data.frame(time= X$t, y = X$logV)
    
    g =g+       geom_point(data = as.data.frame(obs.data),aes(time,x),shape = 21,color="black",fill= colors[j],size=4) + 
      geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors[j])+
      geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
      geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors[j]) +
       xlim(c(NA,Xlim))+
      #  ylim(c(0,13))+
      theme_bw() +
      xlab('Days')+
      ylab('Viral load (log10 copies/mL)')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))#+ggtitle(paste0('Household ', hh.index))
  }
  
  
  return(g)
}


MCMC_plot_predicted_viral_load<- function(Chains, params,xdata, Xlim = Tmax){
  
  G=NULL
  # load(filename)
  #  Chains= rstan::extract(fit)
  S = length(Chains$logm)
  discrete_times = seq(0,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  
  colors = c('#c45441',"#1212ca","#128a18","#ffb20d")
  
  
  
  for(k in 1:xdata$n_households ){
    g=ggplot()
    for(j in 1:xdata$n_infected[k]){
      
      LVL = array(data=NA, dim = c(S, n_discrete_times))
      for (i in 1:S){
        mod =   model_gen_triangle(T_P = Chains$T_P[i,k,j], 
                                   T_U = Chains$T_U[i,k,j],
                                   V_P = Chains$V_P[i,k,j],
                                   t =  discrete_times,
                                   tinfection = Chains$tinfection[i,k,j])
        LVL[i,] = mod$logV
        
      }
      
      med = apply(LVL,2,median)
      lower = apply(LVL,2,quantile025)
      upper = apply(LVL,2,quantile975)
      
      
      P1= params %>% 
        filter(household.id == k) %>% 
        mutate(tinfection = t.infection)
      
      P1=P1[j,]
      
      X = model_gen_triangle(T_P = P1$T_P,
                             T_U = P1$T_U, 
                             V_P = P1$V_P,  
                             t = seq(0, max(discrete_times), by=0.1), 
                             tinfection = P1$t.infection )
      
      X = X %>%
        mutate(logV = log10(V))  %>% 
        mutate(logV = ifelse(logV>detection.limit,logV,detection.limit)) %>%
        #mutate(t=t+ P1$t.infection)  %>% 
        filter(t <=Tmax)
      
      pred.data =  data.frame(time=discrete_times,
                              mean = med, 
                              quantile025=lower, 
                              quantile975=upper) %>% 
        mutate(mean = ifelse(mean>detection.limit,mean,detection.limit)) %>%
        mutate(quantile025 = ifelse(quantile025>detection.limit,quantile025,detection.limit)) %>%
        mutate(quantile975 = ifelse(quantile975>detection.limit,quantile975,detection.limit))
      
      # obs.data = data.frame(time=xdata$time_common,   
      #                           x = xdata$LVLObs[k,xdata$start[k,j]:xdata$end[k,j]])
      obs.data = data.frame(time=xdata$time[k,xdata$start[k,j]:xdata$end[k,j]],   
                            x = xdata$LVLObs[k,xdata$start[k,j]:xdata$end[k,j]])
      real.curve = data.frame(time= X$t, y = X$logV)
      
      g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j]) + 
        geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors[j])+
        geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
        geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors[j]) +
        xlim(c(0,Xlim))+
        #  ylim(c(0,13))+
        theme_bw() +
        xlab('Days')+
        ylab('Log Viral load')+
        theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
              axis.text.y = element_text(size=16),
              text=element_text(size=16))+ggtitle(paste0('Household ', k))
    }
    #   print(g)
    G[[k]] = g
  }
  
  return(G)
}

MCMC_plot_predictive_posterior_check<- function(Chains, params,xdata, Xlim = Tmax){
  
  # load(filename)
  #  Chains= rstan::extract(fit)
  S = length(Chains$logm)
  G=NULL
  
  colors = c('#c45441',"#1212ca","#128a18","#ffb20d")
  Tmax= max(xdata$time_common)
  
  discrete_times = seq(0,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  
  for(k in 1:xdata$n_households ){
    g=ggplot()
    print(k)
    
    for(j in 1:xdata$n_infected[k]){
      
      LVL = array(data=NA, dim = c(S, n_discrete_times))
      for (i in 1:S){
        mod =   model_gen_triangle(T_P = Chains$T_P[i,k,j], 
                                   T_U = Chains$T_U[i,k,j],
                                   V_P = Chains$V_P[i,k,j],
                                   t = discrete_times,
                                   tinfection = Chains$tinfection[i,k,j])
        sigma_noise = Chains$sigma_noise[i]
        
        #    mod$
        #LVL[i,] = mod$logVNI  + rnorm(n=length( mod$logVNI ),mean = 0, sd= sigma_noise)
        LVL[i,] = mod$logV  + rnorm(n=length( mod$logV),mean = 0, sd= sigma_noise)
        
      }
      
      med = apply(LVL,2,median)
      lower = apply(LVL,2,quantile025)
      upper = apply(LVL,2,quantile975)
      
      P1= params %>% 
        filter(household.id == k) %>% 
        mutate(tinfection = t.infection)
      
      P1=P1[j,]
      
      X = model_gen_triangle(T_P = P1$T_P,
                             T_U = P1$T_U, 
                             V_P = P1$V_P,  
                             t = seq(0, Tmax, by=0.1), 
                             tinfection = P1$t.infection )
      
      X = X %>%
        mutate(logV = log10(V))  %>% 
        mutate(logV = ifelse(logV>detection.limit,logV,detection.limit)) %>%
        filter(t <=Tmax)
      
      pred.data =  data.frame(time=discrete_times,
                              mean = med, 
                              quantile025=lower, 
                              quantile975=upper) %>% 
        mutate(mean = ifelse(mean>detection.limit,mean,detection.limit)) %>%
        mutate(quantile025 = ifelse(quantile025>detection.limit,quantile025,detection.limit)) %>%
        mutate(quantile975 = ifelse(quantile975>detection.limit,quantile975,detection.limit))
      
      obs.data = data.frame(time=xdata$time_common,   
                            x = xdata$LVLObs[k,xdata$start[k,j]:xdata$end[k,j]])
      
      real.curve = data.frame(time= X$t, y = X$logV)
      
      g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j], size=3) + 
        geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors[j])+
        geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
        geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors[j]) +
        xlim(c(0,Xlim))+
        #   #  ylim(c(0,13))+
        theme_bw() +
        xlab('Days')+
        ylab('Log Viral load')+
        theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
              axis.text.y = element_text(size=16),
              text=element_text(size=16))+ggtitle(paste0('Household ', k))
      
    }
    G[[k]] = g
  }
  
  return(G)
} 

MCMC_plot_predictive_posterior_symptoms<- function(Chains, params,xdata, Xlim = Tmax, show.symptoms = FALSE, show.fit = TRUE){
  
  # load(filename)
  #  Chains= rstan::extract(fit)
  S = length(Chains$logm)
  G=NULL
  
  colors = c('#c45441',"#1212ca","#128a18","#11a1a8")
  colors = c('#c45441',"#1212ca","#128a18","#ffb20d","#776666")
  #  Tmax= max(xdata$time_common)
  Tmax = xdata$max_time
  discrete_times = seq(0,Tmax+15,by=0.1)
  n_discrete_times = length(discrete_times)
  
  offset.xlim.negative  = -15
  
  for(k in 1:xdata$n_households ){
    g=ggplot()
    print(k)
    
    for(j in 1:xdata$n_infected[k]){
      LVL = array(data=NA, dim = c(S, n_discrete_times))
      
      if(show.fit){ 
        for (i in 1:S){
          mod =   model_gen_triangle(T_P = Chains$T_P[i,k,j], 
                                     T_U = Chains$T_U[i,k,j],
                                     V_P = Chains$V_P[i,k,j],
                                     t = discrete_times+offset.xlim.negative,
                                     tinfection = Chains$tinfection[i,k,j],
                                     normalize.time=FALSE )
          sigma_noise = Chains$sigma_noise[i]
          
          #    mod$
          LVL[i,] = mod$logV  + rnorm(n=length( mod$logV ),mean = 0, sd= sigma_noise)
          
        }
        med = apply(LVL,2,median)
        lower = apply(LVL,2,quantile025)
        upper = apply(LVL,2,quantile975)
        
        pred.data =  data.frame(time=discrete_times+offset.xlim.negative,
                                mean = med, 
                                quantile025=lower, 
                                quantile975=upper) %>% 
          mutate(mean = ifelse(mean>detection.limit,mean,detection.limit)) %>%
          mutate(quantile025 = ifelse(quantile025>detection.limit,quantile025,detection.limit)) %>%
          mutate(quantile975 = ifelse(quantile975>detection.limit,quantile975,detection.limit))
        
      }
      
      P1= params %>% 
        filter(household.id == k) %>% 
        mutate(tinfection = t.infection)
      
      P1=P1[j,]
      # to 
      X = model_gen_triangle(T_P = P1$T_P,
                             T_U = P1$T_U, 
                             V_P = P1$V_P,  
                             t = seq(0, Tmax, by=0.1), 
                             tinfection =0, 
                             normalize.time=FALSE )
      X = X%>% 
        mutate(t = t+P1$t.infection) %>%
        mutate(logV = log10(V))  %>% 
        mutate(logV = ifelse(logV>detection.limit,logV,detection.limit)) %>%
        #mutate(t=t+ P1$t.infection)  %>% 
        filter(t <=Tmax)
      
      
      
      obs.data = data.frame(time=xdata$time[k,xdata$start[k,j]:xdata$end[k,j]],   
                            x = xdata$LVLObs[k,xdata$start[k,j]:xdata$end[k,j]])
      
      real.curve = data.frame(time= X$t, y = X$logV)
      
      # g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j], size=3) + 
      #   geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors[j])+
      #   geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
      #   geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors[j]) +
      #   xlim(c(offset.xlim.negative,Xlim))+
      #   #   #  ylim(c(0,13))+
      #   theme_bw() +
      #   xlab('Days')+
      #   ylab('Log Viral load')+
      #   theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
      #         axis.text.y = element_text(size=16),
      #         text=element_text(size=16))+
      #   geom_vline(xintercept = 0, color = 'black')+
      #   ggtitle(paste0('Household ', k))
      # 
      
      if(show.fit == TRUE){
        
        g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j], size=3) + 
          geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors[j])+
          geom_line(data = real.curve,aes(x=time, y= y),linetype='dashed', color= colors[j])+
          geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors[j]) +
          xlim(c(offset.xlim.negative,Xlim))+
          #   #  ylim(c(0,13))+
          theme_bw() +
          xlab('Days')+
          ylab('Log Viral load')+
          theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
                axis.text.y = element_text(size=16),
                text=element_text(size=16))+
          geom_vline(xintercept = 0, color = 'black')+
          ggtitle(paste0('Household ', k))
      }
      
      if(show.fit == FALSE){
        g =g+ geom_point(data = as.data.frame(obs.data),aes(time,x),color= colors[j], size=3) + 
          geom_line(data = real.curve,aes(x=time, y= y),size=2, color= colors[j])+
          xlim(c(offset.xlim.negative,Xlim))+
          #   #  ylim(c(0,13))+
          theme_bw() +
          xlab('Days')+
          ylab('Log Viral load')+
          theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
                axis.text.y = element_text(size=16),
                text=element_text(size=16))+
          geom_vline(xintercept = 0, color = 'black')+
          ggtitle(paste0('Household ', k))
        
      }
      
      if(show.symptoms)
        g = g+ geom_vline(xintercept = P1$t.symptom, color= colors[j])
    }
    
    for(x in xdata$time[k,xdata$start[k,j]:xdata$end[k,j]]){
      g=g +geom_vline(xintercept = x, color = 'grey', linetype= 'dashed')
    }
    G[[k]] = g
  }
  
  return(G)
  
} 

MCMC_plot_predicted_foi<- function(Chains, params,xdata, Xlim = Tmax){
  
  # load(filename)
  #  Chains= rstan::extract(fit)
  colors_hh = c('#c45441',"#1212ca","#128a18")
  G=NULL
  
  discrete_times = seq(0,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  S = length(Chains$logm)
  h=params$h[1]
  
  for(k in 1:xdata$n_households){
    print(k)
    P1= params %>% 
      filter(household.id == k) %>% 
      mutate(tinfection = t.infection)
    
    foi = array(data=0, dim = c(S, n_discrete_times))
    for (i in 1:S){
      m = Chains$m[i]
      # m =params$m[1]
      for(j in 1:xdata$n_infected[k]){
        mod =   model_gen_triangle(T_P = Chains$T_P[i,k,j], 
                                   T_U = Chains$T_U[i,k,j],
                                   V_P = Chains$V_P[i,k,j],
                                   t = discrete_times,
                                   tinfection = Chains$tinfection[i,k,j]) %>%
          mutate(prob.infection = m*V^P1$h[j] )
        
        foi[i,] = foi[i,] + mod$prob.infection
        
      }
    }
    
    avg = apply(foi,2,mean)
    lower = apply(foi,2,quantile025)
    upper = apply(foi,2,quantile975)
    
    pred.data =  data.frame(time=discrete_times,
                            mean = avg, 
                            quantile025=lower, 
                            quantile975=upper)
    
    
    # the real force of infection
    
    input.foi = 0
    for( j in 1:nrow(P1)){
      X =model_gen_triangle(T_P = P1$T_P[j],
                            T_U = P1$T_U[j], 
                            V_P = P1$V_P[j],  
                            t = discrete_times, 
                            tinfection = P1$t.infection[j] ) %>% 
        mutate(input.foi = P1$m[j]*V^P1$h[j] )
      input.foi = input.foi + X$input.foi
    }
    
    
    g = ggplot() +
      geom_line(data = as.data.frame(pred.data),aes(x=time, y= mean), color= colors_hh[1])+
      geom_ribbon(data = as.data.frame(pred.data), aes(x = time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= colors_hh[1]) +
      geom_line(data = data.frame(time = discrete_times, y = input.foi),aes(x=time, y= y),linetype='dashed', color= colors_hh[1])+
      xlim(c(0,Xlim))+
      #    ylim(c(0,7))+
      theme_bw() +
      xlab('Days')+
      ylab('Force of infection')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+ggtitle(paste0('Household ', k))
    #print(g)
    G[[k]] = g
  }
  
  return(G)
}

MCMC_individual_infection_probability <- function(Chains, params,xdata ){
  
  # load(filename)
  #  Chains= rstan::extract(fit)
  colors_hh = c('#c45441',"#1212ca","#128a18")
  G=NULL
  
  discrete_times = seq(0,Tmax,by=0.1)
  n_discrete_times = length(discrete_times)
  S = length(Chains$logm)
  h=params$h[1]
  m.input = params$m[1]
  
  SAR = array(data=0, dim = c(S, xdata$n_infected_total))
  SAR_input = array(data=0, dim = c(1, xdata$n_infected_total))
  J = 0
  for(k in 1:xdata$n_households){
    #print(k)
    P1= params %>%  filter(household.id == k) %>%  mutate(tinfection = t.infection)
    for(j in 1:xdata$n_infected[k]){
      J=J+1
      SAR_input[J] = 1-exp(-m.input*integral_triangle_model(Tmax,P1[j,],0))
    }
  } 
  
  for (i in 1:S){
    m = Chains$m[i]
    J = 0
    for(k in 1:xdata$n_households){
      for(j in 1:xdata$n_infected[k]){
        J=J+1
        par = data.frame(T_P = Chains$T_P[i,k,j], 
                         T_U = Chains$T_U[i,k,j],
                         V_P = Chains$V_P[i,k,j],
                         h = h)
        SAR[i,J] = 1-exp(-m*integral_triangle_model(100,par,0))
      }
    } 
  }
  
  avg = apply(SAR,2,mean)
  lower = apply(SAR,2,quantile025)
  upper = apply(SAR,2,quantile975)
  
  
  d = data.frame(input = SAR_input[1,], avg = avg, lower=lower, upper = upper) %>% arrange(input) %>% mutate(x =seq(1,xdata$n_infected_total) ) %>% mutate(overestimate )
  
  
  g = ggplot(data=d) +
    geom_segment(aes(x = x, xend= x,y = lower,yend = upper), color="black")+
    geom_point(aes(x = x, y =avg), color="black")+
    geom_point(aes(x=x, y= input), color="red")+
    theme_bw() +
    xlab('Individual') + 
    ylab('Infection probability') +
    theme(axis.text.x = element_text(size=16, angle = 0, vjust=0),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))
  
  print(g)
  return(g)
}

MCMC_plot_generation_interval<- function(Chains,  xdata, Xlim = 20){
  
  S = length(Chains$logm)
  n = 0
  interval = c()
  for(k in 1:xdata$n_households ){
    
    if(xdata$n_infected[k] == 2){
      for (i in 1:S){
        
        interval[n ] = Chains$tinfection[i,k,2]- Chains$tinfection[i,k,1]
        n=n+1
      }
    }  
  }
  print(mean(interval))
  g= data.frame(x = interval) %>%
    ggplot(aes(x=x)) +
    geom_density( colour="darkred",show.legend = FALSE)+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+xlab('Generation interval (days)')
  
  return(g)
}





