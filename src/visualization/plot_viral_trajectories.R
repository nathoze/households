

## plot individual trajectories of viral load -


 

plot_viral_trajectories <- function(model.virus = "triangle", Tmax = 70){ 
  epidemiological.model = "doseresponse" # doesn't matter
  dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  n.times = length(total.time)
  
  all.parameters_not_fixed = NULL
  all.parameters_fixed = NULL

  if(model.virus == "triangle"){
    for(m.index in 1:4){
      all.parameters_not_fixed[[m.index]] = read.csv( paste0("../parameters/parameter_list_triangle_",epidemiological.model,"_m",m.index,".csv"))
      all.parameters_fixed[[m.index]] = read.csv( paste0("../parameters/parameter_list_triangle_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }
  
  if(model.virus == "ode"){
    for(m.index in 1:4){
      all.parameters_not_fixed[[m.index]] = read.csv( paste0("../parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
      all.parameters_fixed[[m.index]] = read.csv( paste0("../parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }
 
  wrapper_trajectory<-function(i){
    
    param.list = all.parameters_not_fixed[[1]]
    par =param.list[i,]
  #  all.titers <- get_viral_load_transmission_probability(par, delay = par$incubation, epid.model = epidemiological.model, Tmax) %>% filter(t>=par$incubation)
    all.titers <- get_viral_load_transmission_probability(par, delay =0, epid.model = epidemiological.model, Tmax, dt = 0.01) 
    
    df = data.frame(i = i, t = all.titers$t, VI = all.titers$VI + all.titers$VNI, par=par)   
    return(df)
    
  } 

  param.list = all.parameters_fixed[[1]]
  par =param.list[1,]
  all.titers <- get_viral_load_transmission_probability(par, delay =0, epid.model = epidemiological.model, Tmax, dt = 0.01) #%>% filter(t>=par$incubation)
 # all.titers <- get_viral_load_transmission_probability(par, delay = par$incubation, epid.model = epidemiological.model, Tmax) %>% filter(t>=par$incubation)
  df = data.frame(i = 1, t = all.titers$t, V = all.titers$VI + all.titers$VNI, VI = all.titers$VI)#, par=par)   

  
  indices = seq(1,200)
  A = indices %>%
    map(wrapper_trajectory) %>%
    bind_rows() %>%
    ggplot()+
    geom_line(aes(x=t, y = log10(VI), group = i), alpha=0.12) +
    geom_line(data = df, aes(x=t, y = log10(V), group = i), linewidth = 1.2, color='red') +
    #geom_line(data = df, aes(x=t, y = log10(V), group = i), linewidth = 1.2, color='blue') +
 #  geom_line(data = df, aes(x=t, y = log10(V)-log10(VI), group = i), linewidth = 1.2, color='blue') +
    ylim(c(6,NA))+
    #  scale_y_log10(limits = c(0.1,NA))+
    xlab('Days') +
    ylab('Log Viral load')+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    xlim(c(0,25))
  print(A)
  
  dev.copy(pdf,paste0("../Results/Viral_load_trajectories.pdf"), width = 7,  height = 4)
  dev.off()
  
}






plot_probability_trajectories <- function(epidemiological.model, fixed.parameters = FALSE ,  Tmax = 70){ 
  
  dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  n.times = length(total.time)
  
  
  M = get_m_values(epidemiological.model)
  mean.large = M$mean.large
  mean.small = M$mean.small
  
  all.parameters = NULL
  if(fixed.parameters == FALSE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("../parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
    }
  }   
  if(fixed.parameters == TRUE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("../parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }   
  
  #  wrapper_trajectory<-function(i){
  #    
  #    param.list = all.parameters[[1]]
  #    par =param.list[i,]
  #    all.titers <- get_viral_load_transmission_probability(par, delay = par$incubation, epid.model = epidemiological.model, Tmax) %>% filter(t>=par$incubation)
  #    
  #    df = data.frame(i = i, t = all.titers$t, VI = all.titers$VI, par=par)   
  #    return(df)
  #    
  #  } 
  # 
  #  
  # indices = seq(1,100)
  #  A = indices %>%
  #    map(wrapper_trajectory) %>%
  #    bind_rows() %>%
  #    ggplot()+
  #    geom_line(aes(x=t, y = VI, group = i), alpha=0.2) +
  #    scale_y_log10(limits = c(0.1,NA))+
  #    xlab('Days') +
  #    ylab('Viral load')+
  #    theme_bw()+
  #    theme(axis.text.x = element_text(size=16),
  #          axis.text.y = element_text(size=16),
  #          text=element_text(size=16))
  #  print(A)
  #  #dev.copy(pdf,paste0("results/Viral_load_trajectories.pdf"), width = 7,  height = 6)
  #  
  #  if(fixed.parameters == TRUE){
  #    dev.copy(pdf,paste0("results/Viral_load_trajectories_fixed_parameters.pdf"), width = 7,  height = 6)
  #  } 
  #  if(fixed.parameters == FALSE){
  #    dev.copy(pdf,paste0("results/Viral_load_trajectories.pdf"), width = 7,  height = 6)
  #  }
  #  
  #  dev.off()
  #  
  #  
  #  
  #  
  #  
  wrapper_probability<-function(i){
    v=value[i]
    if(v ==1){
      m.name =   "large_narrow_m"
    }
    if(v ==2){
      m.name =   "small_wide_m"
    }
    if(v ==3){
      m.name =   "large_wide_m"
    }
    if(v ==4){
      m.name =   "small_narrow_m"
    }
    
    param.list = all.parameters[[v]]
    k = sample(param.list$id,size = 1)
    
    par =param.list[k,]
   #all.titers <- get_viral_load_transmission_probability(par, delay = par$incubation, epid.model = epidemiological.model, Tmax)
     all.titers <- get_viral_load_transmission_probability(par, delay = 0, epid.model = epidemiological.model, Tmax)
    df = data.frame(i = i, m = v, t = all.titers$t, VI = all.titers$VI, 
                    prob.transmission = all.titers$prob.transmission,
                    total.proba.transmission = all.titers$total.prob.transmission )   
    return(df)
    
  }
  
  n.simu=50
  
  indices = seq(1,4*n.simu)
  value =  rep(c(1,2,3,4), each = n.simu)
  
  
  plot.proba = indices %>%
    map(wrapper_probability) %>%
    bind_rows()  %>%
    group_by(m) %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) %>%
    ggplot() +
    geom_line(aes(x=t, y = total.proba.transmission,  colour=m, group = i), alpha=0.4, show.legend = FALSE) +
    facet_grid(cols = vars(m.width), rows=vars(m.value), scales = "free") +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+ ylab("Cumulative transmission probability") +
    xlab("Days")+
    scale_color_manual(values = colors.m.values)
  
  print(plot.proba)
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/Cumulative_transmission_",epidemiological.model, "_fixed_parameters.pdf"), width = 7,  height = 6)
    
  } 
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/Cumulative_transmission_",epidemiological.model, ".pdf"), width = 7,  height = 6)
  }
  dev.off()
  
  n.simu=50
  
  indices = seq(1,4*n.simu)
  value =  rep(c(1,2,3,4), each = n.simu)
  
  plot.proba = indices %>%
    map(wrapper_probability) %>%
    bind_rows()  %>%
    group_by(m) %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) %>%
    ggplot() +
    geom_line(aes(x=t, y = prob.transmission,  colour=m, group = i), alpha=0.4, show.legend = FALSE) +
    facet_grid(cols = vars(m.width), rows=vars(m.value), scales="free") +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+ 
    ylab("Transmission probability") +
    xlab("Days")+
    scale_color_manual(values = colors.m.values)  
  
  print(plot.proba)
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("../Results/Transmission_",epidemiological.model, "_fixed_parameters.pdf"), width = 7,  height = 6)
    
  } 
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("../Results/Transmission_",epidemiological.model, ".pdf"), width = 7,  height = 6)
  }
  dev.off()
  
}


