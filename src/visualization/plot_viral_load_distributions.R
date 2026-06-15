
plot_viral_load_distribution <- function(epidemiological.model, size.dependency, fixed.parameters = FALSE, Tmax= 70, household.size = 5 , n.simu =1000){
  #    Tmax = 70
  #dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  
  n.times = length(total.time)
  integer.times.index = seq(1,n.times, by = round(1/dt))
  
  M = get_m_values(epidemiological.model)
  mean.large = M$mean.large
  mean.small = M$mean.small
  
  all.parameters = NULL
  if(fixed.parameters == FALSE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
    }
  }   
  if(fixed.parameters == TRUE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }   
  
  wrapper_VL_distribution <- function(i){
    print(i)
    if(i ==1){
      m.name =   "large_narrow_m"
    }
    if(i ==2){
      m.name =   "small_wide_m"
    }
    if(i ==3){
      m.name =   "large_wide_m"
    }
    if(i ==4){
      m.name =   "small_narrow_m"
    }
    
    param.list = all.parameters[[i]]
    hh.index = sample(param.list$id,size=household.size,replace = FALSE)
    
    HH.parameters  = param.list[hh.index,]
    HH.parameters$status = 'contact'
    HH.parameters$status[1] = 'index' 
    
    HH.parameters$t.infection = -1
    HH.parameters$t.infection[1] = 0
    
    list.transmission <- generate_household_transmission(HH.parameters, size.dependency = size.dependency, Tmax = Tmax)
    logVL = list.transmission$logVL %>% 
      filter(t %in% seq(0,Tmax))
    
    d=data.frame(m = i, t= logVL$t, logVL = logVL$logVL, VL = 10^logVL$logVL)
    
    d  = d%>%
      group_by(m,t) %>%
      mutate(sumlogVL = sum(logVL)) %>%
      mutate(sumVL = 10^sumlogVL)
    return(d)
    
  }
  
  indices = rep(c(1,2,3,4), each = n.simu)
  
  VL_distribution = indices %>%
    map(wrapper_VL_distribution) %>%
    bind_rows() %>%
    group_by(t,m)
  
  B = VL_distribution %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) 
  
  VL_distribution_observed = B %>% filter(logVL >=2)
  VL_distribution_observed = B %>% mutate(logVL = ifelse(logVL<2,2,logVL))
  VL_distribution_mean = VL_distribution_observed %>% mutate(M = 10^mean(logVL))
  
  
  g = VL_distribution_observed%>% 
    sample_frac(0.25) %>%
    ggplot() + 
    geom_jitter(aes(x=t ,y =logVL), alpha = 0.03) + 
    ylim(c(2,NA)) + 
    geom_line(data=VL_distribution_mean, aes(x=t, y = log10(M), color = m), linewidth=1.2, show.legend = FALSE)+ 
    ylim(c(2,5))+
    # scale_y_log10(limits = c(10,NA))+
    facet_grid(cols = vars(m.width), rows=vars(m.value)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    xlab("Days")+
    ylab('Viral load (log10)')+  
    scale_color_manual(values = colors.m.values.dark) 
  
  print(g)
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.png"), width = 12,  height =10, units="in",res=200)
    dev.off()
  }   
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 12,  height =10, units="in",res= 200)
    dev.off()
  }   
  
  
}

plot_viral_load_distribution_mean_viral_load <- function(epidemiological.model, size.dependency,fixed.parameters = FALSE, Tmax= 70, household.size = 5 , n.simu =1000 ){
  #    Tmax = 70
  #dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  
  n.times = length(total.time)
  integer.times.index = seq(1,n.times, by = round(1/dt))
  M = get_m_values(epidemiological.model)
  mean.large = M$mean.large
  mean.small = M$mean.small
  
  all.parameters = NULL
  if(fixed.parameters == FALSE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
    }
  }   
  if(fixed.parameters == TRUE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }   
  
  wrapper_VL_distribution <- function(i){
    print(i)
    if(i ==1){
      m.name =   "large_narrow_m"
    }
    if(i ==2){
      m.name =   "small_wide_m"
    }
    if(i ==3){
      m.name =   "large_wide_m"
    }
    if(i ==4){
      m.name =   "small_narrow_m"
    }
    
    param.list = all.parameters[[i]]
    hh.index = sample(param.list$id,size=household.size,replace = FALSE)
    
    HH.parameters  = param.list[hh.index,]
    HH.parameters$status = 'contact'
    HH.parameters$status[1] = 'index' 
    
    HH.parameters$t.infection = -1
    HH.parameters$t.infection[1] = 0
    
    list.transmission <- generate_household_transmission(HH.parameters, size.dependency = size.dependency, Tmax = Tmax)
    
    n.inf= list.transmission$n.infection
    
    logVL = list.transmission$logVL %>% 
      filter(t %in% seq(0,Tmax)) %>% 
      group_by(t) %>%
      mutate(logVL = ifelse(logVL<detection.limit,detection.limit, logVL)) %>%
      mutate(meanlogVL = mean(logVL) ) %>%
      mutate(weighted.meanlogVL = ( n.inf*mean(logVL) + (household.size-n.inf)*detection.limit)/household.size  )  
    
    d=data.frame(m = i, t= logVL$t, logVL = logVL$weighted.meanlogVL, VL = 10^logVL$weighted.meanlogVL)
    
    return(d)
    
  }
  indices = rep(c(1,2,3,4), each = n.simu)
  
  VL_distribution = indices %>%
    map(wrapper_VL_distribution) %>%
    bind_rows() %>%
    group_by(t,m)
  
  VL_distribution_observed = VL_distribution %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) %>% 
    mutate(logVL = ifelse(logVL<detection.limit,detection.limit,logVL))%>% 
    mutate(M = 10^mean(logVL))
  
  #  VL_distribution_observed = B %>% filter(logVL >= detection.limit )
  # VL_distribution_observed = B %>% mutate(logVL = ifelse(logVL<detection.limit,detection.limit,logVL))
  #  VL_distribution_mean = VL_distribution_observed %>% mutate(M = 10^mean(logVL))
  
  g = VL_distribution_observed%>%   
    sample_frac(0.25) %>%
    ggplot() + 
    geom_jitter(aes(x=t ,y = logVL), alpha = 0.03) + 
    ylim(c(2,NA)) + 
    geom_line(aes(x=t, y = log10(M), color=m), linewidth=1.2, show.legend = FALSE)+ 
    ylim(c(2,5))+
    #scale_y_log10(limits = c(2,4.5))+
    facet_grid(cols = vars(m.width), rows=vars(m.value)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    xlab("Days")+
    ylab('Viral load (log10)')+  
    scale_color_manual(values = colors.m.values.dark)
  
  print(g)
  # dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent",".pdf"), width = 12,  height =10)
  
  
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.png"), width = 12,  height =10, units="in", res=200)  
    dev.off()
  }   
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 12,  height =10, units="in", res=200)
    dev.off()
  }   
  
}

plot_viral_load_distribution_envelop <- function(epidemiological.model, size.dependency, fixed.parameters = FALSE, Tmax= 70, household.size = 5 , n.simu =1000){
  #    Tmax = 70
  #dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  
  n.times = length(total.time)
  integer.times.index = seq(1,n.times, by = round(1/dt))
  
  M = get_m_values(epidemiological.model)
  mean.large = M$mean.large
  mean.small = M$mean.small
  
  all.parameters = NULL
  if(fixed.parameters == FALSE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
    }
  }   
  if(fixed.parameters == TRUE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }   
  
  wrapper_VL_distribution <- function(i){
    print(i)
    if(i ==1){
      m.name =   "large_narrow_m"
    }
    if(i ==2){
      m.name =   "small_wide_m"
    }
    if(i ==3){
      m.name =   "large_wide_m"
    }
    if(i ==4){
      m.name =   "small_narrow_m"
    }
    
    param.list = all.parameters[[i]]
    hh.index = sample(param.list$id,size=household.size,replace = FALSE)
    
    HH.parameters  = param.list[hh.index,]
    HH.parameters$status = 'contact'
    HH.parameters$status[1] = 'index' 
    
    HH.parameters$t.infection = -1
    HH.parameters$t.infection[1] = 0
    
    list.transmission <- generate_household_transmission(HH.parameters, size.dependency = size.dependency, Tmax = Tmax)
    logVL = list.transmission$logVL %>% 
      filter(t %in% seq(0,Tmax))
    
    d=data.frame(m = i, t= logVL$t, logVL = logVL$logVL, VL = 10^logVL$logVL)
    
    d  = d%>%
      group_by(m,t) %>%
      mutate(sumlogVL = sum(logVL)) %>%
      mutate(sumVL = 10^sumlogVL)
    return(d)
    
  }
  
  indices = rep(c(1,2,3,4), each = n.simu)
  
  VL_distribution = indices %>%
    map(wrapper_VL_distribution) %>%
    bind_rows() %>%
    group_by(t,m)
  
  B = VL_distribution %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) 
  
  VL_distribution_observed = B %>% filter(logVL >=2)
  VL_distribution_observed = B %>% mutate(logVL = ifelse(logVL<2,2,logVL))  %>%
    mutate(M = 10^mean(logVL)) %>%
    mutate(lower =  10^quantile05(logVL))   %>%
    mutate(upper = 10^quantile95(logVL))
  
  
  g = VL_distribution_observed%>% 
    sample_frac(0.25) %>%
    ggplot() + 
    geom_jitter(aes(x=t ,y =logVL), alpha = 0.03) + 
    ylim(c(2,NA)) + 
    geom_line(aes(x=t, y = log10(M), color = m), linewidth=1.2, show.legend = FALSE)+ 
    geom_ribbon(aes(x = t, ymin= log10(lower), ymax = log10(upper), fill = m, color = m), alpha=0.7,show.legend = FALSE) +
    # ylim(c(2,5))+
    # scale_y_log10(limits = c(10,NA))+
    facet_grid(cols = vars(m.width), rows=vars(m.value)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    xlab("Days")+
    ylab('Viral load (log10)')+  
    scale_color_manual(values = colors.m.values.dark) +
    scale_fill_manual(values = colors.m.values) 
  
  print(g)
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.png"), width = 12,  height =10, units="in",res=200)
    dev.off()
  }   
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 12,  height =10, units="in",res= 200)
    dev.off()
  }   
}

plot_viral_load_distribution_mean_envelop <- function(epidemiological.model, size.dependency,fixed.parameters = FALSE, Tmax= 70, household.size = 5 , n.simu =1000 ){
  #    Tmax = 70
  #dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  
  n.times = length(total.time)
  integer.times.index = seq(1,n.times, by = round(1/dt))
  M = get_m_values(epidemiological.model)
  mean.large = M$mean.large
  mean.small = M$mean.small
  
  all.parameters = NULL
  if(fixed.parameters == FALSE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,".csv"))
    }
  }   
  if(fixed.parameters == TRUE){
    for(m.index in 1:4){
      all.parameters[[m.index]] = read.csv( paste0("parameters/parameter_list_",epidemiological.model,"_m",m.index,"_fixed_parameters.csv"))
    }
  }   
  
  wrapper_VL_distribution <- function(i){
    print(i)
    if(i ==1){
      m.name =   "large_narrow_m"
    }
    if(i ==2){
      m.name =   "small_wide_m"
    }
    if(i ==3){
      m.name =   "large_wide_m"
    }
    if(i ==4){
      m.name =   "small_narrow_m"
    }
    
    param.list = all.parameters[[i]]
    hh.index = sample(param.list$id,size=household.size,replace = FALSE)
    
    HH.parameters  = param.list[hh.index,]
    HH.parameters$status = 'contact'
    HH.parameters$status[1] = 'index' 
    
    HH.parameters$t.infection = -1
    HH.parameters$t.infection[1] = 0
    
    list.transmission <- generate_household_transmission(HH.parameters, size.dependency = size.dependency, Tmax = Tmax)
    
    n.inf= list.transmission$n.infection
    
    logVL = list.transmission$logVL %>% 
      filter(t %in% seq(0,Tmax)) %>% 
      group_by(t) %>%
      mutate(logVL = ifelse(logVL<detection.limit,detection.limit, logVL)) %>%
      mutate(meanlogVL = mean(logVL) ) %>%
      mutate(weighted.meanlogVL = ( n.inf*mean(logVL) + (household.size-n.inf)*detection.limit)/household.size  )  
    
    d=data.frame(m = i, t= logVL$t, logVL = logVL$weighted.meanlogVL, VL = 10^logVL$weighted.meanlogVL)
    
    return(d)
    
  }
  indices = rep(c(1,2,3,4), each = n.simu)
  
  VL_distribution = indices %>%
    map(wrapper_VL_distribution) %>%
    bind_rows() %>%
    group_by(t,m)
  
  VL_distribution_observed = VL_distribution %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) %>% 
    mutate(logVL = ifelse(logVL<detection.limit,detection.limit,logVL))%>% 
    mutate(M = 10^mean(logVL)) %>%
    mutate(lower =  10^quantile05(logVL))   %>%
    mutate(upper = 10^quantile95(logVL))
  
  
  g = VL_distribution_observed%>%   
    sample_frac(0.25) %>%
    ggplot() + 
    geom_jitter(aes(x=t ,y = logVL), alpha = 0.03) + 
    geom_line(aes(x=t, y = log10(M), color=m), linewidth=1.2, show.legend = FALSE)+ 
    geom_ribbon(aes(x = t, ymin= log10(lower), ymax = log10(upper), fill = m, color = m), alpha=0.7,show.legend = FALSE) +
    #ylim(c(2,NA)) + 
    ylim(c(2,5))+
    #scale_y_log10(limits = c(2,4.5))+
    facet_grid(cols = vars(m.width), rows=vars(m.value)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    xlab("Days")+
    ylab('Viral load (log10)')+  
    scale_color_manual(values = colors.m.values.dark)+
    scale_fill_manual(values = colors.m.values)
  
  print(g)
  # dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent",".pdf"), width = 12,  height =10)
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_sum_viral_load_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent.png"), width = 12,  height =10, units="in", res=200)  
    dev.off()
  }   
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/ViralLoadDistributions_sum_viral_load_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 12,  height =10)
    dev.off()
    dev.copy(png,paste0("results/ViralLoadDistributions_sum_viral_load_envelop_",epidemiological.model, "_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 12,  height =10, units="in", res=200)
    dev.off()
  }   
  
}