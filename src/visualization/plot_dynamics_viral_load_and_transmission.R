
plot_dynamics_viral_load_and_transmission <- function(epidemiological.model, fixed.parameters = FALSE ,  Tmax = 70){ 
  
  dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  n.times = length(total.time)
  
  m.model =1 
  
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
  
  map_trajectories<-function(i){
    
    param.list = all.parameters[[m.model]]
    par =param.list[i,]
#    all.titers <- get_viral_load_transmission_probability(par, delay = par$incubation, epid.model = epidemiological.model, Tmax) 
    all.titers <- get_viral_load_transmission_probability(par, delay = 0, epid.model = epidemiological.model, Tmax) 
    
    df = data.frame(i = i, t = all.titers$t, VI = all.titers$VI, prob.transmission  = all.titers$prob.transmission)   
    return(df)
    
  } 
  
  map_trajectories_no_delay <-function(i){
    
    param.list = all.parameters[[m.model]]
    par =param.list[i,]
    all.titers <- get_viral_load_transmission_probability(par, delay = 0, epid.model = epidemiological.model, Tmax)
    #%>% 
    # filter(t>=par$incubation)
    
    df = data.frame(i = i, t = all.titers$t, VI = all.titers$VI, prob.transmission  = all.titers$prob.transmission)   
    return(df)
    
  } 
  
  indices = seq(1,1000)
  
  viral.load = indices %>%
    map(map_trajectories_no_delay) %>%
    bind_rows() %>%
    mutate(logVI = log10(VI)) %>%
    group_by(t) %>%
    summarise_at(.vars = "logVI",
                 .funs = c(mean="mean",lower = "quantile05", upper="quantile95"))
  
  
  viral.load.NA = viral.load %>%
    mutate(mean.viral.load = ifelse(mean>2, mean, 2)) %>%
    mutate(lower.viral.load = ifelse(lower>2, lower, 2)) %>% 
    mutate(upper.viral.load = ifelse(upper>2, upper, 2))
  
  P=NULL
  for(m.model in 1:4){
    P[[m.model]] = indices %>%
      map(map_trajectories_no_delay) %>%
      bind_rows() %>%
      group_by(t) %>%
      summarise_at(.vars = "prob.transmission",
                   .funs = c(mean.p="mean",median.p="median",lower.p = "quantile05", upper.p="quantile95")) %>% 
      mutate(m.model = m.model)
    
  }
  Q = rbind(P[[1]],P[[2]],P[[3]],P[[4]])
  
  
  m2 = max(viral.load$mean)
  w =which.max(viral.load$mean)
  
  factor.probability.y=1
  Q2 = Q %>%
    left_join(viral.load.NA , by = c( "t")) %>%
    group_by(m.model) %>%
    mutate( m.width = ifelse(m.model ==1 | m.model==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m.model ==1 | m.model==3, mean.large, mean.small))  %>%
    mutate(m = as.factor(m.model)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small))) %>%
    mutate(factor.probability.y = max(mean.viral.load)/max(upper.p)) %>%
    filter(lower>0) %>% 
    ggplot()+
    geom_ribbon(aes(x = t, ymin=lower.viral.load, ymax=upper.viral.load), alpha=0.2, fill= "black",show.legend = FALSE) +
    geom_line(aes(x= t, y= mean.viral.load),linewidth =1.3, color= "black",show.legend = FALSE)+
    geom_ribbon(aes(x = t, ymin= factor.probability.y*lower.p, ymax = factor.probability.y*upper.p, fill = m, color = m), alpha=0.8,show.legend = FALSE) +
    geom_line(aes(x= t, y= factor.probability.y*median.p, color = m),linetype = 'dashed', linewidth =1.1,show.legend = FALSE)+
    facet_grid(cols = vars(m.width), rows=vars(m.value)) +
    scale_fill_manual(values = colors.m.values)+
    scale_color_manual(values = colors.m.values.dark)+
    scale_y_continuous( "Viral load (log10)", 
                        sec.axis = sec_axis(~ . /factor.probability.y , name = "Transmission probability")    )+
   # xlab('Time since symptom onset (days)') +
    xlab('Days post infection') +
    ylab('Viral load (log10)')+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))
  
  print(Q2)
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/CompareProba_",epidemiological.model,"_fixed_parameters.pdf"), width = 7,  height = 6)
    dev.off()
  }
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/CompareProba_",epidemiological.model,".pdf"), width = 7,  height = 6)
    dev.off()
  }
  
  
  
  
}


