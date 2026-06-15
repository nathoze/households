# Effect of the viral load on the probability of transmission
# Figure 2 of Marc et al. Elife

get_pr_trans <- function(prob.transmission, logV,target){
  A=abs(logV-target)
  w = which.min(A)
  if(min(A)>0.5){
    p = NA
  }else{
    p =  prob.transmission[w]
  }
  return(p)
}

# depends on household size

plot_effect_viral_load_on_transmission  <- function(epidemiological.model, fixed.parameters = FALSE, size.dependency, household.size = 5, Tmax = 70, n.simu=1000){
  
  dt = 0.01
  total.time <- seq(0, Tmax, by = dt)
  n.times = length(total.time)
  df.transmission= NULL
  
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
  
  wrapper_transmission <- function(i){
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
    
    list.transmission <- generate_transmission_individual_based(HH.parameters, size.dependency, Tmax= Tmax)
    events  =cbind(list.transmission$all.events, m=i)
    # if(is.null(list.transmission$viral.load.at.infection)){
    #   d =NULL
    # }else{
    #   d = data.frame(m = i,
    #                  viral.load.at.infection = list.transmission$viral.load.at.infection)
    # }
    
    return(events)
  }
  
  #n.simu = 600
  indices = rep(c(1,2,3,4), each = n.simu)
  #indices = rep(c(1,2), each = n.simu)
  
  A = indices %>%
    map(wrapper_transmission) %>%
    bind_rows() %>%
    # mutate(logV = log10(viral.load.at.infection)) %>%
    mutate(logV = log10(V)) %>%
    filter(!is.na(logV)) %>% 
    group_by(m) %>%
    # group_by(infection) %>%
    # mutate(p3= length(which(((logV-4)<0.0)& infection==0))) %>%
    #  mutate(q3= length(which(((logV-4)<0.0)& infection==1)))  %>%
    mutate(p3= length(which((abs(logV-3.5)<0.5)& infection==0))) %>%
    mutate(q3= length(which((abs(logV-3.5)<0.5)& infection==1)))  %>%
    mutate(p4= length(which((abs(logV-4.5)<0.5)& infection==0))) %>%
    mutate(q4= length(which((abs(logV-4.5)<0.5)& infection==1)))  %>%
    mutate(p5= length(which((abs(logV-5.5)<0.5)& infection==0))) %>%
    mutate(q5= length(which((abs(logV-5.5)<0.5)& infection==1)))  %>%
    mutate(p6= length(which((abs(logV-6.5)<0.5)& infection==0))) %>%
    mutate(q6= length(which((abs(logV-6.5)<0.5)& infection==1)))  %>%
    mutate(p7= length(which(( (logV-7)>0)& infection==0)))  %>%
    mutate(q7= length(which(( (logV-7)>0)& infection==1)))  %>% 
    select(m,p3,q3,p4,q4,p5,q5,p6,q6,p7,q7) %>%
    mutate("3-4"=q3/(p3+q3) )%>%
    mutate("4-5"=q4/(p4+q4)) %>%
    mutate("5-6"=q5/(p5+q5) )%>%
    mutate("6-7"=q6/(p6+q6) )%>%
    mutate(">7"=q7/(p7+q7) ) %>% distinct( .keep_all = TRUE)
  
  
  B =  A %>% select(m,"3-4","4-5","5-6","6-7",">7") %>%
    pivot_longer(cols=c("3-4","4-5","5-6","6-7",">7"),
                 values_to='value') %>%
    mutate( m.width = ifelse(m ==1 | m==4, "sd = 0","sd = 0.85")) %>%
    mutate( m.value = ifelse(m ==1 | m==3, mean.large, mean.small))   %>% 
    mutate(N2 = factor(name, levels =c("3-4","4-5","5-6","6-7",">7") ) ) %>%
    mutate(m = as.factor(m)) %>%
    mutate(m.width = as.factor(m.width))%>%
    mutate(m.value = as.factor(m.value)) %>%
    mutate(m.value = relevel(m.value, as.character(mean.small)))
  
  B = B%>% ggplot()+
    geom_col(aes(x = N2, y=value, fill = m) , colour="black",show.legend = FALSE)+
    facet_grid(cols = vars(m.width), rows=vars(m.value), scales = "free") +
    ylab('Mean probability of transmission') +
    xlab('Viral load (log10)')+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16)) + 
    scale_fill_manual(values = colors.m.values)
  print(B)
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/Probability_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent.pdf") , width = 7,  height = 6)
    dev.off() 
    dev.copy(png,paste0("results/Probability_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent.png") , width = 7,  height = 6, units="in", res=200)
    dev.off() 
  }
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/Probability_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 7,  height = 6)
    dev.off() 
    dev.copy(png,paste0("results/Probability_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 7,  height = 6, units="in", res=200)
    dev.off() 
  }
  
  # return(df.transmission)
} 


