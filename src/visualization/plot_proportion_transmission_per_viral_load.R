# Effect of the viral load on the probability of transmission
# Figure 2 of Marc et al. Elife

get_pr_trans <- function(prob.transmission, logVI,target){
  A=abs(logVI-target)
  w = which.min(A)
  if(min(A)>0.5){
    p = NA
  }else{
    p =  prob.transmission[w]
  }
  return(p)
}

# depends on household size

plot_proportion_transmission_per_viral_load <- function(epidemiological.model, fixed.parameters = FALSE, size.dependency, household.size = 5, Tmax = 70,  n.simu = 600){
  
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
  
  indices = rep(c(1,2,3,4), each = n.simu)
  #indices = rep(c(1,2), each = n.simu)
  
  A = indices %>%
    map(wrapper_transmission) %>%
    bind_rows() %>%
    # mutate(logVI = log10(viral.load.at.infection)) %>%
    mutate(logVI = log10(VI)) %>%
    filter(!is.na(logVI)) %>% 
    filter(infection==1) %>%
    group_by(m) %>%
    # group_by(infection) %>%
    mutate(q3= length(which((abs(logVI-3.5)<0.5))))  %>%
    mutate(q4= length(which((abs(logVI-4.5)<0.5))))  %>%
    mutate(q5= length(which((abs(logVI-5.5)<0.5))))  %>%
    mutate(q6= length(which((abs(logVI-6.5)<0.5))))  %>%
    mutate(q7= length(which((logVI-7)>0)))  %>% 
    select(m,q3,q4,q5,q6,q7) %>%
    mutate("3-4"=q3/(q3+q4+q5+q6+q7) )%>%
    mutate("4-5"=q4/(q3+q4+q5+q6+q7)) %>%
    mutate("5-6"=q5/(q3+q4+q5+q6+q7))%>%
    mutate("6-7"=q6/(q3+q4+q5+q6+q7) )%>%
    mutate(">7"=q7/(q3+q4+q5+q6+q7) ) %>% distinct( .keep_all = TRUE)
  
  
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
    ylim(c(0,0.8))+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16)) + 
    scale_fill_manual(values = colors.m.values)
  print(B)
  
  if(fixed.parameters == FALSE){
    dev.copy(pdf,paste0("results/Transmission_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent.pdf") , width = 7,  height = 6)
    dev.off() 
    dev.copy(png,paste0("results/Transmission_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent.png") , width = 7,  height = 6, units="in", res=200)
    dev.off() 
  }
  if(fixed.parameters == TRUE){
    dev.copy(pdf,paste0("results/Transmission_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.pdf"), width = 7,  height = 6)
    dev.off() 
    dev.copy(png,paste0("results/Transmission_ViralLoad_",epidemiological.model,"_HHsize_",household.size,"_", size.dependency,"_dependent_fixed_parameters.png"), width = 7,  height = 6, units="in", res=200)
    dev.off() 
  }
  
  # return(df.transmission)
} 


