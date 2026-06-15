MCMC_compare_input_estimates_without_symptoms <- function(Chains, xdata, lims = c(25,10,15,15,15)){
  
  # load(filename)
  # Chains= rstan::extract(fit)
  plot_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    # print(X.estimated.lower)
    # print(X.estimated.upper)
    df.times= data.frame(X.estimated = X.estimated, X.input=X.input)
    g = ggplot() +geom_point(data=df.times,aes(y=X.estimated, x = X.input),alpha=0.3)+
      geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower,yend = X.estimated.upper),alpha=0.4)+
      theme_bw() +
      geom_abline(slope=1, intercept = 0, linetype='dashed')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+
      xlab(paste0(label, ' - Input'))+
      ylab(paste0(label, ' - Posterior estimate'))+
      xlim(c(xmin,lim))+ylim(c(xmin,lim))#+coord_equal()
    
    return(g)
  }
  
  G = list()
  P = xdata$all.HH %>% mutate(tinfection = t.infection) 
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection","T_P",'T_U',"V_P" )){
    C = Chains[[var]]
    J=J+1
    label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load')[J]
    
    lim = lims[J]
    G[[J]]=  plot_compare(var)
    
  }
  return(G)
  # grid.arrange(G[[1]],G[[2]],G[[3]],G[[4]],
  #              ncol =2, nrow = 2)
  
  # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
  # dev.off()
}

# 
# MCMC_compare_input_estimates_residuals <- function(Chains, xdata, lims = c(3,3,3,3,5)){
#   
#   # load(filename)
#   # Chains= rstan::extract(fit)
#   plot_compare <- function(var){
#     xmin = 0
#     if(var == "tinfection")
#       xmin = -9
#     if(var == "TS")
#       xmin = -2
#     
#     X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
#     for(k in 1:xdata$n_households ){
#       for(j in 1:xdata$n_infected[k]){
#         
#         P1= P%>% filter(household.id == k)
#         P1=P1[j,]
#         X.estimated = c(X.estimated, median(C[,k,j]))
#         X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
#         X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
#         X.input = c(X.input,  P1[[var]])
#       }
#     }
#     # print(X.estimated.lower)
#     # print(X.estimated.upper)
#     df.times= data.frame(X.estimated = X.estimated, X.input=X.input)
#     g = ggplot() +    geom_hline(yintercept = 0, linetype='dashed')+
#       geom_point(data=df.times,aes(y=X.estimated-X.input, x = X.input),alpha=0.2)+
#       geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower-X.input,yend = X.estimated.upper-X.input),alpha=0.2)+
#       theme_bw() +
#       #  geom_abline(slope=1, intercept = 0, linetype='dashed')+
#       theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
#             axis.text.y = element_text(size=16),
#             text=element_text(size=16))+
#       xlab(paste0(label, ' - Input'))+
#       ylab(paste0(label, ' - Posterior estimate'))+
#       # xlim(c(NA,NA))+
#       ylim(c(-lim,lim))#+coord_equal()
#     
#     return(g)
#   }
#   
#   G = list()
#   P = xdata$all.HH %>% mutate(tinfection = t.infection) 
#   J=0
#   #for(var in c("tinfection","slope_1","slope_2","time_up" )){
#   for(var in c("tinfection","T_P",'T_U',"V_P" )){
#     C = Chains[[var]]
#     J=J+1
#     label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load')[J]
#     
#     lim = lims[J]
#     G[[J]]=  plot_compare(var)
#     
#   }
#   return(G)
#   # grid.arrange(G[[1]],G[[2]],G[[3]],G[[4]],
#   #              ncol =2, nrow = 2)
#   
#   # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
#   # dev.off()
# }
# 




MCMC_compare_input_estimates_residuals <- function(Chains, xdata, lims = c(3,3,3,3,5)){
  
  
  color1= 'red'
  color2 = 'black'
  # load(filename)
  # Chains= rstan::extract(fit)
  plot_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    # print(X.estimated.lower)
    # print(X.estimated.upper)
    df.times= data.frame(X.estimated = X.estimated, X.input=X.input, status = xdata$all.HH$status) %>% mutate(status = as.factor(status))
    
    
    g = ggplot() +    geom_hline(yintercept = 0, linetype='dashed')+
      geom_point(data=df.times,aes(y=X.estimated-X.input, x = X.input, color= status),alpha=0.4)+
      geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower-X.input,yend = X.estimated.upper-X.input, color= status),alpha=0.4)+
      theme_bw() +
      #  geom_abline(slope=1, intercept = 0, linetype='dashed')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+
      xlab(paste0(label, ' - Input'))+
      ylab(paste0(label, ' - Posterior estimate'))+
      scale_color_brewer(palette = "Set1")+
      # xlim(c(NA,NA))+
      ylim(c(-lim,lim))#+coord_equal()
    
    return(g)
  }
  
  G = list()
  P = xdata$all.HH %>% mutate(tinfection = t.infection) 
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection","T_P",'T_U',"V_P" )){
    C = Chains[[var]]
    J=J+1
    label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load')[J]
    
    lim = lims[J]
    G[[J]]=  plot_compare(var)
    
  }
  return(G)
  # grid.arrange(G[[1]],G[[2]],G[[3]],G[[4]],
  #              ncol =2, nrow = 2)
  
  # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
  # dev.off()
}















MCMC_compare_summary_statistics <- function(Chains, xdata){
  
  # load(filename)
  # Chains= rstan::extract(fit)
  summary_statistics_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    
    df.times= data.frame(X.estimated = X.estimated, 
                         X.estimated.lower = X.estimated.lower,
                         X.estimated.upper = X.estimated.upper, 
                         X.input=X.input) %>%
      mutate(covered = X.estimated.lower < X.input & X.estimated.upper>X.input) %>% 
      mutate(bias = X.estimated-X.input)
    
    
    
  }
  
  G = list()
  P = xdata$all.HH %>% mutate(tinfection = t.infection) 
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection","T_P",'T_U',"V_P" )){
    C = Chains[[var]]
    J=J+1
    G[[J]]=  summary_statistics_compare(var)
    
  }
  return(G)
  
}

MCMC_compare_summary_statistics_index_contact <- function(Chains, xdata){
  
  # load(filename)
  # Chains= rstan::extract(fit)
  summary_statistics_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    df.times= data.frame(X.estimated = X.estimated, 
                         X.estimated.lower = X.estimated.lower,
                         X.estimated.upper = X.estimated.upper, 
                         X.input=X.input,
                         status = xdata$all.HH$status) %>% 
      mutate(status = as.factor(status))%>%
      group_by(status)%>%
      mutate(covered = X.estimated.lower < X.input & X.estimated.upper>X.input)  
  }
  
  G = list()
  P = xdata$all.HH %>% mutate(tinfection = t.infection) 
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection","T_P",'T_U',"V_P" )){
    C = Chains[[var]]
    J=J+1
    G[[J]]=  summary_statistics_compare(var)
    
  }
  return(G)
  
}








# 
# 
# MCMC_compare_input_estimates_without_symptoms <- function(Chains, params,xdata, lims = c(25,10,15,5,15)){
#   
#   # load(filename)
#   # Chains= rstan::extract(fit)
#   plot_compare <- function(var){
#     xmin = 0
#     if(var == "tinfection")
#       xmin = -9
#     if(var == "TS")
#       xmin = -2
#     
#     X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
#     for(k in 1:xdata$n_households ){
#       for(j in 1:xdata$n_infected[k]){
#         
#         P1= P%>% filter(household.id == k)
#         P1=P1[j,]
#         X.estimated = c(X.estimated, median(C[,k,j]))
#         X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
#         X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
#         X.input = c(X.input,  P1[[var]])
#       }
#     }
#     print(X.estimated.lower)
#     print(X.estimated.upper)
#     df.times= data.frame(X.estimated = X.estimated, X.input=X.input)
#     g = ggplot() +geom_point(data=df.times,aes(y=X.estimated, x = X.input))+
#       geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower,yend = X.estimated.upper))+
#       theme_bw() +
#       geom_abline(slope=1, intercept = 0, linetype='dashed')+
#       theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
#             axis.text.y = element_text(size=16),
#             text=element_text(size=16))+
#       xlab(paste0(label, ' - Input'))+
#       ylab(paste0(label, ' - Posterior estimate'))+
#       xlim(c(xmin,lim))+ylim(c(xmin,lim))+coord_equal()
#     
#     return(g)
#   }
#   
#   G = list()
#   P = params %>% mutate(tinfection = t.infection) 
#   J=0
#   #for(var in c("tinfection","slope_1","slope_2","time_up" )){
#   for(var in c("tinfection","T_P",'T_U',"V_P" )){
#     C = Chains[[var]]
#     J=J+1
#     label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load')[J]
#     
#     lim = lims[J]
#     G[[J]]=  plot_compare(var)
#     
#   }
#   
#   grid.arrange(G[[1]],G[[2]],G[[3]],G[[4]],
#                ncol =2, nrow = 2)
#   
#   # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
#   # dev.off()
# }

MCMC_compare_input_estimates <- function(Chains, params,xdata, lims = c(25,10,15,5,15)){
  
  # load(filename)
  # Chains= rstan::extract(fit)
  plot_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    
    df.times= data.frame(X.estimated = X.estimated, X.input=X.input)
    g = ggplot() +geom_point(data=df.times,aes(y=X.estimated, x = X.input))+
      geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower,yend = X.estimated.upper))+
      theme_bw() +
      geom_abline(slope=1, intercept = 0, linetype='dashed')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+
      xlab(paste0(label, ' - Input'))+
      ylab(paste0(label, ' - Posterior estimate'))+
      xlim(c(xmin,lim))+ylim(c(xmin,lim))+coord_equal()
    
    return(g)
  }
  
  G = list()
  P = params %>% mutate(tinfection = t.infection) %>% mutate(TS = t.symptom)
  
  #mutate(slope_2 =slope.2) %>% 
  #mutate(time_up =time.up) %>% 
  #  mutate(tinfection = t.infection)
  
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection","T_P",'T_U',"V_P","TS" )){
    C = Chains[[var]]
    J=J+1
    label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load',"Time symptoms")[J]
    
    lim = lims[J]
    G[[J]]=  plot_compare(var)
    
  }
  
  grid.arrange(G[[1]],G[[2]],G[[3]],G[[4]],G[[5]],
               ncol =3, nrow = 2)
  
  # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
  # dev.off()
}
MCMC_compare_input_estimates_binary <- function(Chains, params,xdata, lims = c(25,5,15)){
  
  # load(filename)
  # Chains= rstan::extract(fit)
  plot_compare <- function(var){
    xmin = 0
    if(var == "tinfection")
      xmin = -9
    if(var == "TS")
      xmin = -2
    
    # 
    # if(var=="Delta_I" ){
    #   
    #   X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    #   for(k in 1:xdata$n_households ){
    #     for(j in 1:xdata$n_infected[k]){
    #       
    #       P1= P%>% filter(household.id == k)
    #       P1=P1[j,]
    # 
    #       X.input = c(X.input,  P1[[var]])
    #     }
    #   } 
    #   X.estimated =  median(C )
    #   X.estimated.lower = quantile025(C ) 
    #   X.estimated.upper =  quantile975(C) 
    # }else{
    X.estimated=X.input= X.estimated.lower = X.estimated.upper = c()
    for(k in 1:xdata$n_households ){
      for(j in 1:xdata$n_infected[k]){
        
        P1= P%>% filter(household.id == k)
        P1=P1[j,]
        X.estimated = c(X.estimated, median(C[,k,j]))
        X.estimated.lower = c(X.estimated.lower, quantile025(C[,k,j]))
        X.estimated.upper = c(X.estimated.upper, quantile975(C[,k,j]))
        X.input = c(X.input,  P1[[var]])
      }
    }
    
    #}
    
    
    df.times= data.frame(X.estimated = X.estimated, X.input=X.input)
    g = ggplot() +geom_point(data=df.times,aes(y=X.estimated, x = X.input))+
      geom_segment(data=df.times , aes(x = X.input, xend= X.input,y = X.estimated.lower,yend = X.estimated.upper))+
      theme_bw() +
      geom_abline(slope=1, intercept = 0, linetype='dashed')+
      theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
            axis.text.y = element_text(size=16),
            text=element_text(size=16))+
      xlab(paste0(label, ' - Input'))+
      ylab(paste0(label, ' - Posterior estimate'))+
      xlim(c(xmin,lim))+ylim(c(xmin,lim))+coord_equal()
    
    return(g)
  }
  
  G = list()
  P = params %>% mutate(tinfection = t.infection) %>% mutate(TS = t.symptom) %>% mutate(Delta_I = duration.detectability )  
  
  J=0
  #for(var in c("tinfection","slope_1","slope_2","time_up" )){
  for(var in c("tinfection", "TS","Delta_I" )){
    C = Chains[[var]]
    J=J+1
    print(var)
    label = c('Time infection', "Time symptoms","Duration Detectable")[J]
    
    lim = lims[J]
    G[[J]]=  plot_compare(var)
    
  }
  
  grid.arrange(G[[1]],G[[2]],G[[3]] ,
               ncol =3, nrow = 1)
  
  # dev.copy(pdf,paste0("../Results/Comparison_estimates_",filename,".pdf"), width=12, height =10)
  # dev.off()
}

MCMC_dataframe_comparison <- function(Chains, params,xdata){
  
  estimation_compare <- function(var){
    df = NULL
    C = Chains[[var]]
    
    if( var =="m" || var == "sigma_noise"){
      med =  median(C)
      lower = quantile025(C)
      upper = quantile975(C)
      me = mean(C)
      # input = P1[[var]] 
      
      if(var == "m"){
        input = xdata$m_input
      }
      if(var == "sigma_noise"){
        input = xdata$noise
      }
      
      df= rbind(df,data.frame(median = med,
                              lower =lower,
                              upper = upper, 
                              mean= me,
                              input=input,    
                              hh.id = 0,
                              indiv.id = 0,
                              var= var))
      
    }else{ 
      for(k in 1:xdata$n_households ){
        for(j in 1:xdata$n_infected[k]){
          P1= P%>% filter(household.id == k)
          P1=P1[j,]
          
          med =  median(C[,k,j])
          lower = quantile025(C[,k,j])
          upper = quantile975(C[,k,j])
          me = mean(C[,k,j])
          input = P1[[var]]
          
          df= rbind(df,data.frame(median = med,
                                  lower =lower,
                                  upper = upper,       
                                  mean= me,
                                  input=input,
                                  hh.id = k,
                                  indiv.id = j,
                                  var = var))
        }
      }
    }
    return(df)
  }
  G = NULL
  J=0 
  P = params %>% mutate(tinfection = t.infection)
  
  for(var in c("tinfection","T_P",'T_U',"V_P" ,"m","sigma_noise")){
    J=J+1
    #  label = c('Time infection', 'Time to peak', 'Time to clearance','Peak viral load')[J]
    G = rbind(G, estimation_compare(var) )
  }
  
  return(G)
  
}

MCMC_bias_estimates <- function(Chains, params,xdata){
  
  ## Boxplot of the bias
  Bias = NULL
  for(foi in c(0,1)){
    J=0
    for(var in c("tinfection","slope_1","slope_2","time_up" )){
      C = Chains[[var]]
      J=J+1
      label = c('Time infection', 'Slope 1', 'Slope 2', 'Time to peak')[J]
      
      X.estimated = X.input = c()
      for(k in 1:xdata$n_households ){
        for(j in 1:xdata$n_infected[k]){
          
          P1= P%>% filter(household.id == k)
          P1=P1[j,]
          X.estimated = c(X.estimated, median(C[,k,j]))
          X.input = c(X.input,  P1[[var]])
        }
      }
      Bias =  rbind(Bias,data.frame( var = var, bias = X.estimated-X.input, norm.bias =X.estimated/X.input-1 ))
    }
  }
  
  Bias = Bias %>%mutate(var = as.factor(var))  %>% filter(!is.infinite(norm.bias))
  
  g = ggplot(Bias, aes(x=var, y=norm.bias*100)) + 
    geom_boxplot(outlier.colour=foi, outlier.shape=NA,
                 outlier.size=4)+ylim(c(-50,70))+
    theme_bw() +
    geom_hline(aes( yintercept = 0), linetype='dashed')+
    theme(axis.text.x = element_text(size=16,angle = 0, vjust=0),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+xlab('')+ylab('Bias (%)')
  
  #  print(g)
  #dev.copy(pdf,paste0("../Results/Bias_estimates_",filename,".pdf"), width=12, height =10)
  # dev.off()
  
} 