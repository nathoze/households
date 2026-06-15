# transform_to_antigenic_test <-function(xdata,specificity_ag  , sensitivity_ag , LOD_ag  ){
#   
#   
#   if(LOD_ag<6){
#     print("In principle LOD_ag should be >6")
#   }
#   xdata$Binary_obs  = matrix(as.numeric(ifelse(xdata$LVLObs>LOD_ag,rbernoulli(n=1, p=sensitivity_ag), rbernoulli(n=1, p=1-specificity_ag))), 
#                              nrow = nrow(xdata$LVLObs),
#                              ncol = ncol(xdata$LVLObs))
#   
#   xdata$specificity_ag  = specificity_ag
#   xdata$sensitivity_ag = sensitivity_ag
#   xdata$LOD_ag = LOD_ag
#   
#   dt=0.01
#   L1=rep(0,nrow(xdata$all.HH))
#   for(i in 1:nrow(xdata$all.HH)){
#     discrete.values = get_viral_load(params =   xdata$all.HH[i,] , delay= 2, Tmax = 30, dt = dt)   
#     L1[i] = length(which(log10(discrete.values$VI+discrete.values$VNI)>LOD_ag))*dt
#   }  
#   xdata$all.HH$duration.detectability = L1
#   
#   
#   # xdata$Binary_obs  = ifelse(xdata$LVLObs>LOD_ag,rbernoulli(n=1, p=sensitivity_ag),0)
#   return(xdata)
#   
# }


transform_to_antigenic_test <-function(xdata,prob_symptomatic, specificity_ag =1 , sensitivity_ag_asymptomatic, sensitivity_ag_symptomatic, LOD_ag  ){
  
  if(LOD_ag<6){
    print("In principle LOD_ag should be >6")
  }
  
  observation.table = xdata$LVLObs*0
  symptomatic =xdata$symptomatic
  
  
  xdata$n_infected_binary = rep(0, xdata$n_households)
  xdata$infected_array = array(0, dim=c(xdata$n_households,xdata$max_infected))
    infected = rep(0, xdata$n_infected_total)
# 
#   I=0
#   for(k in 1:xdata$n_households){
#     for(i in 1:xdata$n_infected[k]){
#       I=I+1
#       infected[I] = runif(n=1)<prob_symptomatic
#      }
#   }
#   
  
  #  sensitivity_ag_asymptomatic=0.5
  # sensitivity_ag_symptomatic = 0.7
  I=0
  x=0
  for(k in 1:xdata$n_households){
    for(i in 1:xdata$n_infected[k]){
      
      I=I+1
      U=xdata$start[k,i]:xdata$end[k,i]
      post.symptoms = xdata$time[k,U] > xdata$t_symptoms[I] 
      lvlobs = xdata$LVLObs[k,U] 
      observation.vec = rep(0, length( observation.table[k,U]))
      
      for( l in 1:length(observation.vec)) {
        if(lvlobs[l]>LOD_ag){
          x=x+1
          if(post.symptoms[l] == FALSE || symptomatic[I]==FALSE ){
            observation.vec[l] = runif(n=1)< sensitivity_ag_asymptomatic
          }
          if(post.symptoms[l] == TRUE & symptomatic[I]==TRUE ){
            observation.vec[l] = runif(n=1)< sensitivity_ag_symptomatic
          }
        }
      }
      observation.table[k,U] = observation.vec
      if(sum(observation.vec)>0){
        xdata$n_infected_binary[k] = xdata$n_infected_binary[k] +1
        xdata$infected_array [k,i]  =1
        infected[I] = 1 
      }
    }
  }
  
  xdata$symptomatic = symptomatic
  xdata$infected = infected
  xdata$Binary_obs = observation.table
  
  xdata$n_non_infected_binary = xdata$n_non_infected +  xdata$n_infected -  xdata$n_infected_binary 
  xdata$n_infected_binary_total = sum(xdata$n_infected_binary) 
  
  xdata$sensitivity_ag_asymptomatic = sensitivity_ag_asymptomatic
  xdata$sensitivity_ag_symptomatic = sensitivity_ag_symptomatic
  xdata$specificity_ag = specificity_ag
  # xdata$Binary_obs  = matrix(as.numeric(ifelse(xdata$LVLObs>LOD_ag,rbernoulli(n=1, p=sensitivity_ag), rbernoulli(n=1, p=1-specificity_ag))), 
  #                            nrow = nrow(xdata$LVLObs),
  #                            ncol = ncol(xdata$LVLObs))
  # 
  # xdata$specificity_ag  = specificity_ag
  # xdata$sensitivity_ag = sensitivity_ag
  # xdata$LOD_ag = LOD_ag
  
  dt=0.01
  L1=rep(0,nrow(xdata$all.HH))
  for(i in 1:nrow(xdata$all.HH)){
    discrete.values = get_viral_load(params =   xdata$all.HH[i,] , delay= 2, Tmax = 30, dt = dt)   
    L1[i] = length(which(log10(discrete.values$VI+discrete.values$VNI)>LOD_ag))*dt
  }  
  xdata$all.HH$duration.detectability = L1
  
  
  # xdata$Binary_obs  = ifelse(xdata$LVLObs>LOD_ag,rbernoulli(n=1, p=sensitivity_ag),0)
  return(xdata)
  
}
