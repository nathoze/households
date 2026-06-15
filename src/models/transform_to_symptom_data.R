

transform_to_symptom_data <-function(xdata,prob_symptomatic){
 
  symptomatic =rep(0, xdata$n_infected_total)
  symptomatic_array = array(0, dim=c(xdata$n_households,xdata$max_infected))
  
  t_symptoms=c()
  I=0
  for(k in 1:xdata$n_households){
    for(i in 1:xdata$n_infected[k]){
      I=I+1
      symptomatic[I] = runif(n=1)<prob_symptomatic
      if(i == 1)
        symptomatic[I] = 1
      symptomatic_array[k,i] =  symptomatic[I]
      if( symptomatic[I]  ==1){
        t_symptoms=c(t_symptoms, xdata$t_symptoms[I])
      }
    }
  }
  xdata$symptomatic = symptomatic
  xdata$symptomatic_array = symptomatic_array
  
  xdata$n_symptomatic = rowSums(xdata$symptomatic_array)
  xdata$n_non_symptomatic = xdata$size_households-n_symptomatic
  
  max_households = max(xdata$size_households)
  n_households = xdata$n_households
  index_symptomatic = t(array(data= seq(1,n_households*max_households), dim =c( max_households, n_households)))
  
  is_symptomatic = array(data= 0, dim =c(  n_households,max_households))
  for(k in 1:xdata$n_households){
    for(i in 1:sum(symptomatic_array[k,])){
      is_symptomatic[k,i] = 1
    }
  }
  
  
   
  xdata_2 = list(n_households = n_households,
                 size_households = xdata$size_households, 
                 n_symptomatic = xdata$n_symptomatic, 
                 n_symptomatic = xdata$n_symptomatic, 
                 n_symptomatic_total = sum(xdata$n_symptomatic),
                 n_non_symptomatic=xdata$n_non_symptomatic,
                 symptomatic_array = xdata$symptomatic_array,
                 symptomatic = xdata$symptomatic,
                 max_symptomatic = max(xdata$n_symptomatic),
                 prob_symptomatic = prob_symptomatic,
                 t_symptoms = t_symptoms,
                 n_subjects = xdata$n_subjects,
                 max_households =max_households,
                 index_symptomatic =  index_symptomatic)
  
  
  return(xdata)
  
}




 