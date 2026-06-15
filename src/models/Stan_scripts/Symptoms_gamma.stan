
functions{
  
  //real integral_gamma_model(real t_infection_infector, real t_infection_infected, real t_symptoms_infector, real gamma_shape, real gamma_inverse_scale) {
  real integral_gamma_model(real t_infection_infector, real t_infection_infected, real gamma_shape, real gamma_inverse_scale) {
    real VL_integral;
    VL_integral = 0;
    real t_lag = t_infection_infector;//+t_symptoms_infector;
    real y; 
      y =   t_infection_infected-  t_infection_infector ;

    if(y<=0){
      VL_integral +=  0.001;
    }
    if(y>0){
      VL_integral +=  0.001;
      VL_integral += gamma_cdf(y |  gamma_shape,  gamma_inverse_scale)/gamma_cdf(1000000 |  gamma_shape,  gamma_inverse_scale);        
    }    
    return VL_integral;
   } 
  vector log_gamma_model(real t_infection_infector, real t_infection_infected, real gamma_shape, real gamma_inverse_scale) {
    //vector[n] vi; 
    vector[1] vi; 
    real t_lag = t_infection_infector;//+t_symptoms_infector;
    real y; 
    y = t_infection_infected - t_infection_infector;
  
    vi[1] = 0; 
    if(y>0){  
      vi[1] += gamma_lpdf( y |  gamma_shape,  gamma_inverse_scale);
    }
    return vi;
 } 

  real integral_gamma_model_max_time(real gamma_shape, real gamma_inverse_scale) {
    real VL_integral;

    real y = 1000000;  
    VL_integral  =  gamma_cdf(y |  gamma_shape,  gamma_inverse_scale);     
    return VL_integral;
   } 
 
  // computes fmax on all element of a vector 
  // (to select the highest value between 0 and log10(LVL))
  row_vector new_fmax(vector v, int npoints){
    row_vector[npoints] new_v;
    for(i in 1:npoints)
    new_v[i] = fmax(0, log10(v[i]));
    return new_v;
  }
  
  row_vector left_censored(vector v, int npoints, real min_value){
    row_vector[npoints] new_v;
    for(i in 1:npoints){
      new_v[i] =  v[i];
      if(new_v[i] <min_value){
        new_v[i] = min_value; 
      }
    }
    return new_v ;
  }  
}

data{
  int<lower = 1> n_households;       // Nb households  
  array[n_households] int<lower = 1> size_households;    // Nb subjects for each households
  array[n_households] int<lower = 0> n_infected;    // Nb infected subjects for each households
  array[n_households] int<lower = 0> n_non_infected;    // Nb non infected subjects for each households  
  int<lower = 1> n_subjects ; // Total number of persons summed over all households
  int<lower = 1> NT_max;       // Nb de points totaux maximal (max(NT))
  int max_infected ; // maximal number of infected across households  
  array[n_households, max_infected] int<lower = 0> start;// Indices de debut pour chaque indivs 
  array[n_households, max_infected] int<lower = 0> end;// Indices de fin pour chaque indivs 
  array[n_households,NT_max] real<lower = 0> time;       // all time points of observations (for each individual) 
  int n_time ; // number of time points 
  real Delta_S_star; 
  real Delta_S_eta; 
  real sigma_prior[2];  // Prior observation noise Cauchy scale
  int n_infected_total ; // sum of all infected over all the households
  int<lower=0, upper = n_households> row_infected[n_infected_total];
  int<lower=0, upper = max_infected> col_infected[n_infected_total];
  real<lower=0> max_time;
  real t_symptoms[n_infected_total];
  real symptomatic[n_infected_total];  // 0 or 1 if asymptomatic or symptomatic 
}

parameters{
  // Population-level model parameters
  // These are the parameters for which you specify prior distributions and initial estimates,

  real logm;  // the transmission parameter
  real alpha; // exponent for the household size   
  real<lower= 0> gamma_infection_shape;
  real<lower= 0> incubation;

 }

transformed parameters{
  real gamma_infection_inverse_scale = 0.3;
  // Individual-level model parameters,
  matrix[n_households, max_infected] TS ; // = T_S-T_I
  matrix[n_households, max_infected] tinfection;
  vector[n_infected_total] TI;  // Contains the same informaion as tinfection but with a different 
  
 // matrix[n_households, NT_max] LVL;               // Predicted viral loads
  real<lower=0> m; 
  m = 10^(logm);
  vector[1] lvi;  
  vector[1] tinf ; 
  real LLinfection1;  // All the observations before the time of infection
  real LLinfection2 ;  // At the time of infection
  real LLinfection3 ;  // the non infected
  real LLsymptoms ; // In the likelihood, the contribution of the delay between infection and symptoms
  real loglik;

  real P;
  loglik = 0;

  LLinfection1 = 0;
  LLinfection2 = 0;
  LLinfection3 = 0;
  LLsymptoms = 0; 

  for(k in 1:n_households){
    for(j in 1:max_infected){
      tinfection[k,j] = 0;
      TS[k,j] = 1000;
    }
  } 

  for(k in 1:n_infected_total){
    if(symptomatic[k] == 1){
      TS[row_infected[k],col_infected[k]] = incubation; 
      tinfection[row_infected[k],col_infected[k]] = t_symptoms[k] - TS[row_infected[k],col_infected[k]]; //print(k);
      TI[k] =  tinfection[row_infected[k],col_infected[k]] ; 
    }
  } 

  for(k in 1:n_infected_total){     
    if(symptomatic[k] == 1){
      LLsymptoms +=  normal_lpdf(log( TS[row_infected[k],col_infected[k]])|  log(Delta_S_star), Delta_S_eta); // with input value of the incubation period
     } 
  }
 
     real Int ; 
    for(k in 1:n_households){ 
      for(j in 2:n_infected[k]){  // we don't consider here the first individual who is the index case 
        tinf[1] =  tinfection[k,j];
        for(i in 1:n_infected[k]){ // sum over the possible infectors of individual j
          if(i != j){
           //  Int = integral_gamma_model(tinfection[k,i], tinfection[k,j],TS[row_infected[k],i], gamma_infection_shape, gamma_infection_inverse_scale );//tinf[1],  T_P[k,i], T_U[k,i],V_P[k,i],  tinfection[k,i]); // the cumulative hazard for individual j up to time tinf[1] 
             Int = integral_gamma_model(tinfection[k,i], tinfection[k,j],gamma_infection_shape, gamma_infection_inverse_scale );//tinf[1],  T_P[k,i], T_U[k,i],V_P[k,i],  tinfection[k,i]); // the cumulative hazard for individual j up to time tinf[1] 
            LLinfection1+= -m/((size_households[k]-1)^alpha)*Int;
          }  
       }

       P = 0;
       P = 0.00000001;
       // I could use log_sum_exp 
        for(i in 1:n_infected[k]){ // sum over the possible infectors of individual j
          if(i != j){
            lvi= log_gamma_model(tinfection[k,i], tinfection[k,j],  gamma_infection_shape, gamma_infection_inverse_scale );
            P = P  +  exp(lvi[1]);
          }
        } 
        LLinfection2 += log(m/(size_households[k]-1)^alpha) + log(P);
      }
    }

      for(k in 1:n_households){
          LLinfection3+= -m/((size_households[k]-1)^alpha)*n_non_infected[k]*(n_infected[k]);
       }

    loglik = LLsymptoms+LLinfection1+LLinfection2+LLinfection3;
  }
  
  model{
    
    gamma_infection_shape ~ normal(10, 10) T[0,]; 
    incubation ~ normal(3, 3) T[0,]; 
    alpha ~ normal(0.5,1) ;//T[0,]; 
    logm ~ normal(-3,3);    
  //    target+= LLobs;
      target+=LLinfection1;
      target+=LLinfection2;
      target+=LLinfection3;
      target+=LLsymptoms;
  }

  /*
  generated quantities{
  // matrix[n_households, NT_max] LVLpred;
  //  matrix[n_households, NT_max] LVLpred_small_noise;
  
  //matrix[n_households, n_inf_times_max] LVL_triangle;
  
  //for(k in 1:n_households){
  //  for(j in 1:n_infected[k]){ 
  
  //  LVLpred[k,start[k,j]:end[k,j]] = left_censored( to_vector(normal_rng( LVL[k,start[k,j]:end[k,j]],sigma_noise )), n_time, 6.0);
  //  LVLpred_small_noise[k,start[k,j]:end[k,j]] = left_censored( to_vector(normal_rng( LVL[k,start[k,j]:end[k,j]], 0.000001 )), n_time, 6.0); 
  //    LVL_triangle[k,start_inf_times[k,j]:end_inf_times[k,j]] = left_censored(   to_vector(logVNI_all[k,start_inf_times[k,j]:end_inf_times[k,j]] ), n_discrete_times, 6.0);
  //   }
  // }
  
  }*/
  
