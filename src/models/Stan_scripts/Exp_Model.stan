
functions{

  vector triangle_model_obs(vector t, int n, real T_P, real T_U, real V_P, real tlag) {
    vector[n] vi; 
    real t_detectable = 4.0/V_P*T_P;
    for(i in 1:n){
      if(t[i] < tlag+t_detectable){
        vi[i] = 2.0; 
      }
      else if(t[i] <tlag+T_P+t_detectable){    
        vi[i]  = 2.0 + (V_P)/T_P*(t[i]-tlag-t_detectable); 
      }
      else if(t[i] >=tlag+T_P+t_detectable){
        vi[i] = V_P + 2.0 + (2.0-(V_P+2.0))/T_U*(t[i] - tlag - T_P-t_detectable ); 
      }
    }
    return vi;
  } 

  
  vector triangle_model(vector t, int n, real T_P, real T_U, real V_P, real tinf) {
    vector[n] vi; 
    real t_detectable = 4.0/V_P*T_P;

    for(i in 1:n){
        vi[i] = -2.0; 
      }

    for(i in 1:n){
      if(t[i] <= tinf){
          vi[i] = -2.0;
      }
      if(t[i] > tinf && t[i] < tinf+t_detectable ){
        vi[i] = -2.0 + 4.0*(t[i]-tinf)/t_detectable ; 
      }
      if(t[i] >= tinf+t_detectable && t[i] <tinf+T_P+t_detectable){
          vi[i]  = 2.0 + (V_P)/T_P*(t[i]-tinf-t_detectable); 
      }
      if(t[i] >=tinf+T_P+t_detectable){
        vi[i] = V_P + 2.0 + (2.0-(V_P+2.0))/T_U*(t[i] - tinf - T_P-t_detectable ); 
      }
    }
    return vi;
  } 

real integral_triangle_model(real t,  real T_P, real T_U, real V_P, real tinf) {

   real VL_integral;
    VL_integral = 0.01;
    real ln10 = 2.30258;
    real h=0.5; 
    real t_detectable = 4.0/V_P*T_P;
  
  real A = 10^(2*h)*(1/V_P/(h*ln10))*(10^(h*V_P)-1);

  real t0 = t-tinf;
  real u;

  if(t0>= 0){    
    if(t0 < T_P+t_detectable) {
      VL_integral =VL_integral+10^(-2*h)*(T_P/V_P/(h*ln10))*(10^(h*V_P/T_P*(t0))-1);
      
    }
    if(t0 >= T_P+t_detectable) {
      VL_integral =VL_integral+10^(2*h)*(T_P/V_P/(h*ln10))*(10^(h*V_P)-1);
      
      if( t0 < T_P+T_U+t_detectable) {
         u= T_P+T_U+t_detectable - t0 ;
        VL_integral =VL_integral  - 10^(2*h)*(T_U/V_P/(h*ln10))*(10^(h*V_P/T_U*(u))-1) + 10^(2*h)*(T_U/V_P/(h*ln10))*(10^(h*V_P)-1);
      }
      if(t0 >= T_P+T_U+t_detectable) {
        VL_integral =VL_integral+10^(2*h)*(T_U/V_P/(h*ln10))*(10^(h*V_P)-1);
      }
    }
  }    
    return VL_integral;
} 
 
 real integral_triangle_model_full(real T_P, real T_U, real V_P) {
   real VL_integral;
    VL_integral = 0.01;
    real ln10 = 2.30258;
    real h=0.5; 
  
  VL_integral = 10^(2*h)*(T_U+T_P)/(V_P*h*ln10)*(10^(h*V_P)-1);
 
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
  matrix[n_households, NT_max]  LVLObs;              // All viral load 
  matrix[n_households, NT_max]  censored;    // whether an observation is left censored   
  int n_time ; // number of time points 
  real h; // one of the transmission parameters
  real K; //#100^h
  //int with_foi; // Whether we estimate the individual parameters without accounting for the infections (with_foi=0)  
  real T_P_star; 
  real T_U_star; 
  real V_P_star; 
  real Delta_S_star; 
  real Delta_S_eta; 
  real sigma_prior[2];  // Prior observation noise Cauchy scale
  
  int n_infected_total ; // sum of all infected over all the households
  int<lower=0, upper = n_households> row_infected[n_infected_total];
  int<lower=0, upper = max_infected> col_infected[n_infected_total];
  real first_detectable_time[n_infected_total];
  real<lower=0> max_time;
  real t_symptoms[n_infected_total];
  real symptomatic[n_infected_total];  // 0 or 1 if asymptomatic or symptomatic 
}

parameters{
  // Population-level model parameters
  // These are the parameters for which you specify prior distributions and initial estimates,
  vector[n_infected_total] T_P_eta;
  vector[n_infected_total] T_U_eta;
  vector[n_infected_total] V_P_eta;
  vector<lower=0>[n_infected_total]  TI_tmp_positive; 
  real log_T_P_mu;
  real log_T_U_mu;
  real log_V_P_mu;  
  real<lower=0> log_T_P_sd;
  real<lower=0> log_T_U_sd;
  real<lower=0> log_V_P_sd;
   
  real<lower=0, upper=2> sigma_noise;
  real logm;  // the transmission parameter
  real alpha; // exponent for the household size 
  // Parameters of the incubation period
  real log_Delta_S_mu   ; 
  real<lower=0> Delta_S_sd ; 
 }

transformed parameters{
  // Individual-level model parameters
  matrix[n_households, max_infected] T_P ; //=  exp(T_P_eta )*T_P_star;
  matrix[n_households, max_infected] T_U ;//=  exp(T_U_eta )*T_U_star;
  matrix[n_households, max_infected] V_P ;// =  exp(V_P_eta )*V_P_star;
  matrix[n_households, max_infected] TS ;// the "reconstructed" symptom time 

  matrix[n_households, max_infected] tinfection;   
  vector[n_infected_total] TI;  // Contains the same informaion as tinfection but with a different 

  matrix[n_households, NT_max] LVL;               // Predicted viral loads
  real<lower=0> m; 
  m = 10^(logm);
  //real I;
  vector[1] lvi;  
  vector[1] tinf ; 
  real LLobs ; // the observation model
  real LLinfection1;  // All the observations before the time of infection
  real LLinfection2 ;  // At the time of infection
  real LLinfection3 ;  // the non infected
  real LLsymptoms ; // In the likelihood, the contribution of the delay between infection and symptoms

  real P;
  real loglik;
  LLobs = 0;
  LLinfection1 = 0;
  LLinfection2 = 0;
  LLinfection3 = 0;
  LLsymptoms = 0;
 
  for(k in 1:n_households){
    for(j in 1:max_infected){
      T_P[k,j] = 0;
      T_U[k,j] = 0;
      V_P[k,j] = 0;
      tinfection[k,j] = 0;
      TS[k,j] = 1000;
    }
  } 

  for(k in 1:n_infected_total){
     T_P[row_infected[k],col_infected[k]] = exp(log_T_P_mu+log_T_P_sd*T_P_eta[k])*T_P_star;
      T_U[row_infected[k],col_infected[k]] = exp(log_T_U_mu+log_T_U_sd*T_U_eta[k])*T_U_star;
      V_P[row_infected[k],col_infected[k]] = exp(log_V_P_mu+log_V_P_sd*V_P_eta[k])*V_P_star;
      tinfection[row_infected[k],col_infected[k]] = first_detectable_time[k] - TI_tmp_positive[k]; 
      if(symptomatic[k] == 1){
        TS[row_infected[k],col_infected[k]] =  t_symptoms[k] - tinfection[row_infected[k],col_infected[k]] ;
      }
      TI[k] =  tinfection[row_infected[k],col_infected[k]] ; // another indexation of the time of symptom onset
   } 
 
  for(k in 1:n_households){
  for(j in 1:NT_max){
      LVL[k,j] = 2;
    }
  } 

  for(k in 1:n_households){
  for(j in 1:n_infected[k]){     
      LVL[k,start[k,j]:end[k,j]] = to_row_vector( triangle_model_obs( to_vector(time[k,start[k,j]:end[k,j]]), n_time,  T_P[k,j], T_U[k,j], V_P[k,j], tinfection[k,j])  ) ; // VI is negligeble, so the observable log viral load is the non infectious viruses     
    }
  } 

  for(k in 1:n_households){        
    for(i in 1: end[k,n_infected[k]]){
      if(censored[k,i] == 0){   
        LLobs += normal_lpdf(LVLObs[k, i] |LVL[k,i], sigma_noise ) ;
      }      
      if(censored[k,i] == 1){

        LLobs += normal_lcdf(LVLObs[k, i] |LVL[k,i], sigma_noise ); // left censored

      }    
    }
  }

  for(k in 1:n_infected_total){     
    if(symptomatic[k] == 1){
        LLsymptoms +=  normal_lpdf(log(TS[row_infected[k],col_infected[k]])|  log_Delta_S_mu, Delta_S_sd);  
      } 
    }
  
     real Int ;
     real Int3 ;  
 
    for(k in 1:n_households){
      for(j in 2:n_infected[k]){  // we don't consider here the first individual who is the index case 
        tinf[1] =   tinfection[k,j];
        for(i in 1:n_infected[k]){ // sum over the possible infectors of individual i
          if(i != j){
            Int= integral_triangle_model( tinf[1],  T_P[k,i], T_U[k,i], V_P[k,i], tinfection[k,i]); // the cumulative hazard for individual j up to time tinf[1] 
            LLinfection1+= -m/((size_households[k]-1)^alpha)*Int;
          }  
       }

       P = 0;
        for(i in 1:n_infected[k]){ // sum over the possible infectors of individual i
          if(i != j){
            lvi = triangle_model( tinf, 1, T_P[k,i], T_U[k,i], V_P[k,i], tinfection[k,i]); //   log infectious viral load for individual j
              P = P  + 10^(h*lvi[1]);
          }
        } 
        LLinfection2 += log(m/((size_households[k]-1)^alpha)) + log(P);
      }
    }

      for(k in 1:n_households){
        Int3=0;
        for(j in 1:n_infected[k]){ // sum over the possible infectors of individual i
            Int3+= integral_triangle_model_full( T_P[k,j], T_U[k,j], V_P[k,j]);
          }  
          LLinfection3+= -m/((size_households[k]-1)^alpha)*n_non_infected[k]*(Int3);
       }

     loglik = LLobs+LLsymptoms+LLinfection1+LLinfection2+LLinfection3;
  }
  
  model{    
    log_T_P_mu ~ normal(0, 0.25); 
    log_T_U_mu ~ normal(0, 0.25); 
    log_V_P_mu ~ normal(0, 0.25);  
    log_T_P_sd ~ normal(0, 0.25) T[0,]; 
    log_T_U_sd ~ normal(0, 0.25) T[0,]; 
    log_V_P_sd ~ normal(0, 0.25) T[0,]; 
    alpha ~ normal(0.5,1) ;//T[0,]; 
    log_Delta_S_mu  ~ normal(0, 0.25); 
    Delta_S_sd ~ normal(0, 0.25) T[0,]; 

    for(k in 1:n_infected_total){
      T_P_eta[k] ~ std_normal();
      T_U_eta[k] ~ std_normal();
      V_P_eta[k] ~ std_normal(); 
      TI_tmp_positive[k] ~ normal(4,4)  T[0,]; 
   } 
    
    sigma_noise ~ normal(sigma_prior[1],sigma_prior[2]) T[0,];
   // logm ~ normal(-3,3);    
    logm ~ normal(-5,3);    

    target+= loglik;
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
  
