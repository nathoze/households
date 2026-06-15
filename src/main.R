#rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(odin)
library(grid)
library(gridExtra)
library(png)
library(ggpubr)
library(cowplot)
library(posterior)
library(rstan)
library(labeling)

source("models/functions_household_transmission.R")
source("models/generate_parameters_triangle.R")
source("models/generate_observed_VI_symptoms.R")
source("models/reconstruction_chain_transmission.R")
source("models/run_MCMC.R")
source("visualization/MCMC_predictions.R")

minimal.logRNA = 2.39794
maximal.CT = 40.93733

colors.m.values = c("#21c271", "#a5bfcc","#3e9ecf","#d5ede1")
colors.m.values.dark = c("#116211", "#8193ac","#1e4eaf","#a5ada1")

detection.limit= 2 
infection.threshold = detection.limit
initial.logV.value = -2  
diff.vl = detection.limit - initial.logV.value

#epid.model = "exp" 
#epid.model = "constant" 
model.virus = "triangle"#model.virus ="ode"


Tmax = 70
discrete.dt  = 0.1
discrete.times = seq(0,Tmax,by=discrete.dt)
n.discrete.times = length(discrete.times)

dt = 0.01
total.time <- seq(0, Tmax, by = dt)
indices.discrete = match(discrete.times,total.time)

complete_vector <-function(x, N){
  y = rep(0,N)
  y[1:length(x)] =x 
  return(y)
} 
# From Kissler et al.
# concentration of viral RNA is in copies/ml
logRNA_2_CT<-function(logRNA){
  return((logRNA-log10(250))*(-3.60971)+maximal.CT)
}
CT_2_logRNA<-function(CT){  
  return(log10(250) +(CT - maximal.CT)/(-3.60971))
}
logRNA_2_CT_observable<-function(logRNA){
  CT = (logRNA-log10(250))*(-3.60971)+maximal.CT
  CT[which(CT>maximal.CT)] = maximal.CT
  return(CT)
}
CT_2_logRNA_observable<-function(CT){  
  CT[which(CT>maximal.CT)] = maximal.CT
  logRNA = CT_2_logRNA(CT)
  return(logRNA)
}
logRNA_observable<-function(V){  
  logRNA = log10(V)  
  return(logRNA)
}
logit <- function(p){
  return(log(p/(1-p)))
}
inverse.logit <-function(x){
  return(1/(1+exp(-x)))
}
quantile025 <- function(X){
  return(as.numeric(quantile(X, probs=0.025 )))
}
quantile975 <- function(X){
  return(as.numeric(quantile(X, probs=0.975 )))
}
quantile05 <- function(X){
  return(as.numeric(quantile(X, probs=0.05 )))
}
quantile95 <- function(X){
  return(as.numeric(quantile(X, probs=0.95 )))
}

quantile25 <- function(X){
  return(as.numeric(quantile(X, probs=0.25 )))
}
quantile75 <- function(X){
  return(as.numeric(quantile(X, probs=0.75 )))
}

observed_titers <- function(params, sigma, obs.times){
  df = get_viral_load(params) %>%
    filter(t %in% obs.times) %>% #find the closest time
    mutate(observed.logRNA = rnorm(n = length(obs.times), mean = logRNA, sd = sigma))%>%
    mutate(observed.logRNA = ifelse(observed.logRNA > minimal.logRNA, observed.logRNA, minimal.logRNA))# Adds a left censoring 
  return(df)
}
# observed_titers_three_days <- function(all.titers){
#   df =all.titers %>% 
#     filter(t %in% seq(0,Tmax,3)) %>% #find the closest time
#     mutate(observed.logRNA = rnorm(n= n(),mean = logRNA, sd = sigma)) %>%
#     mutate(observed.logRNA = ifelse(observed.logRNA > minimal.logRNA, observed.logRNA, minimal.logRNA))# Adds a left censoring 
#   return(df)
# }
#  s
## T_P is the time interval between LOD and maximal value (V_P + LOD)
## T_U is the time interval between maximal value (V_P + LOD) and LOD
## The slope s1 = V_P/T_P; the slope s2  = V_P / T_U
## The time to reach the LOD from the initial value -2 is t.detectable, given by -2 + s1*t.detectable = 2, or t.detectable= 4/V_P*T_P =diff.vl/V_P*T_P

model_gen_triangle <-function(T_P, T_U, V_P, t, tinfection = 0, normalize.time = TRUE){
  if(normalize.time){
    t0=t-t[1]
  }else
    t0=t
  n = length(t0)
  logv = rep(-5, n)
  t.detectable = diff.vl/V_P *T_P
  for(i in 1:n){
    if(t0[i]>=tinfection){
      if(t0[i]<T_P+tinfection+t.detectable){
        #     logv[i]= initial.logV.value + (V_P+ diff.vl)/T_P*(t0[i]-tinfection)
        logv[i]= initial.logV.value + V_P/T_P*(t0[i]-tinfection)
      }
      if(t0[i]>=T_P+tinfection+t.detectable){
        logv[i] = V_P+ detection.limit +(-V_P)/T_U*(t0[i] -tinfection- T_P -t.detectable)
      }
    }
  }
  df = data.frame(t= t0, V  = 10^(logv ),  logV  = logv  )
  return(df)
}


integral_triangle_model <- function(t, params, tlag){
  T_P = params$T_P
  T_U = params$T_U
  V_P = params$V_P
  # m = params$m
  h = params$h
  VL_integral = 0
  ln10 = log(10)
  #K =   10^(2*h)
  K =   10^(detection.limit*h)
  t.detectable = diff.vl/V_P *T_P
  tlag = tlag  + t.detectable
  
  if(t< tlag){
    VL_integral = VL_integral+ t*K
  }
  if(t>= tlag){
    VL_integral = VL_integral+ tlag*K
    
    if(t < tlag+T_P) {
      VL_integral = VL_integral+ K*(T_P/V_P/(h*ln10))*(10^(h*V_P/T_P*(t-tlag))-1)
    }
    if(t >= tlag+T_P) {
      VL_integral =VL_integral+ K*(T_P/V_P/(h*ln10))*(10^(h*V_P)-1)
      if(t < tlag+T_P+T_U) {
        # VL_integral=VL_integral+  (K*(T_U/V_P/(h*ln10))*(10^(h*V_P)-1) -K*(T_U/V_P/(h*ln10))*(10^(h*V_P/T_U*( T_U - t+tlag+T_P))-1))
        VL_integral=VL_integral+ K*(T_U/V_P/(h*ln10))*(10^(h*V_P)- 10^(h*V_P/T_U*( T_U - t+tlag+T_P)))
      }  
      if(t >= tlag+T_P+T_U) {
        VL_integral=VL_integral+K*(T_U/V_P/(h*ln10))*(10^(h*V_P)-1)
      }
      
    }
  }    
  return(VL_integral)
  
} 

get_viral_load <-function(params, delay,Tmax, dt){
  
  total.time <- seq(0, Tmax, by = dt)
  n.times = length(total.time)
  
  viral_load=NA
  if(delay <Tmax){
    simulation.time = seq(delay, Tmax, by = dt)
    
    n0 = n.times-length(simulation.time)
    if(model.virus == "ode"){
      x <- model_generator$new(T_ini = params$T_ini, 
                               V_ini = params$V_ini, 
                               beta = params$beta, 
                               kappa= params$kappa,
                               delta = params$delta,
                               phi = params$phi,
                               theta = params$theta, 
                               mu = params$mu,
                               p = params$p,
                               c = params$c,
                               dF = params$dF)
      y <- x$run(t = simulation.time)
      y <- as.data.frame(y)
    }
    
    if(model.virus == "triangle"){
      y <- model_gen_triangle(T_P = params$T_P,
                              T_U = params$T_U,
                              V_P = params$V_P,
                              t = simulation.time)
    }
    
    viral_load = y %>%
      slice(rep(1, n0)) %>%
      bind_rows(.,y)%>%
      # mutate(t = seq(0,by=dt,length.out = n()))%>%
      mutate(t = seq(0,by=dt,length.out = nrow(.)))%>%
      filter(t <= Tmax) %>%
      mutate(logRNA = logRNA_observable(V ))
  }
  
  return(viral_load) 
}

get_viral_load_symptom_onset <-function(params, delay, Tmax, dt, observation.times = NULL){
  
  if(delay <Tmax){
    if(is.null(observation.times) ){
      observation.times = seq(delay, Tmax, by = dt)
    }
    
    viral_load=NA
    
    if(model.virus == "ode"){
      x <- model_generator$new(T_ini = params$T_ini, 
                               V_ini = params$V_ini, 
                               beta = params$beta, 
                               kappa= params$kappa,
                               delta = params$delta,
                               phi = params$phi,
                               theta = params$theta, 
                               mu = params$mu,
                               p = params$p,
                               c = params$c,
                               dF = params$dF)
      y <- x$run(t = observation.times)
      y <- as.data.frame(y)
    }
    
    if(model.virus == "triangle"){
      
      y <- model_gen_triangle(T_P = params$T_P,
                              T_U = params$T_U,
                              V_P = params$V_P,
                              t = observation.times)
    }
    
    viral_load = y %>%
      #   slice(rep(1, n0)) %>%
      bind_rows(.,y)%>%
      #   mutate(t = seq(0,by=dt,length.out = n()))%>%
      #      mutate(t = seq(0,by=dt,length.out = n()))%>%
      filter(t <= Tmax) %>%
      mutate(logRNA = logRNA_observable(VI))# %>%
  }
  return(viral_load) 
}

make_household <- function(param.list, household.size){
  
  hh.index = sample(param.list$id, size=household.size,replace = FALSE)
  HH.parameters= param.list[hh.index,]
  HH.parameters$status = 'contact'
  HH.parameters$status[1] = 'index'
  
  HH.parameters$t.infection = -1
  HH.parameters$t.infection[1] = 0
  
  return(HH.parameters)
}


