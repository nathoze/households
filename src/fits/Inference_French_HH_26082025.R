setwd("/Users/nathanaelhoze/Documents/Work/Scientific_Projects/Projet_Orchestra/CT_Model/src")
rm(list = ls())
source("main.R")
set.seed(1954)
# 
fit.epid.model = "exp"
alpha=1
vp=6
gi = 5
ds= 5
tp = gi/(1+4/vp)
########### Fit with the full dataset + reconstruction of the chains of transmission  #############

n.h = 100

for(sar in c(0.15, 0.30, 0.45)){
  print(sar)
  #  for(I in 1:100){
  for(I in 1:100){
    print(I)
    # filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
    # data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
    
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    
    output.file.ct = paste0('ct_',filename)
    filename.1 = paste0("../results/results_",filename,".rds")
    
    if(!file.exists(filename.1)){
      if(fit.epid.model == "exp"){
        res = run_inference(xdata, init, rstan_model,iter_warmup=1500)
      }
      
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1 )
    }
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = fit.epid.model)
        save('CT', file = paste0('../results/',output.file.ct) )
      }  
    }
  }
}



########### Fit with the symptoms only + reconstruction of the chains of transmission  #############


n.h = 100

for(sar in c(0.15, 0.30,0.45)){
  print(sar)
  #  for(I in 1:100){
  for(I in 1:100){
    print(I)
    # filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
    # data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
    
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    
    output.file.ct = paste0('ct_symptoms',filename)
    filename.1 = paste0("../results/symptoms_results_",filename,".rds")
    
    if(!file.exists(filename.1)){
      if(fit.epid.model == "exp"){
        res=run_inference_symptoms(xdata,init, rstan_model = rstan_model_symptoms_gammma,iter_warmup = 1500)
      }
      
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1 )
    }
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = reconstruct_transmission_chain_symptoms(xdata, Chains = Chains,hh.index  )
        save('CT', file = paste0('../results/',output.file.ct) )
      }  
    }
  }
}


######## 28 OCT 2025 ##########
# Fit of the constant model to the 100 simulations of 100 households, only for SAR = 30% -----

n.h = 100
sar = 0.45 


n.h = 100

for(sar in c(0.45)){
  print(sar)
  #  for(I in 1:100){
  for(I in 1:100){
    print(I)
    # filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
    # data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
    
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    
    output.file.ct = paste0('ct_constant_',filename)
    filename.1 = paste0("../results/results_constant_",filename,".rds")
    
    if(!file.exists(filename.1)){
      res = run_inference(xdata, init, rstan_model_constant,iter_warmup=1500)
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1 )
    }
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = "constant")
        save('CT', file = paste0('../results/',output.file.ct) )
      }  
    }
  }
}


########### Fit with the   dataset  with only positivity results + reconstruction of the chains of transmission  + NO SYMPTOMS FOR THE CONTACTS #############


n.h = 100

# for(sar in c(0.15,0.30, 0.45)){
for(sar in c(0.45)){
  print(sar)
  #  for(I in 1:100){
  for(I in 2:100){
    print(I)
    # filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
    # data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
    
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    contact.id = which(xdata$col_infected>1)
    xdata$symptomatic[contact.id] =0
    
    output.file.ct = paste0('ct_no_symptoms_',filename)
    filename.1 = paste0("../results/results_no_symptoms_",filename,".rds")
    
    if(!file.exists(filename.1)){
      # if(fit.epid.model == "exp"){
      res = run_inference(xdata, init, rstan_model = rstan_model_positivity_3,iter_warmup=1500)
      # }
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1 )
    }
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = reconstruct_transmission_chain(xdata, Chains = Chains,hh.index, model = fit.epid.model)
        save('CT', file = paste0('../results/',output.file.ct) )
      }  
    }
  }
}






########### Fit with the full dataset + reconstruction of the chains of transmission  + NO SYMPTOMS FOR THE CONTACTS #############




# 
# 
# n.h = 100
# 
# for(sar in c(0.15, 0.30, 0.45)){
#   print(sar)
#   for(I in 1:100){
#     print(I)
#     
#     filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
#     data.filename = paste0("../data/simulated_data/",filename)
#     
#     load(data.filename)
#     
#     contact.id = which(xdata$col_infected>1)
#     xdata$symptomatic[contact.id] =0
#     
#     
#     
#     output.file.ct = paste0('ct_no_symptoms_',filename)
#     filename.1 = paste0("../results/results_no_symptoms_",filename,".rds")
#     
#     if(!file.exists(filename.1)){
#       if(fit.epid.model == "exp"){
#         res = run_inference(xdata, init, rstan_model = rstan_model_exp,iter_warmup=1500)
#       }
#       
#       Chains = res$Chains
#       DIC = res$DIC
#       pD = res$pD
#       print(DIC)
#       save('Chains','DIC','pD', file =filename.1 )
#     }
#     
#     if(file.exists(filename.1)){
#       load(file = filename.1)
#       
#       CT=NULL
#       for(hh.index in 1:xdata$n_households){
#         print(hh.index)
#         CT[[hh.index]] = reconstruct_transmission_chain(xdata, Chains = Chains,hh.index, model = fit.epid.model)
#         save('CT', file = paste0('../results/',output.file.ct) )
#       }  
#     }
#   }
# }
# 
# 
# 


### xxxxxx 31 mars  -----
### all viral loads, no symptoms


for(sar in c(0.45)){
  print(sar)
  #  for(I in 1:100){
  for(I in 1:100){
    print(I)
    # filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
    # data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
    
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    contact.id = which(xdata$col_infected>1)
    xdata$symptomatic[contact.id] =0
    
    
    output.file.ct = paste0('ct_viral_loads_no_symptoms_',filename)
    filename.1 = paste0("../results/results_viral_loads_no_symptoms_",filename,".rds")
    
    if(!file.exists(filename.1)){
      if(fit.epid.model == "exp"){
        res = run_inference(xdata, init, rstan_model = rstan_model_exp,iter_warmup=1500)
      }
      
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1 )
    }
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = reconstruct_transmission_chain(xdata, Chains = Chains,hh.index, model = fit.epid.model)
        save('CT', file = paste0('../results/',output.file.ct) )
      }  
    }
  }
}
  



### Fit the model for the data with positivity ----- 

n.h = 100

for(sar in c(0.15, 0.30, 0.45)){
  print(sar)
  for(I in 1:100){
    print(I)
    
    filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
    data.filename = paste0("../data/simulated_data/",filename)
    load(data.filename)
    
    output.file.ct = paste0('ct_',filename)
    filename.1 = paste0("../results/results_positivity/results_positivity_",filename,".rds")
    
    
    
    init = function(){
      list(logm = runif(n=1, min=-7, max=0),
           alpha =  runif(n=1, min=0, max=1),
           # log_Delta_S_mu =  2,
           # Delta_S_sd =  0.5,
           # log_T_P_mu = 0,
           # log_T_U_mu =  0,
           # log_V_P_mu = 0,        
           # log_T_P_sd =  1,
           # log_T_U_sd = 1,
           # log_V_P_sd =  1,
           # T_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
           # T_U_eta = as.array(rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
           # V_P_eta = as.array( rnorm( n = xdata$n_infected_total, mean=0, sd=0.0)),
           #  TI_tmp =  as.array(runif( n = xdata$n_infected_total,  min = -10, max=30)),
           # TI_tmp_positive =  as.array(runif( n = xdata$n_infected_total,  min = 0, max=10)),
           sigma_noise = 0.8)
    }
    
    
    
    
    res = run_inference_positivity(xdata, init, rstan_model_positivity_2,iter_warmup=1500)
    Chains = res$Chains
    DIC = res$DIC
    pD = res$pD
    print(DIC)
    save('Chains','DIC','pD', file =filename.1 )
    
    if(file.exists(filename.1)){
      load(file = filename.1)
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = reconstruct_transmission_chain_positivity(xdata, Chains = Chains,hh.index)
        save('CT', file = paste0('../results/results_positivity/',output.file.ct) )
      }  
    }
  }
}



rstan_model_positivity_2 = stan_model("models/Stan_scripts/Positivity_Model_2.stan") 
res = run_inference_positivity(xdata, init, rstan_model_positivity_2,iter_warmup=1500)
Chains = res$Chains
DIC = res$DIC
pD = res$pD
print(DIC)



############ OLD ##############


########### 3 #############

for(gi in c(3,5,7)){
  for(ds in c(3,5,7)){
    for(sar in c(0.15,0.3,0.45)){
      tp = gi/(1+4/vp)
      
      filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
      data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
      load(data.filename)
      
      
      filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
      output.file.ct = paste0('ct_',filename,"_fit_", fit.epid.model,".rds")
      
      filename.1 = paste0("../results/results_",filename,"_exp_fit_", fit.epid.model, ".rds")
      # if(!file.exists(filename.1)){
      if(file.exists(filename.1)){
        load(file = filename.1)
        
        CT=NULL
        for(hh.index in 1:xdata$n_households){
          print(hh.index)
          CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = fit.epid.model)
          save('CT', file = paste0('../results/',output.file.ct) )
        }  
      }
    }
  }
}




########### 4 #############

for(gi in c(3,5,7)){
  for(ds in c(3,5,7)){
    for(sar in c(0.15,0.3,0.45)){
      tp = gi/(1+4/vp)
      
      filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
      data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
      load(data.filename)
      output.file.ct = paste0('ct_',filename,"_fit_", fit.epid.model,"_positive.rds")
      
      filename.1 = paste0("../results/results_",filename,"_exp_fit_", fit.epid.model, "_positive.rds")
      if(!file.exists(filename.1)){
        # if(fit.epid.model == "exp"){
        res = run_inference(xdata, init, rstan_model_positive,iter_warmup=1500)
        # }
        # if(fit.epid.model == "constant"){
        #   res = run_inference(xdata, init, rstan_model_constant,iter_warmup=1500)
        # }
        Chains = res$Chains
        DIC = res$DIC
        pD = res$pD
        print(DIC)
        save('Chains','DIC','pD', file =filename.1 )
        # }
      }
    }
  }
} 
########### 5 #############

for(gi in c(3,5,7)){
  for(ds in c(3,5,7)){
    for(sar in c(0.15,0.3,0.45)){
      tp = gi/(1+4/vp)
      
      filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
      data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
      load(data.filename)
      output.file.ct = paste0('ct_',filename,"_fit_", fit.epid.model,"_positive.rds")
      
      filename.1 = paste0("../results/results_",filename,"_exp_fit_", fit.epid.model, "_positive.rds")
      # if(!file.exists(filename.1)){
      if(file.exists(filename.1)){
        load(file = filename.1)
        
        CT=NULL
        for(hh.index in 1:xdata$n_households){
          print(hh.index)
          CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = fit.epid.model)
          save('CT', file = paste0('../results/',output.file.ct) )
        }  
      }
    }
  }
}




########### 6+7 #############

for(gi in c(3,5,7)){
  for(ds in c(3,5,7)){
    for(sar in c(0.15,0.3,0.45)){
      
      filename= paste0("200_households_11042025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
      data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
      load(data.filename)
      output.file.ct = paste0('ct_symptoms_',filename,"_fit_", fit.epid.model,".rds")
      
      filename.1 = paste0("../results/results_symptoms_",filename,"_exp_fit_", fit.epid.model, ".rds")
      
      res=run_inference_symptoms(xdata,init, rstan_model = rstan_model_symptoms_gammma,iter_warmup = 1500)
      
      Chains = res$Chains
      DIC = res$DIC
      pD = res$pD
      print(DIC)
      save('Chains','DIC','pD', file =filename.1  )
      
      CT=NULL
      for(hh.index in 1:xdata$n_households){
        print(hh.index)
        CT[[hh.index]] = chain_transmission_reconstruction_symptoms(xdata, Chains = Chains,hh.index)
        #  save('CT', file = paste0('../results/',output.file.ct) )
      }
      
    }
  }
} 
# 






# 
# #### Understanding why the fit with the constant model is so bad
# 
# gi=5
# ds=5
# sar= 0.30
# 
# filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
# load(data.filename)
# 
# 
# filename= paste0("200_households_05032025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# output.file.ct = paste0('ct_',filename,"_fit_exp.rds")
# filename.1 = paste0("../results/results_",filename,"_exp_fit_exp.rds")
# 
# res = run_inference(xdata, init, rstan_model_essai,iter_warmup=1500)
# 
# Chains = res$Chains
# DIC = res$DIC
# pD = res$pD
# print(DIC)
# save('Chains','DIC','pD', file =filename.1 )
# 
# 
# 
# 
# g1 = MCMC_plot_predicted_viral_load_hhindex(Chains = Chains,xdata, Xlim = 20, 149)
# 
# 
# 
# 
# 
# gi=5
# ds=5
# sar= 0.30
# alpha=1
# 
# filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
# load(data.filename)
# 
# 
# filename= paste0("200_households_05032025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# output.file.ct = paste0('ct_',filename,"_fit_constant.rds")
# filename.1 = paste0("../results/results_",filename,"_exp_fit_constant.rds")
# 
# res = run_inference(xdata, init, rstan_model_constant,iter_warmup=1500)
# 
# Chains = res$Chains
# DIC = res$DIC
# pD = res$pD
# print(DIC)
# save('Chains','DIC','pD', file =filename.1 )
# 
# print(g1)
# CT=NULL
# for(hh.index in 1:xdata$n_households){
#   print(hh.index)
#   CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = "constant")
#   save('CT', file = paste0('../results/',output.file.ct) )
# }  
# 
# 
# 
# 
# 
# 
# 
# rstan_model_constant_essai = stan_model("models/Stan_scripts/Constant_Model_essai.stan") 
# 
# gi=5
# ds=5
# sar= 0.30
# alpha=1
# 
# filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
# load(data.filename)
# 
# 
# filename= paste0("200_households_05032025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# output.file.ct = paste0('ct_',filename,"_fit_constant.rds")
# filename.1 = paste0("../results/results_",filename,"_exp_fit_noinfectionmodel.rds")
# 
# res = run_inference(xdata, init, rstan_model_constant_essai,iter_warmup=1500)
# 
# Chains = res$Chains
# DIC = res$DIC
# pD = res$pD
# print(DIC)
# save('Chains','DIC','pD', file =filename.1 )
# 
# g1 = MCMC_plot_predicted_viral_load_hhindex(Chains = Chains,xdata, Xlim = 20, 149)
# MCMC_plot_predicted_viral_load_hhindex(Chains = Chains,xdata, Xlim = 20, 67)
# MCMC_plot_predicted_viral_load_hhindex(Chains = Chains,xdata, Xlim = 20, 68)
# MCMC_plot_predicted_viral_load_hhindex(Chains = Chains,xdata, Xlim = 20, 67)
# 
# plot(exp(Chains$log_Delta_S_mu))
# 
# 
# print(g1)
# CT=NULL
# for(hh.index in 1:xdata$n_households){
#   print(hh.index)
#   CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index, model = "constant")
#   save('CT', file = paste0('../results/',output.file.ct) )
# }  
#  
# 
# 
# load("../results/results_200_households_05032025_gi_7_incubation_7_sar_0.45_alpha_1_exp_fit_noinfectionmodel.rds")
# 
# # #### check files ----
# # path="/Users/nathanaelhoze/Documents/Work/Scientific_Projects/Projet_Orchestra/CT_Model/"
# # fit.epid.model="exp"
# # alpha=1
# # df=NULL
# # I=0
# # for(gi in c(3,5,7)){
# #   for(ds in c(3,5,7)){
# #     for(sar in c(0.15,0.3,0.45)){
# #       
# #       filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# #       data.filename = paste0(path,"data/simulated_data/",filename,"_exp.rds")
# #       load(data.filename)
# #       output.file.ct = paste0(path,"results/ct_",filename,"_fit_", fit.epid.model,".rds")
# #       
# #       filename.1 = paste0(path,"results/results_",filename,"_exp_fit_", fit.epid.model, ".rds")
# #       load(file = filename.1)
# #       print(DIC)
# #       print(pD)
# #       I=I+1
# #       CT=0
# #       if(file.exists(output.file.ct)){
# #         CT=1
# #       }
# #       df = rbind(df, data.frame(index = I, gi=gi, ds=ds, sar=sar, pD=pD, DIC=DIC, CT=CT))
# #     }
# #   }
# # } 
# # 
# # 
# # D= df %>% filter(pD<100)
# # 
# # rstan_model = stan_model("models/Stan_scripts/Exp_Model.stan") ### Infer the time of symptoms, use it to improve the inference of the infection time.
# # 
# # 
# # for(a in 1:nrow(D)){
# #   
# #   print(a)
# #   gi = D$gi[a]
# #   ds = D$ds[a]
# #   sar = D$sar[a]
# #   
# #   
# #   filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# #   data.filename = paste0("../data/simulated_data/",filename,"_exp.rds")
# #   load(data.filename)
# #   output.file.ct = paste0('ct_',filename,"_fit_", fit.epid.model,".rds")
# #   
# #   filename.1 = paste0("../results/results_",filename,"_exp_fit_", fit.epid.model, ".rds")
# #   if(fit.epid.model == "exp"){
# #     res = run_inference(xdata, init, rstan_model,iter_warmup=1500)
# #   }
# #   if(fit.epid.model == "constant"){
# #     res = run_inference(xdata, init, rstan_model_constant,iter_warmup=1500)
# #   }
# #   Chains = res$Chains
# #   DIC = res$DIC
# #   pD = res$pD
# #   print(DIC)
# #   save('Chains','DIC','pD', file =filename.1 )
# #   
# #   
# # }
# # 
# # 
# # 
# # 
# # 
# # path="/Users/nathanaelhoze/Documents/Work/Scientific_Projects/Projet_Orchestra/CT_Model/"
# # fit.epid.model="exp"
# # alpha=1
# # df=NULL
# # I=0
# # for(gi in c(3,5)){
# #   for(ds in c(3,5,7)){
# #     for(sar in c(0.45)){
# #       
# #       filename= paste0("200_households_26022025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha)
# #       data.filename = paste0(path,"data/simulated_data/",filename,"_exp.rds")
# #       load(data.filename)
# #       output.file.ct = paste0(path,"results/ct_",filename,"_fit_", fit.epid.model,".rds")
# #       
# #       filename.1 = paste0(path,"results/results_",filename,"_exp_fit_", fit.epid.model, ".rds")
# #       load(file = filename.1)
# #       print(DIC)
# #       print(pD)
# #       
# #       
# #       I=I+1
# #       CT=0
# #       if(!file.exists(output.file.ct) & pD>100){
# # 
# #           CT=NULL
#           for(hh.index in 1:xdata$n_households){
#             print(hh.index)
#             CT[[hh.index]] = chain_transmission_reconstruction(xdata, Chains = Chains,hh.index)
#             save('CT', file = output.file.ct) 
#           }
#         # df = rbind(df, data.frame(index = I, gi=gi, ds=ds, sar=sar, pD=pD, DIC=DIC, CT=CT))
#         
#       }
#     }
#   }
# } 
# 
# 
