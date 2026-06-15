setwd("/Users/nathanaelhoze/Documents/Work/Scientific_Projects/Projet_Orchestra/CT_Model/src")
rm(list = ls())
source("main.R")
source("data/generate_datasets.R")

alpha = 1
vp = 6
n.h = 100
tmp_generate_datasets_sampling_scheme
for(I in 1:3){
  for(gi in c(5)){
    for(ds in c(5)){
      for(sar in c(0.15,0.3,0.45)){
        print(sar)
        tp = gi/(1+4/vp)
        filename= paste0("Households_26082025_gi_", gi,"_incubation_",ds,"_sar_",sar,"_alpha_", alpha, "_NHH_", n.h, "_index_", I,".rds")
        tmp_generate_datasets_sampling_scheme(n.household = n.h,
                                              V_P_star = 6,
                                              T_P_star = tp,
                                              T_U_star = 5,  
                                              T_P_eta =  0.15,
                                              T_U_eta =  0.15,
                                              V_P_eta =   0.1,
                                              Delta_S_star = ds, 
                                              Delta_S_eta = 0.3, 
                                              filename = filename,
                                              insee.size = TRUE,
                                              epid.model = "exp",
                                              size.parameter.alpha = 1, 
                                              prob_symptomatic = 1,
                                              min.offset.symptoms = 0,
                                              max.offset.symptoms = 2,
                                              SAR.target = sar)
        # }
      }
    }
  }
}
