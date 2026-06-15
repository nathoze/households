
source("main.R")

## Plot all the summary statistics for a specific model of transmission and a given household size
## by default simulations are run on 1000 simulations

for(model in c( "doseresponse","logit","constant")){
  for(household.size in c(3,5)){
    for(size.dependency in c("frequency","density")){
      for(fixed.parameters in c("TRUE","FALSE")){
        simulate_households_plot(epidemiological.model = model, household.size = household.size, size.dependency = size.dependency, fixed.parameters = fixed.parameters)
      }
    }
  }
}

# Viral load in the household. We represent the value for each individual
# for all household size, size dependency, and parameters

 
for(model in c( "doseresponse","logit","constant")){
  for(household.size in c(3,5)){
    for(size.dependency in c("frequency","density")){
      for(fixed.parameters in c("TRUE","FALSE")){
        plot_viral_load_distribution(epidemiological.model = model, household.size = household.size, size.dependency = size.dependency, fixed.parameters = fixed.parameters)
        plot_viral_load_distribution_mean_viral_load(epidemiological.model = model, household.size = household.size, size.dependency = size.dependency, fixed.parameters = fixed.parameters)
        plot_viral_load_distribution_mean_envelop(epidemiological.model = model, household.size = household.size, size.dependency = size.dependency, fixed.parameters = fixed.parameters)
      }
    }
  }
} 

# plot the dynamics of viral load, probability of transmission, and cumulative probability of transmission
# this is the person-to-person probability of transmission (not in households)
 
plot_viral_trajectories()

plot_probability_trajectories("doseresponse",fixed.parameters = FALSE)
plot_probability_trajectories("logit",fixed.parameters = FALSE)
plot_probability_trajectories("constant",fixed.parameters = FALSE)

plot_probability_trajectories("doseresponse",fixed.parameters = TRUE)
plot_probability_trajectories("logit",fixed.parameters = TRUE)
plot_probability_trajectories("constant",fixed.parameters = TRUE)

# # Probability that a contact at a given viral load leads to an infection (cf Figure 2 Marc et al. elife)
# Mean predicted probability of transmission stratified by viral load;
# P =( all individual contacts at a given viral load that results in a tranmission)/ (all individual contacts at a given viral load)
# It is obtained from individual-based simulations in HH

 
#for(model in c( "doseresponse","logit","constant")){
for(model in c( "doseresponse")){
  for(household.size in c(2,3,5)){
    # for(size.dependency in c("frequency","density")){
    for(size.dependency in c("density")){
      #    for(fixed.parameters in c("TRUE","FALSE")){
      for(fixed.parameters in c("FALSE")){
        plot_effect_viral_load_on_transmission(epidemiological.model = model, 
                                               household.size = household.size, size.dependency = size.dependency, fixed.parameters = fixed.parameters)
      }
    }
  }
}




## At which titer value do we have an infection ?  What is the proportion of infections that happen at a given viral load ?
 
#for(model in c( "doseresponse","logit","constant")){
  for(model in c( "doseresponse")){
    #for(model in c( "doseresponse")){
  for(household.size in c(2,3,5)){
    # for(size.dependency in c("frequency","density")){
    for(size.dependency in c("density")){
      for(fixed.parameters in c("TRUE","FALSE")){
        #  for(fixed.parameters in c("FALSE")){
        plot_proportion_transmission_per_viral_load(epidemiological.model = model,
                                                    household.size = household.size,
                                                    size.dependency = size.dependency,
                                                    fixed.parameters = fixed.parameters,
                                                    n.simu = 1200)
      }
    }
  }
}





# Dynamics of Viral load and transmission probability (cf Figure 3 Marc et al. elife)
 plot_dynamics_viral_load_and_transmission("logit",fixed.parameters = FALSE)
plot_dynamics_viral_load_and_transmission("constant",fixed.parameters = FALSE)
plot_dynamics_viral_load_and_transmission("doseresponse",fixed.parameters = FALSE)

plot_dynamics_viral_load_and_transmission("logit",fixed.parameters = TRUE)
plot_dynamics_viral_load_and_transmission("constant",fixed.parameters = TRUE) 
plot_dynamics_viral_load_and_transmission("doseresponse",fixed.parameters = TRUE)
 

