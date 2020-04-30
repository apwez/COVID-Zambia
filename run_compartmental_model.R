## Setup ####

require(socialmixr)
require(magrittr)
require(stringr)
require(reshape2)
require(dplyr)
require(ggplot2)
require(truncnorm)

source("asyptomatic_age.R")

age.limits = c(0,15,35,55,75)
prop_symptomatic <- c(0.12625, 0.3726667, 0.614, 0.6611429, 0.6611429)

delta.t <- 1/1
time <- seq(1,300,by = delta.t)
t_March25 <- as.numeric(as.Date("2020-03-13") - as.Date("2020-03-01"))
t_May25 <- as.numeric(as.Date("2020-05-25") - as.Date("2020-03-01"))
t_June15 <- as.numeric(as.Date("2020-06-15") - as.Date("2020-03-01"))
t_July1 <- as.numeric(as.Date("2020-07-01") - as.Date("2020-03-01"))

nsim <- 10
start_index <- seq(1, nsim*length(time)+1, by = length(time))

all_prelim_info <- setup_seir_model(stoch = TRUE, 
                                    R0 = 2.5, 
                                    c_scale_vec = 1, 
                                    prop_symptomatic = prop_symptomatic,
                                    sd.dw = 0.1)
Ncomp = all_prelim_info$Ncomp
ICs = all_prelim_info$ICs
params = list(C = all_prelim_info$C, 
              W = all_prelim_info$W, 
              beta0 = all_prelim_info$beta0, 
              beta1 = all_prelim_info$beta1, 
              phase = all_prelim_info$phase, 
              mu = all_prelim_info$mu, 
              v = all_prelim_info$v, 
              N=all_prelim_info$N, 
              gamma=all_prelim_info$gamma, 
              sigma = all_prelim_info$sigma, 
              prop_symptomatic=all_prelim_info$prop_symptomatic,
              sd.dw = all_prelim_info$sd.dw)


cnames.allsim <- c('run_index', 'time',
                   paste0("S", 1:Ncomp),
                   paste0("E", 1:Ncomp),
                   paste0("A", 1:Ncomp),
                   paste0("I", 1:Ncomp),
                   paste0("R", 1:Ncomp),
                   paste0("incid_A", 1:Ncomp),
                   paste0("incid_I", 1:Ncomp),
                   "R0",
                   "Reff")

## ---- Scenario A. Moderate Social Distancing with Constant Effectiveness: ---- ####
## This scenario has mildly restrictive social distancing from March 25 until May 25, similar 
## to that observed during the 1918 influenza pandemic (44-65% reduction in transmission)

all_sim <- matrix(NA,1,(Ncomp*7)+4)
colnames(all_sim) <- cnames.allsim

for(n in 1:nsim){
  print(n)
  R0vec <- rep(runif(1, min = 2, max = 3), length(time))
  c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
  
  eff1 <- runif(1, 0.44, 0.65)
  c_scale_mat[t_March25:(t_May25-1), 1:(Ncomp)] <- (1-eff1)
  # 
  # eff2 <- runif(1, 0.658, 0.858) # MD model uses uniform draws for all effectiveness params
  # c_scale_mat[t_March30:t_July1, 1:(Ncomp-2)] <- (1-eff2)
  
  tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                              R0vec = R0vec, Ncomp = Ncomp, 
                              ICs = ICs, params = params, 
                              time = time, delta.t = delta.t, 
                              c_scale_mat = c_scale_mat)
  run_index = rep(n, nrow(tmp))
  tmp <- cbind(run_index, tmp)
  all_sim <- rbind(all_sim, tmp)
}
all_sim <- all_sim[-1,]
write.csv(all_sim, file="scenarioA_moderate.csv")







# 
# 
# ## ---- Scenario B. Moderate Social Distancing with Degrading Effectiveness --- ####
# ## This scenario has mildly restrictive social distancing from March 13-29, 
# ## similar to that observed during the 1918 influenza pandemic, followed by 
# ## a statewide stay-at-home policy from March 30 - July 1, similar to that in
# ## effect in London, UK. During this stay-at-home period, the effectiveness 
# ## of social distancing decays by 10% every 2 weeks.
# 
# all_sim <- matrix(NA,1,(Ncomp*7)+4)
# colnames(all_sim) <- cnames.allsim
# 
# for(n in 1:nsim){
#   R0vec <- rep(runif(1, min = 2, max = 3), length(time))
#   c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
#   
#   eff1 <- runif(1, 0.44, 0.65)
#   c_scale_mat[t_March13:(t_March30-1), 1:(Ncomp-2)] <- (1-eff1)
#   
#   eff2 <- rtruncnorm(1, a = 0.658, b = 0.858, mean = 0.7615, sd=0.2)
#   c_scale_mat[t_March30:(t_April13-1), 1:(Ncomp-2)] <- 1-eff2
#   c_scale_mat[t_April13:(t_April27-1), 1:(Ncomp-2)] <- 1-eff2*0.9
#   c_scale_mat[t_April27:(t_May11-1), 1:(Ncomp-2)] <- 1-eff2*0.8 #this is how the MD model was run - not successive 10% reductions
#   c_scale_mat[t_May11:(t_May25-1), 1:(Ncomp-2)] <- 1-eff2*0.7
#   c_scale_mat[t_May25:(t_June15-1), 1:(Ncomp-2)] <- 1-eff2*0.6
#   c_scale_mat[t_June15:t_July1, 1:(Ncomp-2)] <- 1-eff2*0.5
#   
#   tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
#                               R0vec = R0vec, Ncomp = Ncomp, 
#                               ICs = ICs, params = params, 
#                               time = time, delta.t = delta.t, 
#                               c_scale_mat = c_scale_mat)
#   run_index = rep(n, nrow(tmp))
#   tmp <- cbind(run_index, tmp)
#   all_sim <- rbind(all_sim, tmp)
# }
# all_sim <- all_sim[-1,]
# write.csv(all_sim, file="output_20200413/scenarioB_mod_degrading.csv")
# 
# ## ---- Scenario C. Moderate Social Distancing Maintained when Schools Re-open April 27. ---- ####
# ## This scenario has mildly restrictive social distancing from March 13-29, 
# ## similar to that observed during the 1918 influenza pandemic, followed by 
# ## a statewide stay-at-home policy from March 30 - July 1, similar to that in 
# ## effect in London, UK. The effectiveness of social distancing decays by 10% 
# ## every 2 weeks and is further reduced 18% when Maryland schools are open from 
# ## April 27 - June 15.
# 
# all_sim <- matrix(NA,1,(Ncomp*7)+4)
# colnames(all_sim) <- cnames.allsim
# 
# for(n in 1:nsim){
#   R0vec <- rep(runif(1, min = 2, max = 3), length(time))
#   c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
#   
#   eff1 <- runif(1, 0.44, 0.65)
#   c_scale_mat[t_March13:(t_March30-1), 1:(Ncomp-2)] <- (1-eff1)
#   
#   eff2 <- rtruncnorm(1, a = 0.658, b = 0.858, mean = 0.7615, sd=0.2)
#   c_scale_mat[t_March30:(t_April13-1), 1:(Ncomp-2)] <- 1-eff2
#   c_scale_mat[t_April13:(t_April27-1), 1:(Ncomp-2)] <- 1-eff2*0.9
#   c_scale_mat[t_April27:(t_May11-1), 1:(Ncomp-2)] <- 1- (eff2*0.8 - 0.18) # school closures were constant, absolute 18% decrease in effectiveness
#   c_scale_mat[t_May11:(t_May25-1), 1:(Ncomp-2)] <- 1- (eff2*0.7 - 0.18)
#   c_scale_mat[t_May25:(t_June15-1), 1:(Ncomp-2)] <- 1- (eff2*0.6 - 0.18)
#   c_scale_mat[t_June15:t_July1, 1:(Ncomp-2)] <- 1-eff2*0.5
#   
#   tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
#                               R0vec = R0vec, Ncomp = Ncomp, 
#                               ICs = ICs, params = params, 
#                               time = time, delta.t = delta.t, 
#                               c_scale_mat = c_scale_mat)
#   run_index = rep(n, nrow(tmp))
#   tmp <- cbind(run_index, tmp)
#   all_sim <- rbind(all_sim, tmp)
# }
# all_sim <- all_sim[-1,]
# write.csv(all_sim, file="output_20200413/scenarioC_mod_degrading_scl.csv")
# 
# 
# ## ---- Scenario D. All restrictions lifted April 27 --- ####
# ## This scenario has mildly restrictive social distancing from March 13-29, 
# ## similar to that observed during the 1918 influenza pandemic, followed by 
# ## a statewide stay-at-home policy from March 30 - April 27, similar to that 
# ## in effect in London, UK. After April 28, Maryland schools reopen and all 
# ## interventions are lifted. During the stay-at-home period, the effectiveness 
# ## of social distancing decays by 10% every 2 weeks.
# 
# all_sim <- matrix(NA,1,(Ncomp*7)+4)
# colnames(all_sim) <- cnames.allsim
# 
# for(n in 1:nsim){
#   R0vec <- rep(runif(1, min = 2, max = 3), length(time))
#   c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
#   
#   eff1 <- runif(1, 0.44, 0.65)
#   c_scale_mat[t_March13:(t_March30-1), 1:(Ncomp-2)] <- (1-eff1)
#   
#   eff2 <- rtruncnorm(1, a = 0.658, b = 0.858, mean = 0.7615, sd=0.2)
#   c_scale_mat[t_March30:(t_April13-1), 1:(Ncomp-2)] <- 1-eff2
#   c_scale_mat[t_April13:(t_April27-1), 1:(Ncomp-2)] <- 1-eff2*0.9
#   
#   tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
#                               R0vec = R0vec, Ncomp = Ncomp, 
#                               ICs = ICs, params = params, 
#                               time = time, delta.t = delta.t, 
#                               c_scale_mat = c_scale_mat)
#   run_index = rep(n, nrow(tmp))
#   tmp <- cbind(run_index, tmp)
#   all_sim <- rbind(all_sim, tmp)
# }
# all_sim <- all_sim[-1,]
# write.csv(all_sim, file="output_20200413/scenarioD_total_reopen.csv")
