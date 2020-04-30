rm(list=ls())
### code adapted from Alex Becker's COVID age structured discrete time model.
### This code will incorporate age structured mixing between groups and produce estimates of the number of 
### infected individuals per day (and cumulative incident infections). We have adapted the code to add in 
### separate healthcare workers and homeless populations. 

#setwd('~/Dropbox/COVID-BaltimoreCity//')

require(socialmixr)
require(magrittr)
require(stringr)
require(reshape2)
require(dplyr)
require(ggplot2)
require(truncnorm)

### Functions: this document has 5 different functions to clean up the age structured data, make a mixing matrix, 
### rescale the mixing matrix, run the discrete time simulation, and set up all of the parameters/etc. for the discrete time simulation.
# make_age_structure_matrix(age_data_file, homeless_n, healthcare_n)
# make_polymod_matrix()
# rescale_age_matrix(Ncomp, W, Zam_pop)
# sair_step(stoch = F, Ncomp, ICs, params, time, delta.t)
# setup_seir_model(stoch)

### general things
# Zam_pop: the population of Zambia by age
# W: age specific mixing rates
# Ncomp: number of compartments (ages + healthcare + homeless)
# SAIR with asymtompatics separate from symtomatics
# currently there is a bit to estimate severity that will not be needed later
# gamma <- 1/6.5 ## infectious period? make a range?
# R0 <- 2.5 ## make a range?
# prop_symptomatic <- 0.20 ## will need to update/change

theme_set(theme_classic(base_size=12))

## ---- loading demographic data ---- ####
### load [location] age structured demographic data. This code assumes that the first row will be the total 
### population and the following rows are the ages broken down by age groups.
## If you do not want homeless or healthcare workers, just set these values to zero
## assumes the additional of healthcare workers and homeless does not decrease from their respective age group

make_age_structure_matrix <- function(age_data_file){
  ## Pull out population data
  n_age_classes <- nrow(age_data_file)-1
  Zam_pop <- age_data_file[2:nrow(age_data_file),]$Estimate
  Ncomp <- length(Zam_pop) 
  return(list(Ncomp = Ncomp, Zam_pop = Zam_pop))  
}

## ---- making polymod matrix ---- ####
### sets up polymod matrix using data from Zambia -- it will make separate age categories that are set up to minimize missing age data. However, children are not in the survey so instead the UK children value is used. 

make_polymod_matrix <- function(age.limits = c(0,15,35,55,75)){
  ## setup polymod matrix
  ## use age classes in the data
  ## for now, hard code it into the function, but can change later
  # zambia_survey <- get_survey('https://doi.org/10.5281/zenodo.2548692') ## - does not have any children 
  # saveRDS(zambia_survey, 'zambia.rds')
  polymod_zambia <- readRDS('zambia.rds')
  age_mix <- contact_matrix(polymod_zambia, age.limits = age.limits, estimated.contact.age = 'mean', missing.contact.age = 'sample')$matrix
  
  ## currently using UK data for the youngest age group 
  data(polymod)
  age_mix_uk <- contact_matrix(polymod, countries = "United Kingdom", age.limits = age.limits)$matrix
  age_mix[1,] = age_mix_uk[1,]

  W <- matrix(NA,nrow(age_mix), ncol(age_mix))
  
  rownames(W) <- c(rownames(age_mix))
  colnames(W) <- c(colnames(age_mix))
  W[1:nrow(age_mix),1:ncol(age_mix)] <- age_mix
  
  return(W)
}

## ---- rescaling age matrix ---- ####
## we also need to adjust our matrix to make sure R0 = beta/gamma
## we want the max eigen value of our input matrix to be 1 such that we can say R0 is just beta / gamma
## the following formula gives the R0 calculation for age strucutred matrix models
## ref : http://www.sherrytowers.com/towers_feng_2012.pdf page 242 right side (matrix is called C_ij)

rescale_age_matrix <- function(Ncomp, W, Zam_pop, c_scale_vec){
  A <- matrix(0,Ncomp,Ncomp)
  for (ii in 1:Ncomp){
    for (jj in 1:Ncomp){
      A[ii,jj] <- W[ii,jj]*Zam_pop[ii]/Zam_pop[jj]
    }
  }
  ## compute spectral radius / scaling parameter
  r <- eigen(A)
  lam <- r$values
  alpha <- max(Re(lam))
  W <- W / alpha
  ## now the matrix is rescaled have R0  = 1, so beta0 can be scaled to be real transmission
  C <- matrix(c_scale_vec,nrow(W),ncol(W)) ## a special contact matrix to be used to rescale health facility worker contact rates - when set to 1, it is turned off 
  return(list(W = W, C = C))
}

## ---- main S(A)IR function ---- ####
sair_step <- function(stoch = F, stoch.init = F, sNcomp, ICs, params, time, delta.t){
  C = params$C
  W = params$W
  beta0 = params$beta0
  beta1 = params$beta1
  phase = params$phase
  mu = params$mu
  v = params$v
  N=params$N
  sigma = params$sigma
  gamma=params$gamma
  prop_symptomatic=params$prop_symptomatic
  sd.dw = params$sd.dw
  
  ## set up a matrix to store values in by variable and time
  ## each X[it,] is the variable at one hour
  x <- matrix(NA,length(time),Ncomp * 7)
  
  ## set initial conditions
  if(stoch.init){
    Ninit <- sample(10:60, 1)
    Ninit_byage <- rmultinom(1, Ninit, prob = N/sum(N))[,1]
    Ninit_asy <- round(Ninit_byage * (1-prop_symptomatic))
    Ninit_sym <- Ninit_byage - Ninit_asy
    ICs <- c(S = N, 
             E = rep(0, Ncomp),
             A = Ninit_asy,
             I = Ninit_sym,
             R = rep(0, Ncomp),
             incid_A = rep(0, Ncomp),
             incid_I = rep(0, Ncomp))
    x[1,] <- round(ICs)
    
  }else{ x[1,] <- round(ICs) }
  
  S <- x[,1:Ncomp]; ## susceptible individuals
  E <- x[,(Ncomp+1):(2*Ncomp)]; ## exposed individuals 
  A <- x[,(2*Ncomp+1):(3*Ncomp)]; ## asymptomatic individuals
  I <- x[,(3*Ncomp+1):(4*Ncomp)];## symp individuals
  R <- x[,(4*Ncomp+1):(5*Ncomp)] ## recovered individuals
  
  ## incidence
  incid_A <- x[,(5*Ncomp+1):(6*Ncomp)];
  incid_I <- x[,(6*Ncomp+1):(7*Ncomp)];
  
  ## seasonal transmission
  seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365 - phase))
  
  for(it in 1:(length(time) - 1)){
    #  WI <- C%*%W%*%(A[it,] + I[it,])
    WI <- (C*W)%*%(A[it,] + I[it,])
    
    WI[!is.finite(WI)] <- 0
    births <-rep(0,Ncomp)
    births[1] <- mu
    deaths <- rep(v,Ncomp)
    
    ## add stochasticity to FOI
    if(stoch == T){
      dw <- rtruncnorm(Ncomp, a=0, mean = 1, sd = sd.dw)
    }else{
      dw <- 1
    }
    ## declare transitions in model
    foi_prob <- 1 - exp( - seas[it] * WI/N * dw * delta.t)
    exposed_prob <- 1 - exp( - sigma * delta.t)
    inf_prob <- 1 - exp( - gamma * delta.t)
    death_prob <- 1 - exp( - deaths * delta.t)
    
    ## stochastic formulation of the model
    if(stoch == T){
      new_exp <- rbinom(n = Ncomp, size = round(S[it,]), prob = foi_prob)
      new_inf <- rbinom(n = Ncomp, size = round(E[it,]) , prob = exposed_prob)
      new_infA <- round( (1-prop_symptomatic)*new_inf )
      new_infI <- new_inf - new_infA
      new_rec_A <- rbinom(n = Ncomp, size = round(A[it,]), prob = inf_prob)
      new_rec_I <- rbinom(n = Ncomp, size = round(I[it,]), prob = inf_prob)
      
      S[it + 1, ] <- S[it,] +  births*delta.t - new_exp - rbinom(n = Ncomp, size = round(S[it,]), prob = death_prob)
      E[it + 1, ] <- E[it,] +  new_exp - new_inf - rbinom(n = Ncomp, size = round(E[it,]), prob = death_prob )
      A[it + 1, ] <- A[it,] +  new_infA - new_rec_A - rbinom(n = Ncomp, size = round(A[it,]), prob = death_prob)
      I[it + 1, ] <- I[it,] +  new_infI - new_rec_I - rbinom(n = Ncomp, size = round(I[it,]), prob = death_prob)
      R[it + 1, ] <- R[it,] +  new_rec_I + new_rec_A - rbinom(n = Ncomp, size = round(R[it,]), prob = death_prob)
      
      ## make incidence the new number of daily individuals becoming infected
      incid_A[it, ] <- new_infA
      incid_I[it, ] <- new_infI
    }
    
    ## deterministic equations
    if(stoch == F){
      S[it + 1, ] <- S[it,] + delta.t * (births - seas[it] * WI * S[it,] * dw / N  - deaths*S[it,])
      E[it + 1, ] <- E[it,] + delta.t * (seas[it] * WI * S[it,] * dw / N - deaths*E[it,] - sigma*E[it,])
      A[it + 1, ] <- A[it,] + delta.t * ( (1 - prop_symptomatic)*sigma*E[it,]   - A[it,]*(gamma - deaths))
      I[it + 1, ] <- I[it,] + delta.t * (  prop_symptomatic*sigma*E[it,] - I[it,]*(gamma - deaths) )
      R[it + 1, ] <- R[it,] + delta.t * (A[it,]*gamma+ I[it,]*gamma - R[it,]* deaths)
      incid_A[it,] <-  (1 - prop_symptomatic)*(seas[it] * WI * S[it,] * dw / N)
      incid_I[it,] <- prop_symptomatic*(seas[it] * WI * S[it,] * dw / N)
    }
  }
  out <- data.frame(cbind(time,S,E,A,I,R,incid_A, incid_I))
  names(out) <- c('time',names(ICs))
  ## output is the number in each class per time point per age-category+homeless+healthcare workers
  return(out)
}

## ---- main SEAIR function with variable R0 input ---- ####
sair_step_variableR0 <- function(stoch = F, stoch.init = F, R0vec, Ncomp, ICs, params, time, delta.t, init.min = 10, init.max=60, init.dist = NULL, c_scale_mat=NULL){
  
  C = params$C
  W = params$W
  beta0 = params$beta0
  beta1 = params$beta1
  phase = params$phase
  mu = params$mu
  v = params$v
  N=params$N
  sigma = params$sigma
  gamma=params$gamma
  prop_symptomatic=params$prop_symptomatic
  sd.dw <- params$sd.dw
  
  ## set up a matrix to store values in by variable and time
  ## each X[it,] is the variable at one hour
  x <- matrix(NA,length(time),Ncomp * 7)
  
  ## set initial conditions
  if(stoch.init){
    Ninit <- sample(init.min:init.max, 1)
    if(is.null(init.dist)){pinit <- N / sum(N)}else{pinit <- init.dist}
    Ninit_byage <- rmultinom(1, Ninit, prob = pinit)[,1]
    Ninit_asy <- round(Ninit_byage * (1-prop_symptomatic))
    Ninit_sym <- Ninit_byage - Ninit_asy
    ICs <- c(S = N, 
             E = rep(0, Ncomp),
             A = Ninit_asy,
             I = Ninit_sym,
             R = rep(0, Ncomp),
             incid_A = rep(0, Ncomp),
             incid_I = rep(0, Ncomp))
    x[1,] <- round(ICs)
    
  }else{ x[1,] <- round(ICs) }
  
  S <- x[,1:Ncomp]; ## susceptible individuals
  E <- x[,(Ncomp+1):(2*Ncomp)]; ## exposed individuals 
  A <- x[,(2*Ncomp+1):(3*Ncomp)]; ## asymptomatic individuals
  I <- x[,(3*Ncomp+1):(4*Ncomp)];## symp individuals
  R <- x[,(4*Ncomp+1):(5*Ncomp)] ## recovered individuals
  
  ## incidence
  incid_A <- x[,(5*Ncomp+1):(6*Ncomp)]
  incid_I <- x[,(6*Ncomp+1):(7*Ncomp)]
  
  # recalculate beta0 for each R0 value
  beta0 <- R0vec * (gamma + v) * (sigma + v) / sigma 
  seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365 - phase))
  R0 <- vector(length=length(time))
  Re <- vector(length=length(time))
  
  for(it in 1:(length(time) - 1)){
    
    # calculate and store the R0 value for time-specific beta0 - to confirm the new value is correct
    R0.mat <- matrix(0,Ncomp,Ncomp)
    for (i in 1:Ncomp){
      for (j in 1:Ncomp){
        R0.mat[i,j] <- W[i,j]*N[i]/N[j]* beta0[it] * sigma / ( (sigma + v) * (v + gamma))
      }
    }
    R0[it] <- Re(eigen(R0.mat)$values[1])
    
    # proceed with model step, as in sair_step
    # if c_scale_mat is provided, recalculate C matrix
    if(!is.null(c_scale_mat)){
      C <- matrix(c_scale_mat[it,],nrow(W),ncol(W))
    }
    WI <- (C*W)%*%(A[it,] + I[it,])
    WI[!is.finite(WI)] <- 0
    
    # calculate and store the effective R value (if any value of C is non-1)
    if(sum(C!=1)>0){
      Wscal <- C*W
      Re.mat <- matrix(0,Ncomp,Ncomp)
      for (i in 1:Ncomp){
        for (j in 1:Ncomp){
          Re.mat[i,j] <- Wscal[i,j]*N[i]/N[j]* beta0[it] * sigma / ( (sigma + v) * (v + gamma))
        }
      }
      Re[it] <- Re(eigen(Re.mat)$values[1])
    }else{
      Re[it] <- R0[it]
    }
    
    
    births <-rep(0,Ncomp)
    births[1] <- mu
    deaths <- rep(v,Ncomp)
    
    ## add stochasticity to FOI
    if(stoch == T){
      dw <- rtruncnorm(Ncomp, a=0, mean = 1, sd = sd.dw)
    }else{
      dw <- 1
    }
    
    ## declare transitions in model
    foi_prob <- 1 - exp( - seas[it] * WI/N * dw * delta.t)
    exposed_prob <- 1 - exp( - sigma * delta.t)
    inf_prob <- 1 - exp( - gamma * delta.t)
    death_prob <- 1 - exp( - deaths * delta.t)
    
    ## stochastic formulation of the model
    if(stoch == T){
      new_exp <- rbinom(n = Ncomp, size = round(S[it,]), prob = foi_prob)
      new_inf <- rbinom(n = Ncomp, size = round(E[it,]) , prob = exposed_prob)
      new_infA <- round( (1-prop_symptomatic)*new_inf )
      new_infI <- new_inf - new_infA
      new_rec_A <- rbinom(n = Ncomp, size = round(A[it,]), prob = inf_prob)
      new_rec_I <- rbinom(n = Ncomp, size = round(I[it,]), prob = inf_prob)
      
      S[it + 1, ] <- S[it,] +  births*delta.t - new_exp - rbinom(n = Ncomp, size = round(S[it,]), prob = death_prob)
      E[it + 1, ] <- E[it,] +  new_exp - new_inf - rbinom(n = Ncomp, size = round(E[it,]), prob = death_prob )
      A[it + 1, ] <- A[it,] +  new_infA - new_rec_A - rbinom(n = Ncomp, size = round(A[it,]), prob = death_prob)
      I[it + 1, ] <- I[it,] +  new_infI - new_rec_I - rbinom(n = Ncomp, size = round(I[it,]), prob = death_prob)
      R[it + 1, ] <- R[it,] +  new_rec_I + new_rec_A - rbinom(n = Ncomp, size = round(R[it,]), prob = death_prob)
      
      ## make incidence the new number of daily individuals becoming infected
      incid_A[it, ] <- new_infA
      incid_I[it, ] <- new_infI
    }
    
    ## deterministic equations to check -- does not currently work 
    if(stoch == F){
      S[it + 1, ] <- S[it,] + delta.t * (births - seas[it] * WI * S[it,] * dw / N  - deaths*S[it,])
      E[it + 1, ] <- E[it,] + delta.t * (seas[it] * WI * S[it,] * dw / N - deaths*E[it,] - sigma*E[it,])
      A[it + 1, ] <- A[it,] + delta.t * ( (1 - prop_symptomatic)*sigma*E[it,]   - A[it,]*(gamma - deaths))
      I[it + 1, ] <- I[it,] + delta.t * (  prop_symptomatic*sigma*E[it,] - I[it,]*(gamma - deaths) )
      R[it + 1, ] <- R[it,] + delta.t * (A[it,]*gamma+ I[it,]*gamma - R[it,]* deaths)
      incid_A[it,] <-  (1 - prop_symptomatic)*(seas[it] * WI * S[it,] * dw / N)
      incid_I[it,] <- prop_symptomatic*(seas[it] * WI * S[it,] * dw / N)
    }
  }
  
  out <- data.frame(cbind(time,S,E,A,I,R,incid_A, incid_I,R0, Re))
  names(out) <- c('time',names(ICs),"R0", "Reff")
  ## output is the number in each class per time point per age-category+homeless+healthcare workers
  return(out)
  
}

## ---- function to set up and organize the mixing data, inital conditions, parameters, etc. ---- ####

setup_seir_model <- function(stoch, R0, c_scale_vec,
                             gamma=1/6.5, sigma=1/5.2,
                             phase=0, beta1=0, mu=0, v=0,
                             prop_symptomatic, sd.dw=0.05){
  #data <- read.csv('Zambia_AgePopulation.csv')
  data <- read.csv('Lusaka_AgePopulation.csv')
  age_data <- make_age_structure_matrix(data)
  Zam_pop = age_data$Zam_pop
  Ncomp = age_data$Ncomp
  W <- make_polymod_matrix()
  rescale_mixing <- rescale_age_matrix(Ncomp, W, Zam_pop, c_scale_vec)
  W <- rescale_mixing$W
  C <- rescale_mixing$C
  
  ## set prop_symtomatic
  prop_symptomatic <- prop_symptomatic
  if(length(prop_symptomatic)!=Ncomp){warning("prop_symptomatic is the wrong length")}
  
  ## set initial conditions
  ICs <- c(S = Zam_pop, 
           E = rep(0,length(Zam_pop)),
           A = rep(8,length(Zam_pop)),
           I = rep(2,length(Zam_pop)), 
           R = Zam_pop,
           incid_A = rep(0,Ncomp),
           incid_I = rep(0,Ncomp))
  
  ## set the R compartment
  ## crudely just set anything negative to be zero but we can fine tune this more depending on what we assume S0 does
  ICs[(4*Ncomp+1):(5*Ncomp)] <- Zam_pop -  ICs[1:Ncomp] - ICs[(Ncomp+1):(2*Ncomp)] - ICs[(2*Ncomp+1):(3*Ncomp)] - ICs[(3*Ncomp+1):(4*Ncomp)]
  ICs[(4*Ncomp+1):(5*Ncomp)] [ICs[(4*Ncomp+1):(5*Ncomp)]  < 0 ] <- 0
  
  ## population sizes by demographic data
  N <- Zam_pop
  
  ## units in days!! 
  gamma <- gamma ## infectious period
  sigma <- sigma ## latent period
  phase <- phase ## when should seasonal forcing peak?
  mu <- mu ## set births to be zero currently
  v <- v ## set natural death rate to be zero currently
  #R0 <- 2.5 ## make a range?   ## R0 = beta * sigma / ((sigma + v) * (v + gamma)) for SEAIR model with split proportion into A-I and only A and I contributing to infection
  beta0 <- R0 * (gamma + v) * (sigma + v) / sigma ## set beta based on that value
  beta1 <- beta1 ## seasonal forcing should be modest here
  sd.dw <- sd.dw
  
  ## now check to make sure the R0 we get is the R0 we put in
  ## using the same formula as above
  R0.mat <- matrix(0,Ncomp,Ncomp)
  for (i in 1:Ncomp){
    for (j in 1:Ncomp){
      R0.mat[i,j] <- W[i,j]*Zam_pop[i]/Zam_pop[j]* beta0 * sigma / ( (sigma + v) * (v + gamma))
    }
  }
  # for (i in 1:Ncomp){
  #   for (j in 1:Ncomp){
  #     R0.mat[i,j] <- W[i,j]*Zam_pop[i]/Zam_pop[j]* beta0 * (gamma + v) * (sigma +v )/sigma
  #   }
  # }
  print(eigen(R0.mat)$values[1]) ## just a check
  return(list(C = C, W = W, beta0 = beta0, beta1 = beta1, 
              phase = phase, mu = mu, v = v, ICs = ICs, 
              Ncomp = Ncomp, N=N, gamma=gamma, sigma = sigma,
              prop_symptomatic=prop_symptomatic, sd.dw=sd.dw))
}
