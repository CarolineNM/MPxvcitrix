#rm(list=ls())
library(R2jags)
library(runjags)
library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)
options(scipen=999)

###Difference between this and fixed inits 1 is the phi prior

textstring <- "
model {
  
  beta ~ dnorm(1.3, 0.05) T(0,)  # Truncated to prevent very high beta values
 
  kappa ~ dbeta(7, 5)     # Prior for kappa (mixing matrix parameter).Assuming higher assortative mixing
  
  #phi_tran ~ dnorm(0, 0.05) T(0,)  # Half-normal prior for 1/sqrt(phi)
  #phi <- 1/sqrt(phi_tran)         ##moderate prior that does not favour high overdispersion
  
  
  #delta_inv<-6.26 ## duration in E comp 6.26(5.55-6.97)
  delta_inv~dgamma(300.65298,47.98256) ## duration in E comp 6.26(5.55-6.97)
  delta <-1/delta_inv   ##daily rate of transition from E to P
  
  epsilon_inv<-2     ##duration in P comp 2days
  epsilon <- 1/epsilon_inv ##daily rate of transition from P to A and I
  
  #m<-0.065    ##proportion of Asymptomatic individuals 6.5%(3.84%-10.8%)
  m~dbeta(14.3,201.01)    ##proportion of Asymptomatic individuals 6.5%(3.84%-10.8%)
  
  #theta_invall<-21 ##duration in all I comps 21(7-28)
  theta_invall~dgamma(42.3,1.99) ##duration in all I comps 21(7-28)
  theta_inv<-theta_invall/3    ##duration in each I comp
  theta <- 1/theta_inv    ##daily rate of transition I1-I2-I3-R
  
  #omega_invall<-14     ###duration in all A comps 14(7-21)days
  omega_invall~ dgamma(20.6,1.45)     ###duration in all A comps 14(7-21)days
  omega_inv<-omega_invall/3     ###duration in each A comp
  omega<-1/omega_inv     ##daily rate of transition A1-A2-A3-R

  
   # Reporting fraction and lag
   
  report_frac~dbeta(7.62,2.37)  # 78%(48-96) reporting fraction dist.
  #report_frac <- 0.78  # 78%(48-96) reporting fraction fixed
  #report_frac <- 1  # 78%(48-96) reporting fraction fixed
  report_lag <- 7     # 7-day reporting lag
  
  
  # Initial conditions for each group
  for (g in 1:3) {
    S[g, 1] <- N[g] - E_start[g] - P_start[g] - A1_start[g] - A2_start[g] - A3_start[g] - I1_start[g] - I2_start[g] - I3_start[g] - R_start[g]
    E[g, 1] <- E_start[g]
    P[g, 1] <- P_start[g]
    A1[g, 1] <- A1_start[g]
    A2[g, 1] <- A2_start[g]
    A3[g, 1] <- A3_start[g]
    I1[g, 1] <- I1_start[g]
    I2[g, 1] <- I2_start[g]
    I3[g, 1] <- I3_start[g]
    R[g, 1] <- R_start[g]
    C[g, 1]<-0     ##initialize reported cases comp
  }

  # Define the assortative mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_ass[g, g_] <- ifelse(g == g_, 1, 0)
    }
  }
  
  # Define the proportional mixing matrix once outside of the loop
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix_prp[g, g_] <- (contact[g_] * N[g_]) / sum(contact * N)
    }
  }
  
  # Initialize prevalence for t = 1
  for (g in 1:3) {
    prevalence[g, 1] <- (I1[g, 1] + I2[g, 1] + I3[g, 1]) / N[g]
  }

  # Initialize mixing matrix (mix) for t = 1
  for (g in 1:3) {
    for (g_ in 1:3) {
      mix[g, g_, 1] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
    }
  }

 
 # Calculate initial lambda_det for t = 1 based on the initial conditions
for (g in 1:3) {
  lambda_det[g, 1] <- sum(beta * contact[g] * mix[g, 1:3, 1] * prevalence[1:3, 1])
}

 ####calculating new cases to include in likelihood


    total_new_cases[1] <- 0 # Initial people in Icomp
    

  # Dynamics of the model (deterministic-like)
  for (t in 2:T) {
  
  # Final mixing matrix for three groups only
    for (g in 1:3) {
      for (g_ in 1:3) {
        mix[g, g_, t] <- kappa * mix_ass[g, g_] + (1 - kappa) * mix_prp[g, g_]
      }
    }

    # Calculate prevalence once per time step (deterministic)
    for (g in 1:3) {
      prevalence[g, t] <- (I1[g, t-1] + I2[g, t-1] + I3[g, t-1]) / N[g]
    }

    # Calculate force of infection (lambda_det) for each group
    for (g in 1:3) {
      lambda_det[g, t] <- sum(beta * contact[g] * mix[g, 1:3, t] * prevalence[1:3, t])
    }

  
    # Deterministic transitions instead of binomial draws
    for (g in 1:3) {
      new_Es[g, t-1] <- lambda_det[g, t] * S[g, t-1]     # Deterministic S to E transition
      new_Ps[g, t-1] <- delta * E[g, t-1]             # Deterministic E to P transition
      new_A1s[g, t-1] <- epsilon * m * P[g, t-1]      # Deterministic P to A1
      new_A2s[g, t-1] <- omega * A1[g, t-1]           # Deterministic A1 to A2
      new_A3s[g, t-1] <- omega * A2[g, t-1]           # Deterministic A2 to A3
      new_RAs[g, t-1] <- omega * A3[g, t-1]           # Deterministic A3 to R
      new_I1s[g, t-1] <- epsilon * (1 - m) * P[g, t-1] # Deterministic P to I1
      new_I2s[g, t-1] <- theta * I1[g, t-1]           # Deterministic I1 to I2
      new_I3s[g, t-1] <- theta * I2[g, t-1]           # Deterministic I2 to I3
      new_RIs[g, t-1] <- theta * I3[g, t-1]           # Deterministic I3 to R
      

        
      # Update compartments deterministically at each time step
      S[g, t] <- max(0, S[g, t-1] - new_Es[g, t-1])
      E[g, t] <- max(0, E[g, t-1] + new_Es[g, t-1] - new_Ps[g, t-1])
      P[g, t] <- max(0, P[g, t-1] + new_Ps[g, t-1] - new_A1s[g, t-1] - new_I1s[g, t-1])
      A1[g, t] <- max(0, A1[g, t-1] + new_A1s[g, t-1] - new_A2s[g, t-1])
      A2[g, t] <- max(0, A2[g, t-1] + new_A2s[g, t-1] - new_A3s[g, t-1])
      A3[g, t] <- max(0, A3[g, t-1] + new_A3s[g, t-1] - new_RAs[g, t-1])
      I1[g, t] <- max(0, I1[g, t-1] + new_I1s[g, t-1] - new_I2s[g, t-1])
      I2[g, t] <- max(0, I2[g, t-1] + new_I2s[g, t-1] - new_I3s[g, t-1])
      I3[g, t] <- max(0, I3[g, t-1] + new_I3s[g, t-1] - new_RIs[g, t-1])
      R[g, t] <- max(0, R[g, t-1] + new_RAs[g, t-1] + new_RIs[g, t-1])
      
      
         # Calculate reported cases based on all infectious compartments with reporting lag
        # Calculate reported cases based on new_I1s only
          new_reported_cases[g, t] <- step(t - report_lag) * report_frac * 
                            new_I1s[g, max(1, t - report_lag)]

    }
    
       # Total new reported cases at time t
        total_new_cases[t] <- sum(new_reported_cases[1:3, t])
  }
  

   for (t in (burn_in_timesteps + 1):(T_obs + burn_in_timesteps)) {
     cases_obs[t - burn_in_timesteps] ~ dpois(total_new_cases[t] + 1e-6) ##poisson likelihood
     # Posterior predictive distribution
    cases_pred[t - burn_in_timesteps] ~ dpois(total_new_cases[t] + 1e-6)
   }
}"

dat=data.frame(setting="Vancouver",
               Sex_groups=c("0-4","4-22","22-100"),
               mean_rate=c(1.31,9.13,39.84),
               mean_rate_dy=c(1.31,9.13,39.84)/180,
               Prop=c(0.60,0.35,0.05),
               Pop=c(0.60,0.35,0.05)*26100)


contact <- dat$mean_rate / 180 ##convert the 6months sexual contact rates to daily contact rates
#contact <- dat$mean_rate / 26 ##convert the 6months sexual contact rates to weekly contact rates
N <- c(15660, 9135, 1305) ##population sizes for the different groups
burn_in_timesteps <- 30
T=169+burn_in_timesteps
##for each grp
##these initial conditions were the output of Negb21(Mainmod_negbin3)
E_start=c(8,3,14)
P_start=c(2,0,5)
A1_start=c(1,0,3)
A2_start=c(0,0,0)
A3_start=c(0,0,0)
I1_start=c(3,0,6)
I2_start=c(0,0,0)
I3_start=c(0,0,0)
R_start=c(0,0,0)


#I1_start+I2_start+I3_start

case_dat=read_excel("data/case_data_V2.xlsx",sheet="cases")
cases_obsb = case_dat$total_cases 
T_obsb=length(cases_obsb)

# Data list for JAGS
dataList <- list(N = N, T = T, contact = contact, E_start = E_start, 
                 P_start = P_start,A1_start=A1_start,
                 A2_start = A2_start,A3_start=A3_start,
                 I1_start=I1_start,I2_start = I2_start,
                 I3_start=I3_start, R_start = R_start,
                 cases_obs=cases_obsb,
                 T_obs=T_obsb,
                 burn_in_timesteps=burn_in_timesteps)


inits_list <- list(
  list(beta = 1.39, kappa = 0.9),  # Initial values for chain 1
  list(beta = 2.01, kappa = 0.83)   # Initial values for chain 2
)

#Run the model with different initial values for each chain
out_negb3<- run.jags(textstring, data = dataList,
                     monitor = c("beta", "kappa","cases_pred"),
                     sample = 60000, adapt = 10000, burnin = 20000, thin = 20,
                     n.chains = 2, inits = inits_list,
                     summarise = FALSE)

dailyout_negb3<- as.mcmc.list(out_negb3)
#current_outputb<-current_output[[1]] 
save(dailyout_negb3,file="data_25/dailyout_negb3.RData") 