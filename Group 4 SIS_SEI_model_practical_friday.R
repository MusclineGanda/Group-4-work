##############################
## SIS–SEI vector–host modelling ##
##############################
library(deSolve)
library(tidyverse)
library(plotly)


# -------------------------
# Parameters 
# -------------------------
biting = 0.3
seeking = 0.1
infectivity = 0.45
alpha = biting * seeking * infectivity
k = 8 # mosquito to human ratio
l=1/7 #latent period for mosquitoes

params <- list(
  alpha     = alpha,      # transmission parameter from mosquitoes to humans
  beta = 0.5,      # P(human infection per bite from  an infectious mosquito)
  gamma = 1/25,    # human recovery rate
  mu_d  = 1/14, # mosquito death rate (1/lifespan)
  mu_b = 1/14, # based on our assumption the birth rate= death rate
  Nh = 15000,     # total number of humans and we assume that the system is closed
  Nm = 15000*k,       # mosquito:human ratio  (Nm = k * Nh)
  l= 1/7#latent period for mosquitoes
)


# -------------------------
# Model state variables
# Sh, Ih : susceptible & infectious humans
# Sm, Em, Im : susceptible,Exposed & infectious mosquitoes
# -------------------------

state0 <- c(
  Sh = params$Nh - 25,  # 25 people are infected
  Ih = 25,
  Sm = params$Nm - 0.12*params$Nm,  # 12% infectious mosquitoes
  Em = 0,
  Im = 0.12*params$Nm
  
  
)

# -------------------------
# ODE system
# -------------------------
sis_sei <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih # should remain constant 
    Nm <- Sm + Em + Im  # should remain constant 
    
    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + gamma * Ih
    dIh <-  (alpha * (Im / Nh)   * Sh) - gamma * Ih
    
    # Mosquitoes
    births<-mu_b*Nm
    dSm <- births - (beta * (Ih / Nh)* Sm) - mu_d* Sm
    dEm <- (beta * (Ih / Nh)* Sm)-(l*Em)-(mu_d*Em)
    dIm <-  (l* Em) - mu_d* Im
    
    newInfections <- alpha * (Im/Nh) * Sh
    
    list(c(dSh, dIh, dSm,dEm, dIm), aux = c(newInfections = newInfections,Nh, Nm))
  })
}



# -------------------------
# Question 2: Simulating the ODE at 365 without intervention
# -------------------------
times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out <- as.data.frame(ode(y = state0, times = times, func = sis_sei, parms = params))

cum_infections_baseline <- tail(out$Ih, 1)

View(out)
### plotting

subplot(
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human Population"),
  
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Em, name = "E_m") %>%  
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito Population"),
  
  nrows = 1
)
# 2a.Calculate the prevalence of the disease in the population at the end of the simulation period?
# Overall number of cases/Total population *100
# Calculate prevalence at the end of the simulation (after 365 days)

final_infected_humans <- out$Ih[366]  # Number of infected humans at the end
prevalence<-(final_infected_humans/15000)*100
cat("Prevalence after 365 days:", prevalence, "%\n")

####Ans Prevalence after 365 days: 59.18367 % #####

#2b.What is your interpretation of the computed prevalence?



#The prevalence at the end of the simulation period is  59.18%  meaning that over half of the human population was infected without malaria control interventions 


# Question 3: Update the model structure: We now consider the following intervention scenarios in the model, i.	Scenario 1: Vector control — Increases the mortality rate of mosquitoes by 60%. Assume the vector control program has an annual budget of 80,000 USD to implement the intervention.

biting = 0.3
seeking = 0.1
infectivity = 0.45
alpha = biting * seeking * infectivity
k = 8 # mosquito to human ratio
l=1/7 #latent period for mosquitoes

params_vector_ctrl <- list(
  alpha     = alpha,      # transmission parameter from mosquitoes to humans
  beta = 0.5,      # P(human infection per bite from  an infectious mosquito)
  gamma = 1/25,    # human recovery rate
  mu_d  = 1/14+(0.6*(1/14)), # increase in mosquito death rate by 60%
  mu_b = 1/14, # based on our assumption the birth rate= death rate
  Nh = 15000,     # total number of humans and we assume that the system is closed
  Nm = 15000*k,       # mosquito:human ratio  (Nm = k * Nh)
  l= 1/7#latent period for mosquitoes
)


# -------------------------
# Model state variables
# Sh, Ih : susceptible & infectious humans
# Sm, Em, Im : susceptible,Exposed & infectious mosquitoes
# -------------------------

statevctrl <- c(
  Sh = params_vector_ctrl$Nh - 25,  # 25 people are infected
  Ih = 25,
  Sm = params_vector_ctrl$Nm - 0.12*params_vector_ctrl$Nm,  # 12% infectious mosquitoes
  Em = 0,
  Im = 0.12*params_vector_ctrl$Nm
  
  
)

# -------------------------
# ODE system
# -------------------------
sis_sei_vector_ctrl <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih # should remain constant 
    Nm <- Sm + Em + Im  # should remain constant 
    
    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + gamma * Ih
    dIh <-  (alpha * (Im / Nh)   * Sh) - gamma * Ih
    
    # Mosquitoes
    births<-mu_b*Nm
    dSm <- births - (beta * (Ih / Nh)* Sm) - mu_d* Sm
    dEm <- (beta * (Ih / Nh)* Sm)-(l*Em)-(mu_d*Em)
    dIm <-  (l* Em) - mu_d* Im
    
    newInfections <- alpha * (Im/Nh) * Sh
    
    list(c(dSh, dIh, dSm,dEm, dIm), aux = c(newInfections = newInfections,Nh, Nm) )
  })
}

# -------------------------
times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out_with_vector_ctrl <- as.data.frame(ode(y = statevctrl, times = times, func = sis_sei_vector_ctrl, parms = params_vector_ctrl))

View(out_with_vector_ctrl)

cum_infections_vector <- tail(out_with_vector_ctrl$Ih,1)
infections_averted_mosq_vec <- cum_infections_baseline - cum_infections_vector

# -----------------------------
# Plotting cases against baseline
# -----------------------------
cases_after365_df <- data.frame(
  Intervention = c("No Intervention", "After Vector Control"),
  Cases = c(cum_infections_baseline, cum_infections_vector)
)

ggplot(cases_after365_df, aes(x=Intervention, y=Cases, fill=Intervention)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Cum_Cases After 365 Days",
       y="Number of Cum Cases",
       x="Intervention")


### plotting

subplot(
  plot_ly(out_with_vector_ctrl, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human Population"),
  
  plot_ly(out_with_vector_ctrl, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Em, name = "E_m") %>%  
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito Population"),
  
  nrows = 1
)


final_infected_humans_after_vectorctrl <- out_with_vector_ctrl$Ih[366]  # Number of infected humans at the end
prevalence_after_vetor_ctrl<-(final_infected_humans_after_vectorctrl/15000)*100
cat("Prevalence after 365 days:", prevalence_after_vetor_ctrl, "%\n")

#Number of cases averted by vector control
infections_averted_mosq_vec = out$Ih[366] - out_with_vector_ctrl$Ih[366]
cat("Number of cases averted after vector control:",infections_averted_mosq_vec)

#The prevalence at the end of the simulation period decreased significantly after introducing vector control where death rate for mosquitoes was increased by 60%, approximately 8877 cases were averted by vector control hence worth the investment 


#Scenario 2: Introducing treatment and reducing infectious period by 40% through treatment with 0.9 efficacy
biting = 0.3
seeking = 0.1
infectivity = 0.45
alpha = biting * seeking * infectivity
k = 8 # mosquito to human ratio
l=1/7 #latent period for mosquitoes

params_treatment <- list(
  alpha     = alpha,      # transmission parameter from mosquitoes to humans
  beta = 0.5,      # P(human infection per bite from  an infectious mosquito)
  gamma_untreated = 1/25, # human recovery rate without treatment
  gamma_treatment = (0.9*1/15) + (0.1*1/25), #recovery rate after treatment increased by 40% at 0.9 treatment eff
  mu_d  = 1/14, # death rate
  mu_b = 1/14, # based on our assumption the birth rate= death rate
  Nh = 15000,     # total number of humans and we assume that the system is closed
  Nm = 15000*k,       # mosquito:human ratio  (Nm = k * Nh)
  l= 1/7,#latent period for mosquitoes
  efi=0.9, #efficacy of treatment
  ratio_access_T=0.75
  
)

# -------------------------
# Model state variables
# Sh, Ih_treated , Ih_untreated: susceptible, infectious humans but treated, infectious humans but not treated humans
# Sm, Em, Im : susceptible,Exposed & infectious mosquitoes
# -------------------------

state_treatment <- c(
  Sh = params_treatment$Nh - 25,  # 25 people are infected
  Ihtreated=0.75*25,
  Ihuntreated=0.25*25,
  Sm = params_treatment$Nm - 0.12*params_treatment$Nm,  # 12% infectious mosquitoes
  Em = 0,
  Im = 0.12*params_treatment$Nm
  
  
)

# -------------------------
# ODE system
# -------------------------
sis_sei_treatment <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ihtreated + Ihuntreated # should remain constant 
    Nm <- Sm + Em + Im  # should remain constant 
    
    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + (params_treatment$gamma_treatment *  Ihtreated)+(params_treatment$gamma_untreated*Ihuntreated)
    dIh_treated <-  0.75*(alpha * (Im / Nh)   * Sh) - ((params_treatment$gamma_treatment *  Ihtreated))
    dIh_untreated <- 0.25*(alpha * (Im / Nh)   * Sh) - ((params_treatment$gamma_untreated *  Ihuntreated))
    
    # Mosquitoes
    births<-mu_b*Nm
    dSm <- births - (beta * (Ihtreated / Nh)* Sm) - (beta * (Ihuntreated / Nh)* Sm)- mu_d* Sm
    dEm <- (beta * (Ihtreated / Nh)* Sm) + (beta * (Ihuntreated / Nh)* Sm)-(l*Em)-(mu_d*Em)
    dIm <-  (l* Em) - mu_d* Im
    
    newInfections <- 0.75*(alpha * (Im / Nh)   * Sh) + 0.25*(alpha * (Im / Nh)   * Sh)
    
    list(c(dSh, dIh_treated, dIh_untreated, dSm,dEm, dIm),  aux = c(newInfections = newInfections, Nh, Nm))
  })
}



# -------------------------
times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out_with_treatment<- as.data.frame(ode(y = state_treatment, times = times, func = sis_sei_treatment, parms = params_treatment))

view(out_with_treatment)

cum_infections_treatment <- tail(out_with_treatment$Ihtreated+out_with_treatment$Ihuntreated, 1)
infections_averted_bytreatment <- cum_infections_baseline - cum_infections_treatment
cat("Number of cases averted after Treatment:",infections_averted_bytreatment)

# Number of cases averted after Treatment: 1351.276

# -----------------------------
# Plotting cases against baseline
# -----------------------------
cases_after365_df <- data.frame(
  Intervention = c("No Intervention", "After Treatment Intervention"),
  Cases = c(cum_infections_baseline, cum_infections_treatment)
)

ggplot(cases_after365_df, aes(x=Intervention, y=Cases, fill=Intervention)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Cum_Cases After 365 Days",
       y="Number of Cum Cases",
       x="Intervention")


### plotting

subplot(
  plot_ly(out_with_treatment, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_humans") %>%
    add_lines(y = ~Ihtreated, name = "I_htreated") %>%
    add_lines(y = ~Ihuntreated, name = "I_treated") %>%
    layout(title = "Human Population"),
  
  plot_ly(out_with_treatment, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Em, name = "E_m") %>%  
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito Population"),
  
  nrows = 1
)


final_infected_humans_with_treatment <- out_with_treatment$Ihtreated[366]+out_with_treatment$Ihuntreated[366]  # Number of infected humans at the end after treatment
prevalence_after_treatment<-(final_infected_humans_with_treatment/15000)*100
cat("Prevalence after 365 days:", prevalence_after_treatment, "%\n")
#Prevalence after 365 days: 50.17516 %


#Cost Analysis and more value for money assessment

cat("Number of cases averted after treatment:",infections_averted_bytreatment)
cat("Number of cases averted after vector control:",infections_averted_mosq_vec)

##  Number of cases averted after treatment: 1351.276       ##########
##  Number of cases averted after vector control: 8877.534   #########


Total_cost_of_treatment<- sum(out_with_treatment$aux.newInfections)*0.75*15
cat("Total cost of treatment intervention:",Total_cost_of_treatment)

## Total cost of treatment intervention: 1670117
## we then conclude that vector control is cheaper based on the assumption that treated individuals do not fully recover hence more expensive that USD 80 000 for vector control which avert almost all cases. we also noticed that vector control is more effective when looking at malaria control , reducing the vector has a significant effect on the transmission of malaria to human beings through the reduction in the force of transmission.

sum(out_with_treatment$aux.newInfections)
view(out_with_treatment)
# -----------------------------
# Compare and plot cases averted
# -----------------------------
cases_averted_df <- data.frame(
  Intervention = c("Vector Control", "Treatment"),
  CasesAverted = c(infections_averted_mosq_vec, infections_averted_bytreatment)
)

ggplot(cases_averted_df, aes(x=Intervention, y=CasesAverted, fill=Intervention)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Cases Averted By Intervention",
       y="Number of Cases Averted",
       x="Intervention")


#Scenario 3, combined Interventions effect

biting = 0.3
seeking = 0.1
infectivity = 0.45
alpha = biting * seeking * infectivity
k = 8 # mosquito to human ratio
l=1/7 #latent period for mosquitoes

# Parameters
params_treatment_vector <- list(
  alpha = alpha,                        # Transmission parameter from mosquitoes to humans
  beta = 0.5,                           # P(human infection per bite from an infectious mosquito)
  gamma = 1/25,                         # Human recovery rate without treatment
  gamma_treatment = 1/25 * (1 - 0.4),   # Recovery rate after treatment increased by 40%
  mu_d = 1/14 * (1 + 0.6),              # Death rate of mosquitoes increased by 60%
  mu_b = 1/14,                          # Birth rate of mosquitoes
  Nh = 15000,                           # Total number of humans (system is closed)
  Nm = 15000 * k,                       # Mosquito:human ratio (Nm = k * Nh)
  l = 1/7,                              # Latent period for mosquitoes
  efi = 0.9,                            # Efficacy of treatment
  ratio_access_T = 0.75                 # Proportion of infectious population accessing treatment
)

# -------------------------
# Model state variables
# Sh, Ih , Th: susceptible, infectious humans, Treated humans
# Sm, Em, Im : susceptible,Exposed & infectious mosquitoes
# -------------------------
state_treatment_vector <- c(
  Sh = params_treatment$Nh - 25,  # 25 people are infected
  Ihtreated=0.75*25,
  Ihuntreated=0.25*25,
  Sm = params_treatment$Nm - 0.12*params_treatment$Nm,  # 12% infectious mosquitoes
  Em = 0,
  Im = 0.12*params_treatment$Nm
  
)

# -------------------------
# ODE system
# -------------------------
sis_sei_treatment_vector <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ihtreated + Ihuntreated # should remain constant 
    Nm <- Sm + Em + Im  # should remain constant 
    
    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + (params_treatment$gamma_treatment *  Ihtreated)+(params_treatment$gamma_untreated*Ihuntreated)
    dIh_treated <-  0.75*(alpha * (Im / Nh)   * Sh) - ((params_treatment$gamma_treatment *  Ihtreated))
    dIh_untreated <- 0.25*(alpha * (Im / Nh)   * Sh) - ((params_treatment$gamma_untreated *  Ihuntreated))
    
    # Mosquitoes
    births<-mu_b*Nm
    dSm <- births - (beta * (Ihtreated / Nh)* Sm) - (beta * (Ihuntreated / Nh)* Sm)- mu_d* Sm
    dEm <- (beta * (Ihtreated / Nh)* Sm) + (beta * (Ihuntreated / Nh)* Sm)-(l*Em)-(mu_d*Em)
    dIm <-  (l* Em) - mu_d* Im
    
    newInfections <- alpha * (Im/Nh) * Sh
    
    list(c(dSh, dIh_treated, dIh_untreated, dSm,dEm, dIm),  aux = c(newInfections = newInfections, Nh, Nm))
  })
}

times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out_with_treatment_vector<- as.data.frame(ode(y = state_treatment_vector, times = times, func = sis_sei_treatment_vector, parms = params_treatment_vector))

view(out_with_treatment_vector)

### plotting

subplot(
  plot_ly(out_with_treatment_vector, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_humans") %>%
    add_lines(y = ~Ihtreated, name = "Infected Treated") %>%
    add_lines(y = ~Ihuntreated, name = "Infected Un Treatment") %>%
    layout(title = "Human Population"),
  
  plot_ly(out_with_treatment_vector, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Em, name = "E_m") %>%  
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito Population"),
  
  nrows = 1
)

View(out_with_treatment_vector)


final_infected_humans_with_both <- out_with_treatment_vector$Ihtreated[366]+out_with_treatment_vector$Ihuntreated[366]  # Number of infected humans at the end after treatment
prevalence_after_both<-(final_infected_humans_with_both/15000)*100
cat("Prevalence after 365 days:", prevalence_after_both, "%\n")

#Number of cases averted by treatment
infections_averted_byboth = out$Ih[366] - (out_with_treatment_vector$Ihtreated[366]+out_with_treatment_vector$Ihuntreated[366])
cat("Number of cases averted after both vector control and treatment:",infections_averted_byboth)
##Number of cases averted after both vector control and treatment: 8877.548

#Policy
# since almost all cases were averted by vector control ,there is therefore no need for treatment if we are to look at the cost, because the same cases are averted using both interventions. 

