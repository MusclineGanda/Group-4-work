#########Code##################


library(deSolve)
library(tidyverse)
library(plotly)

# Parameters
biting = 0.3
seeking = 0.1
infectivity = 0.45
alpha = biting * seeking * infectivity
k = 8 # mosquito to human ratio
params <- list(
  alpha = alpha, # transmission parameter from mosquitoes to humans
  beta = 0.5, # P(human infection per bite from an infectious mosquito)
  gamma = 1/25, # human recovery rate
  mu_m = 1/10, # mosquito death rate (1/lifespan)
  Nh = 15000, # total number of humans
  Nm = 15000*k # mosquito to human ratio (Nm = k * Nh)
)

# Model state variables
state0 <- c(
  Sh = params$Nh - 25, # 0.167% initially infectious
  Ih = 25,
  Sm = params$Nm - 0.12 * params$Nm, # 12% infectious mosquitoes
  Im = 0.12 * params$Nm
)

# ODE system
sis_si <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih # should remain constant
    Nm <- Sm + Im # should remain constant
    
    # Humans
    dSh <- -(alpha * (Im / Nh) * Sh) + gamma * Ih
    dIh <- (alpha * (Im / Nh) * Sh) - gamma * Ih
    
    # Mosquitoes
    births <- mu_m * Nm
    dSm <- births - (beta * (Ih / Nh) * Sm) - mu_m * Sm
    dIm <- (beta * (Ih / Nh) * Sm) - mu_m * Im
    
    list(c(dSh, dIh, dSm, dIm))
  })
}

# Solve
times <- seq(0, 365, by = 1) # 1 year and daily time steps
out <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params))

# Calculate prevalence
####  2a. Calculate the prevalence of the disease in the population at the end of the simulation period?


# Calculate point prevalence

point_prevalence <- out$Ih[nrow(out)] / params$Nh
print(point_prevalence)



# Calculate period prevalence (assuming the entire simulation period)
period_prevalence <- mean(out$Ih / params$Nh)
print(period_prevalence)


### 2b. What is your interpretation of the computed prevalence?

if (prevalence < 0.1){
  message ("The Prevalence is very low in the Population")
  
} else if (prevalence < 0.5){
  message ("There is a moderate Prevalence in the Population")
  
} else if (prevalence > 0.5 ){
  message ( "There is High Prevalence in the Population, disease is Highly endemic and Widespread" )
}



# Plot
subplot(
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human population"),
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Human & Mosquito population"),
  nrows = 1
)



# Intervention scenarios
# Scenario 1: Vector control
params_vector_control <- params
params_vector_control$alpha <- params_vector_control$alpha * (1-0.6) # 60% reduction in host seeking rate

out_vector_control <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params_vector_control))


# Plot Scenario 1: Vector Control
subplot(
  plot_ly(out_vector_control, x = ~time) %>% 
    add_lines(y = ~Sh, name = "S_h") %>% 
    add_lines(y = ~Ih, name = "I_h") %>% 
    layout(title = "Human population (Vector Control)"),
  plot_ly(out_vector_control, x = ~time) %>% 
    add_lines(y = ~Sm, name = "S_m") %>% 
    add_lines(y = ~Im, name = "I_m") %>% 
    layout(title = "Mosquito population (Vector Control)"),
  nrows = 1
)





# Scenario 2: Vaccination
sis_si_vaccination <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih + Vh # should remain constant
    Nm <- Sm + Im # should remain constant
    
    # Humans
    dSh <- -(alpha * (Im / Nh) * Sh) + (gamma * Ih) + (1/30 * Vh) - (0.3 * 0.65*Sh)
    dIh <- (alpha * (Im / Nh) * Sh) - gamma * Ih
    dVh <- (0.65*0.3*Sh) - (1/30 * Vh)
    
    # Mosquitoes
    births <- mu_m * Nm
    dSm <- births - (beta * (Ih / Nh) * Sm) - mu_m * Sm
    dIm <- (beta * (Ih / Nh) * Sm) - mu_m * Im
    
    list(c(dSh, dIh, dVh, dSm, dIm))
  })
}
state0_vaccination <- c(
  Sh = (params$Nh - 25) * (1-0.3),
  Ih = 25,
  Vh = (params$Nh -25) * 0.3,
  Sm = params$Nm - 0.12 * params$Nm,
  Im = 0.12 * params$Nm
)

print(state0_vaccination)


params_vaccination <- params
out_vaccination <- as.data.frame(ode(y = state0_vaccination, times = times, func = sis_si_vaccination, parms = params_vaccination))

# Plot Scenario 2: Vaccination
subplot(
  plot_ly(out_vaccination, x = ~time) %>% 
    add_lines(y = ~Sh, name = "S_h") %>% 
    add_lines(y = ~Ih, name = "I_h") %>% 
    add_lines(y = ~Vh, name = "V_h") %>% 
    layout(title = "Human population (Vaccination)"),
  plot_ly(out_vaccination, x = ~time) %>% 
    add_lines(y = ~Sm, name = "S_m") %>% 
    add_lines(y = ~Im, name = "I_m") %>% 
    layout(title = "Mosquito population (Vaccination)"),
  nrows = 1
)

View(out_vaccination)


# Scenario 3: Combined Intervention

# Combined intervention (Vector control + Vaccination)
sis_si_combined <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih + Vh # should remain constant
    Nm <- Sm + Im # should remain constant
    
    # Humans
    dSh <- -(alpha * (Im / Nh) * Sh) + (gamma * Ih) + (1/30 * Vh) - (0.3 * 0.65*Sh)
    dIh <- (alpha * (Im / Nh) * Sh) - gamma * Ih
    dVh <- (0.65*0.3*Sh) - (1/30 * Vh)
    
    # Mosquitoes
    births <- mu_m * Nm
    dSm <- births - (beta * (Ih / Nh) * Sm) - mu_m * Sm
    dIm <- (beta * (Ih / Nh) * Sm) - mu_m * Im
    
    list(c(dSh, dIh, dVh, dSm, dIm))
  })
}

state0_combined <- c(
  Sh = (params$Nh - 25) - ((params$Nh - 25)* (0.3)),
  Ih = 25,
  Vh = (params$Nh -25) * 0.3,
  Sm = params$Nm - 0.12 * params$Nm,
  Im = 0.12 * params$Nm
)

print(state0_combined)


params_combined <- params
params_combined$alpha <- params_combined$alpha * (1-0.6) # 60% reduction in host seeking rate

out_combined <- as.data.frame(ode(y = state0_combined, times = times, func = sis_si_combined, parms = params_combined))

# Plot Scenario 3: Combined intervention
subplot(
  plot_ly(out_combined, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    add_lines(y = ~Vh, name = "V_h") %>%
    layout(title = "Human population (Combined)"),
  plot_ly(out_combined, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito population (Combined)"),
  nrows = 1
)

View(out_combined)



### This plot shows the dynamics of the human and mosquito populations under the combined intervention scenario, where both vector control and vaccination are implemented simultaneously. 


## 3b Model Assumptions:

## 1. Closed population: Both human and mosquito populations are assumed to be closed, meaning there is no migration or emigration.
## 2. Homogeneous mixing: Humans and mosquitoes interact randomly and uniformly.
## 3. Constant birth and death rates: Mosquito birth and death rates are constant.
## 4. No vertical transmission: Mosquitoes are not born infected.
## 5. Recovery and immunity: Humans recover from infection and become susceptible again after 25 days.
## 6. Vaccination: Vaccination provides temporary immunity, and individuals become susceptible again after 30 days.


### 3c.	 Could you give an example of an intervention that you know has a similar mode of action to the interventions described here?


## Scenario 1: The intervention is similar to the protection ITN can over which prevent human-to-mosquitoes contact and also eliminate mosquitoes that come in contact with the ITN

## Scenario 2: The intervention in scenario 2 is smilar to Seasonal Malaria Chemoprevention which offers protection for children over 28 days





## 7.	Quantify intervention impact: Run simulations for each intervention and in combination and then calculate infections averted compared to the baseline at the end of the simulation period. What is your interpretation of the results.

# Calculate infections averted
infections_averted_vector_control <- out$Ih[nrow(out)] - out_vector_control$Ih[nrow(out_vector_control)]
infections_averted_vaccination <- out$Ih[nrow(out)] - out_vaccination$Ih[nrow(out_vaccination)]
infections_averted_combined <- out$Ih[nrow(out)] - out_combined$Ih[nrow(out_combined)]

print(infections_averted_vector_control)
print(infections_averted_vaccination)
print(infections_averted_combined)


# Calculate infections averted
infections_averted_vector_control <- out$Ih[nrow(out)] - out_vector_control$Ih[nrow(out_vector_control)]
infections_averted_vaccination <- out$Ih[nrow(out)] - out_vaccination$Ih[nrow(out_vaccination)]
infections_averted_combined <- out$Ih[nrow(out)] - out_combined$Ih[nrow(out_combined)]

# Create a dataframe for plotting
infections_averted_df <- data.frame(
  Intervention = c("Vector Control", "Vaccination", "Combined"),
  Infections_Averted = c(infections_averted_vector_control, infections_averted_vaccination, infections_averted_combined)
)

# Plot the infections averted using plotly

plot_ly(infections_averted_df, x = ~Intervention, y = ~Infections_Averted, type = 'bar', color = ~Intervention, colors = c("Vector Control" = "blue", "Vaccination" = "orange", "Combined" = "green")) %>%
  layout(title = "Infections Averted by Intervention", xaxis = list(title = "Intervention"), yaxis = list(title = "Infections Averted"))