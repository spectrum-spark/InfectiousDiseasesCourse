
#############################################################
## MODELLING VECTOR-BORNE DISEASES PRACTICAL SESSION PART 1##
#############################################################

# A model for dengue transmission - Part 1
# Some R code to numerically solve a set of ordinary differential equations
# for a simple model of dengue. Results are then plotted.

# Library imports
library(dplyr)
library(deSolve)
library(ggplot2)

# Model parameters
parameters <- c(
  bites = 0.63, # number of bites per mosquito per day
  T_HM = 0.26, # probability of transmission - mosquito to human
  T_MH = 0.26, # Probability of transmission - human to mosquito 
  infectious_period = 5, # days
  mu_M = 1/14 # Mosquito birth and death rate
)

# Time window
start_date = as.Date("2008-11-02")
end_date = as.Date("2009-05-31")
times = seq(start_date, end_date, by = 1)

# Initial conditions
human_population <- 150000 # population of Cairns in 2008
initial_human_susceptible <- 149998
initial_human_infectious <- 2
initial_human_recovered <- human_population - initial_human_susceptible - initial_human_infectious
mosquito_population <- 225000 
initial_mosquito_susceptible <- 224990
initial_mosquito_infectious <- mosquito_population - initial_mosquito_susceptible

# Define model compartments, gather into a single variable to capture the state of the system and set to initial values
state <- c(Susceptible_human = initial_human_susceptible,
           Infectious_human = initial_human_infectious,
           Recovered_human = initial_human_recovered,
           Susceptible_mosquito = initial_mosquito_susceptible,
           Infectious_mosquito = initial_mosquito_infectious)


# Model function
dengue_base <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # Calculate the total population sizes
    Total_human_population <- Susceptible_human + Infectious_human + Recovered_human
    Total_mosquito_population <- Susceptible_mosquito + Infectious_mosquito
    
    # Calculate the average force of infection imposed on each susceptible human
    force_of_infection_on_human <- bites * T_HM * Infectious_mosquito / Total_human_population
    
    # Calculate the average force of infection imposed on each susceptible mosquito
    force_of_infection_on_mosquito <- bites * T_MH * Infectious_human / Total_human_population
    
    # Calculate the (net) instantaneous change in each compartment
    Susceptible_human_change <- -force_of_infection_on_human * Susceptible_human
    Infectious_human_change <- force_of_infection_on_human * Susceptible_human - Infectious_human / infectious_period
    Recovered_human_change <- Infectious_human / infectious_period
    Susceptible_mosquito_change <- mu_M * Total_mosquito_population - (force_of_infection_on_mosquito + mu_M) * Susceptible_mosquito 
    Infectious_mosquito_change <- force_of_infection_on_mosquito * Susceptible_mosquito - mu_M * Infectious_mosquito
    
    # Return net changes as list
    return(list(
      c(
        Susceptible_human_change,
        Infectious_human_change,
        Recovered_human_change,
        Susceptible_mosquito_change,
        Infectious_mosquito_change
      )
    ))
  })
}

# Solve model
out_dengue_base <- ode(y = state, times = as.numeric(times - times[1]), func = dengue_base, parms = parameters)

# Plot solution
plot(out_dengue_base, 
     main = c("Susceptible humans", "Infectious humans", "Recovered humans",
              "Susceptible mosquitoes", "Infectious mosquitoes"), 
     xlab = "Time", ylab = c("population size"))

# Convert ODE output to a tibble
dengue_base_infections <- tibble::as_tibble(out_dengue_base[, c("time", "Susceptible_human", "Infectious_human", "Recovered_human")])

# Create columns for new infections and cumulative new infections, 
# Select every 7th row and calculate the difference
dengue_base_infections <- dengue_base_infections  |>
  mutate(Susceptible_human_change = -c(0, diff(Susceptible_human)),
         Cumulative_human_infections = cumsum(Susceptible_human_change)) |>
  filter(time %% 7 == 0) |>
  select(time, Cumulative_human_infections) |>
  mutate(Weekly_human_infections = c(0, diff(Cumulative_human_infections)))    


# Plot solution
ggplot(dengue_base_infections, aes(x = time, y = Weekly_human_infections)) +
  geom_line() +
  labs(
    x = "Time",
    y = "Number of infections",
    title = "Weekly human new infections"
  ) +
  theme_bw()

