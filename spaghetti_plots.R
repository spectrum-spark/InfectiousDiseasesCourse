# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)

# Function for one simulation using Gillespie algorithm
simulate_SIR <- function(beta, gamma, S0, I0, R0, max_time) {
  time <- 0
  S <- S0
  I <- I0
  R <- R0
  
  result <- data.frame(time = time, S = S, I = I, R = R)
  
  while (time < max_time && I > 0) {
    rate_infection <- beta * S * I / (S + I + R)
    rate_recovery <- gamma * I
    rate_total <- rate_infection + rate_recovery
    
    if (rate_total == 0) break
    
    # Time to next event
    time <- time + rexp(1, rate_total)
    
    # Determine which event occurs
    if (runif(1) < rate_infection / rate_total) {
      # Infection
      S <- S - 1
      I <- I + 1
    } else {
      # Recovery
      I <- I - 1
      R <- R + 1
    }
    
    result <- rbind(result, data.frame(time = time, S = S, I = I, R = R))
  }
  
  return(result)
}

# Parameters
n_iter <- 100
beta <- 0.2
gamma <- 0.1
S0 <- 990
I0 <- 10
R0 <- 0
max_time <- 100
N <- S0 + I0 + R0

# Run simulations
all_runs <- lapply(1:n_iter, function(i) {
  sim <- simulate_SIR(beta, gamma, S0, I0, R0, max_time)
  sim$run <- i
  return(sim)
})

df_all <- bind_rows(all_runs)

# Convert to long format for ggplot
df_long <- df_all |>
  pivot_longer(cols = c("S", "I", "R"), names_to = "Compartment", values_to = "Count") |>
  mutate(Compartment = factor(Compartment, levels = c("S", "I", "R")))

# Spaghetti plot
ggplot(df_long, aes(x = time, y = Count, group = interaction(run, Compartment), color = Compartment)) +
  geom_line(alpha = 0.1) +
  labs(title = "Stochastic SIR Simulations",
       x = "Time (Days)",
       y = "Number of Individuals") +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha=1)))

cowplot::save_plot("sir_simulations.png", last_plot(), base_height = 3, base_width = 6, bg = "white")



# ODE model definition
sir_ode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# Solve ODE
ode_times <- seq(0, max_time, by = 0.1)
init_state <- c(S = S0, I = I0, R = R0)
parameters <- c(beta = beta, gamma = gamma)
ode_out <- as.data.frame(ode(y = init_state, times = ode_times, func = sir_ode, parms = parameters))

# Convert ODE to long format
ode_long <- ode_out %>%
  pivot_longer(cols = c("S", "I", "R"), names_to = "Compartment", values_to = "Count")

# Plot with ODE solution overlay
ggplot() +
  geom_line(data = df_long, aes(x = time, y = Count, group = interaction(run, Compartment), color = Compartment), alpha = 0.1) +
  geom_line(data = ode_long, aes(x = time, y = Count, color = Compartment), size = 1.5) +
  labs(title = "Stochastic SIR Simulations and ODE Solution",
       x = "Time (Days)",
       y = "Number of Individuals") +
  theme_minimal()

cowplot::save_plot("sir_simulations_ode.png", last_plot(), base_height = 3, base_width = 6, bg = "white")
