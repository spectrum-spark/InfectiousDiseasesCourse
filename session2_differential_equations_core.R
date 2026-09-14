# Session 2 cumulative checkpoint: Parts 1 and 2 core material

#| echo: true
#| warning: false
library(deSolve)
library(tidyverse)

# Define model parameters and store values in a labelled vector
parameters <- c(r = 2) # Growth rate

# Specify the initial conditions, define state variables and assign initial values
initial_population <- 1 # Start with an initial population size N(0) = 1
state <- c(N = initial_population) # Specify the state variables and initialize them to their initial values

# We specify the system of differential equations as a function in R
# Here we specify a function with three arguments:
# 1. t : time, which is the dependent variable in our case
# 2. state : a vector of state variables
# 3. parameters : a vector of parameters (ie fixed variables)
population_exponential_growth <- function(t, state, parameters) {
    # Initialize an environment with the elements of state and parameters as local variables
    with(as.list(c(state, parameters)), {
        dNdt <- r * N # The model equation
        list(c(dNdt)) # Return the population growth rate dNdt
    })
}

# Start with an initial time of t=0 and solve up to t=5 with time increments of 0.01
times <- seq(from = 0, to = 5, by = 0.01)

# Here we specify four arguments of the ode function:
# 1. y : the state variables
# 2. times : the time lattice over which we solve the differential equation
# 3. func : our differential equation given in terms of our state variables and model parameters
# 4. parms : our model parameters
out <- ode(y = state, times = times, func = population_exponential_growth, parms = parameters)

head(out)

ggplot(data = tibble::as_tibble(out), aes(x=time, y=N)) +
    geom_line() +
    labs(x = "Time, t", y="Population size, N(t)")

# Create a coarse time lattice for display purposes
times_coarse <- seq(from = 0, to = 5, by = 1)
# Calculate the exact solution given the model parameters and initial conditions
exact_population <- initial_population * exp(parameters["r"] * times_coarse)

exact_solution <- tibble::tibble(time = times_coarse, exact_population)

ggplot() +
    geom_line(data = tibble::as_tibble(out), aes(x=time, y=N, colour="Numeric")) +
    geom_point(data = exact_solution, aes(x=time, y=exact_population, colour="Exact")) +
    labs(x = "Time, t", y="Population size, N(t)", colour="") +
    scale_color_manual(values=c("red","black"))

sir_model <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        N <- S + I + R
        dSdt <- -beta * S * I / N
        dIdt <- beta * S * I / N - gamma * I
        dRdt <- gamma * I

        return(list(c(dSdt, dIdt, dRdt)))
    })
}

parameters <- c(beta = 5, gamma = 1)
initcond <- c(S = 199, I = 1, R = 0)

times <- seq(0, 8, by = 0.1)

out <- ode(y = initcond, times = times, func = sir_model, parms = parameters)


# ggplot likes "long" data, so we pivot our data from wide (1 column per state)
# to long (1 row per state-value pair).
plot_data <- tibble::as_tibble(out) |>
    pivot_longer(!time, names_to = "state") |>
    mutate(state = factor(state, levels=c("S","I","R")))

# x-axis: Time
# y-axis: value of the state
# Colour the lines by state
# Change the type of line depending on state
ggplot(data = plot_data, aes(x=time, y=value, colour=state, linetype=state)) + 
    geom_line(linewidth = 0.8) + # Plot the lines
    geom_point(data = plot_data |> slice(seq(1, n(), 10)), size = 2) + #Overlay points on the line curves
    ggtitle("SIR Model")


