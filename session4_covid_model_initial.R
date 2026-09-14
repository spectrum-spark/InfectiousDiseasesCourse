# Session 4 cumulative checkpoint: initial SEIR model

#| message: false
#| warning: false

library(dplyr)
library(deSolve)
library(ggplot2)

first_wave <- readr::read_csv(
    "first_wave_TH.csv",
    col_types = readr::cols(
        Date = readr::col_date(),
        Cases = readr::col_double(),
        Cumulative_cases = readr::col_double()
    )
)

summary(first_wave)

parameters <- c(
    R0 = 4,
    latent_period = 5,
    infectious_period = 6
)

# Initial conditions
Total_population <- 6.6e7
Initial_exposed <- 0
Initial_infected <- 20
Initial_recovered <- 0
Initial_susceptible <- Total_population - Initial_exposed - Initial_infected - Initial_recovered

# State variables
state <- c(
    Susceptible = Initial_susceptible,
    Exposed = Initial_exposed,
    Infectious = Initial_infected,
    Recovered = Initial_recovered
)

# Time window
start_date <- as.Date("2020-03-07")
end_date <- as.Date("2020-03-26")
times <- seq(start_date, end_date, by = 1)

# Model function
covid_base <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        # Calculate the total population size
        Total_population <- Susceptible + Exposed + Infectious + Recovered

        # Calculate the average force of infection imposed on each susceptible individual
        force_of_infection <- R0 * Infectious / (Total_population * infectious_period)

        # Calculate the net (instantaneous) change in each state variable
        Susceptible_change <- -force_of_infection * Susceptible
        Exposed_change <- force_of_infection * Susceptible - Exposed / latent_period
        Infectious_change <- Exposed / latent_period - Infectious / infectious_period
        Recovered_change <- Infectious / infectious_period

        # Return net changes as list
        return(list(
            c(
                Susceptible_change,
                Exposed_change,
                Infectious_change,
                Recovered_change
            )
        ))
    })
}

# Solve model
out <- ode(y = state, times = as.numeric(times - times[1]), func = covid_base, parms = parameters)
# Plot solution
par(mar = c(1, 1, 1, 1)) # reduce the margins of the plot in order to fit it in the panel
plot(out)


