# Session 9 cumulative checkpoint: initial intervention model

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

ggplot(first_wave, aes(x = Date, y = Cases)) +
    geom_col(fill = "blue") +
    ylab("Daily cases") +
    xlab("") +
    ggtitle("Thailand's First Wave, Jan-Jun 2020") +
    theme_bw()

fitting_start_date <- as.Date("2020-03-07") # start of fitting window
emergency_start_date <- as.Date("2020-03-26") # Start of Emergency declaration
curfew_start_date <- as.Date("2020-04-03") # Start of Curfew
fitting_end_date <- as.Date("2020-05-13") # end of fitting window

outbreak <- first_wave |>
    filter(Date >= fitting_start_date, Date <= fitting_end_date)

parameters <- c(
    R0 = 3.970588, # from the NM fitting session
    latent_period = 5, # Rounded from Trauer, J.M., Lydeamore, M.J., Dalton, G.W. et al. Understanding how Victoria, Australia gained control of its second COVID-19 wave. Nat Commun 12, 6266 (2021). https://doi.org/10.1038/s41467-021-26558-4
    infectious_period = 6, # Rounded from Trauer et al
    emergency_efficacy = 0.70, # guess, assuming Reff = (1-efficacy) * R0
    curfew_efficacy = 0.80 # guess
)

# Initial conditions
Total_population <- 6.6e7 # Roughly population of Thailand
Initial_exposed <- 0 # Simplifying assumption
Initial_infectious <- 30.356868 # from the NM fitting session
Initial_recovered <- 0 # simplifying assumption
Initial_susceptible <- Total_population - Initial_exposed - Initial_infectious - Initial_recovered


# State variables
state <- c(
    Susceptible = Initial_susceptible,
    Exposed = Initial_exposed,
    Infectious = Initial_infectious,
    Recovered = Initial_recovered
)

# Time window
times <- seq(fitting_start_date, fitting_end_date, by = 1)

covid_intervention <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

        # Calculate the total population size
        Total_population <- Susceptible + Exposed + Infectious + Recovered

        # Calculate intervention efficacy
        if (t < as.numeric(emergency_start_date - fitting_start_date)) {
            intervention_efficacy <- 0
        } else if (t >= as.numeric(emergency_start_date - fitting_start_date) &&
            t < as.numeric(curfew_start_date - fitting_start_date)) {
            intervention_efficacy <- emergency_efficacy
        } else {
            intervention_efficacy <- curfew_efficacy
        }


        # Calculate the effective reproduction number in the presence of interventions
        Reff <- (1 - intervention_efficacy) * R0 # assuming everyone starts susceptible

        # Calculate the average force of infection imposed on each susceptible individual
        force_of_infection <- Reff * Infectious / (Total_population * infectious_period)

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

# Wrapper function to solve model and tidy up output
solve_intervention_model <- function(initial_state = state,
                                     model_times = times,
                                     model_function = covid_intervention,
                                     model_parameters = parameters) {
    out <- ode(
        y = initial_state,
        times = as.numeric(model_times - model_times[1]),
        func = model_function,
        parms = model_parameters
    )


    # Calculate the prevalence, incidence and cumulative incidence (for comparison with data)
    out <- tibble::as_tibble(out) |>
        mutate(
            Incidence = Exposed / model_parameters["latent_period"],
            Cumulative_incidence = cumsum(Incidence) + Incidence[1],
            Population = Susceptible + Exposed + Infectious + Recovered,
            Prevalence = (Exposed + Infectious)/Population,
            Date = model_times
        )

    return(tibble::as_tibble(out))
}

# Run the model for out initial set of parameters
out_init <- solve_intervention_model()


# Make sure to add observational uncertainty
out_init <- out_init |>
    mutate(
        lower50 = qpois(p = 0.25, lambda = Incidence), # 50% confidence interval (i.e., 25 - 75th centiles)
        upper50 = qpois(p = 0.75, lambda = Incidence),
        lower95 = qpois(p = 0.025, lambda = Incidence), # 95% confidence interval (i.e., 2.5 - 97.5th centiles)
        upper95 = qpois(p = 0.975, lambda = Incidence)
    )

# Plot initial estimate
ggplot(outbreak) +
    geom_col(aes(x = Date, y = Cases), width = 1, fill = "dodgerblue2", colour = "blue") +
    geom_ribbon(
        data = out_init[-1, ], aes(x = Date, ymin = lower50, ymax = upper50),
        fill = "firebrick2", colour = "firebrick2", alpha = 0.8
    ) +
    geom_ribbon(
        data = out_init[-1, ], aes(x = Date, ymin = lower95, ymax = upper95),
        fill = "firebrick2", colour = "firebrick2", alpha = 0.5
    ) +
    ylab("Daily cases") +
    xlab("") +
    ggtitle("Thailand's First Wave, 2020") +
    theme_bw()


