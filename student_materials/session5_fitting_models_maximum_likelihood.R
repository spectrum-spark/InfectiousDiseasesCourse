# Session 5 cumulative checkpoint: maximum-likelihood fitting

#| message: false
# SPARK Modelling Short course

#########################
## FITTING MODELS IN R ##
#########################
library(dplyr)
library(deSolve)
library(ggplot2)

# Data imports and filtering
first_wave <- readr::read_csv(
    "first_wave_TH.csv",
    col_types = readr::cols(
        Date = readr::col_date(),
        Cases = readr::col_double(),
        Cumulative_cases = readr::col_double()
    )
)
# Time window
start_date <- as.Date("2020-03-07")
end_date <- as.Date("2020-03-26")

# Filter data to capture period prior to interventions. This has already been done in the last session.
uncontrolled_period <- first_wave |>
    filter(Date >= start_date, Date <= end_date)

# Plot the filtered data
ggplot(uncontrolled_period) +
    geom_col(aes(x = Date, y = Cases), width = 1, fill = "dodgerblue2", colour = "blue") +
    ylab("Daily cases") +
    xlab("") +
    ggtitle("Uncontrolled first wave of COVID-19 in Thailand 1st to 24th March, 2020") +
    theme_bw()

## Define model as per the previous session
# Time window
times <- seq(start_date, end_date, by = 1)

# Model parameters
parameters <- c(
    R0 = 4,
    latent_period = 5,
    infectious_period = 6
)

# Initial conditions
Total_population <- 6.6e7 # Population of Thailand
Initial_exposed <- 0
Initial_infectious <- 20 # Initial infectious seed
Initial_recovered <- 0
Initial_susceptible <- Total_population - Initial_exposed - Initial_infectious - Initial_recovered

# State variables
state <- c(
    Susceptible = Initial_susceptible,
    Exposed = Initial_exposed,
    Infectious = Initial_infectious,
    Recovered = Initial_recovered
)


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

solve_base_model <- function(initial_state = state,
                             model_times = times,
                             model_function = covid_base,
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
            Prevalence = Exposed + Infectious,
            Incidence = Exposed / model_parameters["latent_period"],
            Cumulative_incidence = cumsum(Incidence) + Incidence[1],
            Population = Susceptible + Exposed + Infectious + Recovered,
            Date = model_times
        )

    return(out)
}

# now we use the model developed in the last session and solve it
# for an initial set of parameters, an initial state and using times
# from just 7th - 26th  of March 2020
out_init <- solve_base_model(
    initial_state = state,
    model_times = times,
    model_function = covid_base,
    model_parameters = parameters
)

# Plot the filtered data and the first model "guess"
ggplot(uncontrolled_period) +
    geom_col(aes(x = Date, y = Cases), width = 1, fill = "dodgerblue2", colour = "blue") +
    geom_point(data = out_init, aes(x = Date, y = Incidence), size = 2, colour = "firebrick2") +
    ylab("Daily cases") +
    xlab("") +
    ggtitle("Thailand's First Wave, Jan-Jul 2020") +
    theme_bw()

# find the sum of the residuals squared for day 2 until day 24
SSQ_initial <- sum((uncontrolled_period$Cases[-1] - out_init$Incidence[-1])^2)
SSQ_initial

# transform the parameters
# The search / optimization algorithm we employ searches over
# the range (-Inf, Inf) for each variable. In our application
# we are only interested in solutions in the range (0, Inf)
# for R0 and I(0). Therefore, it
# would be a good idea to apply a transformation to our
# search space so that we are not wasting time exploring
# infeasible regions of parameter space.
initial_transformed_parameters <- log(c(
    "R0" = 4,
    "Initial_infectious" = 20
))

#### Fitting routine ####

# We have selected R0 and I(0) as our free parameters
# We need a function that accepts R0 and I(0) as arguments,
# solves the model using these inputs, and calculates the
# sum of squares of the data given these parameters

SSQ_function <- function(transformed_parameters,
                         data = uncontrolled_period$Cases[-1],
                         state_base = state,
                         model_times = times,
                         model_function = covid_base,
                         model_parameters = parameters) {

    # Untransform parameters
    R0 <- exp(transformed_parameters["R0"])
    Initial_infectious <- exp(transformed_parameters["Initial_infectious"])

    # Calculate updated susceptible population
    Initial_susceptible <- state_base["Susceptible"] + state_base["Infectious"] - Initial_infectious

    # Overwrite baseline parameters with proposed parameters
    model_parameters["R0"] <- R0
    state_base["Susceptible"] <- Initial_susceptible
    state_base["Infectious"] <- Initial_infectious

    # Solve model with updated parameters
    out <- solve_base_model(
        state_base,
        model_times,
        model_function,
        model_parameters
    )
    return(
        sum((uncontrolled_period$Cases[-1] - out$Incidence[-1])^2)
    )
}

SSQ_initial_from_function <- SSQ_function(initial_transformed_parameters)
SSQ_initial_from_function
SSQ_initial

#### Optimal parameters ####

# Use the optim function to determine the parameters that
# minimise the sum of squared residuals
# We will use the Nelder-Mead optimization solver
optim_NM <- optim(
    par = initial_transformed_parameters,
    fn = SSQ_function,
    control = list(maxit = 500),
    method = "Nelder-Mead",
    hessian = TRUE
)

# Check for convergence (always code 0 for NM)
optim_NM$convergence # 0 - converged; 1 - failed to converge

# Inspect solution
optim_NM$par

# Back-transform parameters
optimum_parameters <- exp(optim_NM$par)
optimum_parameters

# Inspect the optimal sum of squared residuals
optimum_SSQ <- optim_NM$value
optimum_SSQ

#### Optimal model solution ####

# Create the optimal parameter and state vectors
optimal_parameters <- parameters
optimal_parameters["R0"] <- optimum_parameters["R0"]

optimal_initial_state <- state
optimal_initial_state["Susceptible"] <- state["Susceptible"] + state["Infectious"] - optimum_parameters["Initial_infectious"]
optimal_initial_state["Infectious"] <- optimum_parameters["Initial_infectious"]

# Solve the model given the optimal parameters and initial conditions
optimal_solution <- solve_base_model(
    initial_state = optimal_initial_state,
    model_times = times,
    model_function = covid_base,
    model_parameters = optimal_parameters
)


# Plot the optimal solution
ggplot(first_wave) +
    geom_col(aes(x = Date, y = Cases), width = 1, fill = "dodgerblue2", colour = "blue") +
    geom_point(data = optimal_solution, aes(x = Date, y = Incidence), size = 2, colour = "firebrick2") +
    ylab("Daily cases") +
    xlab("") +
    ggtitle("Fit to unmitigated period") +
    theme_bw()

#### Fitting routine ####

# We have selected R0 and I(0) as our free parameters
# We need a function that accepts R0 and I(0) as arguments,
# solves the model using these inputs, and calculates the
# negative log-likelihood of the data given these parameters


negative_log_likelihood <- function(transformed_parameters,
                                    data = uncontrolled_period$Cases[-1],
                                    state_base = state,
                                    model_times = times,
                                    model_function = covid_base,
                                    model_parameters = parameters) {

    # Untransform parameters
    R0 <- exp(transformed_parameters["R0"])
    Initial_infectious <- exp(transformed_parameters["Initial_infectious"])
    # Calculate updated susceptible population
    Initial_susceptible <- state_base["Susceptible"] + state_base["Infectious"] - Initial_infectious

    # Overwrite baseline parameters with proposed parameters
    model_parameters["R0"] <- R0
    state_base["Susceptible"] <- Initial_susceptible
    state_base["Infectious"] <- Initial_infectious

    # Solve model with updated parameters
    out <- solve_base_model(
        state_base,
        model_times,
        model_function,
        model_parameters
    )

    return(-sum(dpois(
        x = data,
        lambda = out$Incidence[-1],
        log = TRUE
    )))
}

#### Optimal parameters ####

# Use the optim function to determine the parameters that
# minimize the negative log-likelihood (i.e., maximize
# the likelihood)
# We will use the Nelder-Mead optimization solver
optim_NM <- optim(
    par = initial_transformed_parameters,
    fn = negative_log_likelihood,
    control = list(maxit = 500),
    method = "Nelder-Mead",
    hessian = TRUE
)

# Check for convergence (always code 0 for NM)
optim_NM$convergence # 0 - converged; 1 - failed to converge

# Inspect solution
optim_NM$par

# Back-transform parameters
optimum_parameters <- exp(optim_NM$par)
optimum_parameters

optimal_parameters <- parameters
optimal_parameters["R0"] <- optimum_parameters["R0"]

optimal_initial_state <- state
optimal_initial_state["Susceptible"] <- state["Susceptible"] + state["Infectious"] - optimum_parameters["Initial_infectious"]
optimal_initial_state["Infectious"] <- optimum_parameters["Initial_infectious"]

optimal_solution <- solve_base_model(
    initial_state = optimal_initial_state,
    model_times = times,
    model_function = covid_base,
    model_parameters = optimal_parameters
)


