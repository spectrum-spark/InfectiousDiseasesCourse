# Mathematical Modelling of Infectious Diseases
###############################################################
## INTRODUCTION TO MATHEMATICAL MODELLING PRACTICAL SESSION  ##
###############################################################

# Record simulated epidemic results.
result <- matrix(NA, nrow = 10, ncol = 5)
result[, 1] <- 1:10 # week

# Replace these values with your own S, I and R counts.
result[, 2] <- c(19, 17, 13, 9, 5, 2, 1, 1, 1, 1) # S
result[, 3] <- c(1, 2, 4, 4, 4, 3, 1, 0, 0, 0) # I
result[, 4] <- c(0, 1, 3, 7, 11, 15, 18, 19, 19, 19) # R

# Check that S + I + R equals the total population each week.
result[, 5] <- rowSums(result[, 2:4])
result[, 5]

library(ggplot2)

result_data <- result |>
  tibble::as_tibble(.name_repair = ~ c(
    "week", "susceptible", "infectious", "recovered", "total"
  ))

ggplot(result_data, aes(x = week, y = infectious)) +
  geom_point(colour = "red") +
  labs(
    x = "Time in weeks",
    y = "Infectious people",
    title = "My simulated epidemic"
  )
