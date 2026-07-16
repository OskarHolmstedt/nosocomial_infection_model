# Basic reproduction number for the spatial SIS model.
# Sums per-contact transmission rates across spatial scales,
# weighted by number of contacts at each level, divided by recovery rate.
r_0 <- function(theta) {
  b  <- theta$beta[1:3]
  ss <- theta$spatial_sizes
  n  <- c(ss[1] - 1,      # same room, excluding self
          ss[2] - ss[1],  # same ward, different room
          ss[3] - ss[2])  # different ward
  sum(b * n) / theta$P$Dis[[1L]]
}

# Compute prevalence statistics from a simulate_outbreak() result.
#
# Returns a list with:
#   $series         data.frame of I, S, N, prevalence, I_to_S per time step
#   $mean_prevalence  time-averaged I/N
#   $mean_I_to_S      time-averaged I/S ratio
#   $sis_endemic    list with R0, endemic_prev, endemic_I_to_S from the
#                   deterministic SIS equilibrium (1 - 1/R0)
outbreak_prevalence <- function(OutbreakData) {
  theta  <- OutbreakData$theta
  State  <- OutbreakData$State
  n_beds <- theta$spatial_sizes[1] * theta$spatial_sizes[2]

  # Restrict to bed rows only (exclude room-aggregate nodes)
  State <- State[seq_len(n_beds), , drop = FALSE]

  N <- nrow(State)
  I <- as.numeric(colSums(State))
  S <- N - I

  series <- data.frame(
    t          = seq_along(I) - 1L,
    I          = I,
    S          = S,
    N          = N,
    prevalence = I / N,
    I_to_S     = ifelse(S > 0, I / S, NA_real_)
  )

  R0           <- r_0(theta)
  endemic_prev <- max(0, 1 - 1 / R0)

  list(
    series          = series,
    mean_prevalence = mean(series$prevalence),
    mean_I_to_S     = mean(series$I_to_S, na.rm = TRUE),
    sis_endemic     = list(
      R0             = R0,
      endemic_prev   = endemic_prev,
      endemic_I_to_S = if (endemic_prev < 1) endemic_prev / (1 - endemic_prev) else Inf
    )
  )
}

# Plot simulated prevalence over time against the SIS endemic equilibrium.
plot_prevalence <- function(res) {
  library(ggplot2)

  sis_prev <- res$sis_endemic$endemic_prev
  sim_prev <- res$mean_prevalence
  T        <- max(res$series$t)

  ggplot(res$series, aes(x = t, y = prevalence)) +
    geom_line(colour = "#2c7bb6", linewidth = 0.5) +
    geom_hline(yintercept = sis_prev, colour = "#d7191c",
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = sim_prev, colour = "#1a9641",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = T * 0.7, y = sis_prev + 0.03,
             label = sprintf("SIS endemic = %.0f%%", sis_prev * 100),
             colour = "#d7191c", size = 3.5) +
    annotate("text", x = T * 0.7, y = sim_prev - 0.03,
             label = sprintf("Sim. mean = %.0f%%", sim_prev * 100),
             colour = "#1a9641", size = 3.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      x     = "Day",
      y     = "Prevalence (I / N)",
      title = sprintf("R\u2080 = %.1f \u2014 simulated prevalence vs. SIS endemic",
                      res$sis_endemic$R0)
    )
}
