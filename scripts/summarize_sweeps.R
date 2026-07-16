source(here::here("scripts", "_load_project.R"))

run_validation <- function(n_reps, T, theta, N_samples = 5000, burn_in = 1000, vary = NULL) {
  results <- parallel::mclapply(seq_len(n_reps), function(i) {
    OD  <- simulate_outbreak(T, theta)
    fit <- mcmc(OD, N_samples = N_samples, burn_in = burn_in)
    sc  <- ancestry_score(OD$ObsRec$Anc2, fit$anc, OD$ObsRec$Adm, OD$ObsRec$PTest)
    data.frame(
      rep           = i,
      n_cases       = nrow(OD$ObsRec),
      mode_accuracy = sc$mode_accuracy,
      mean_p_true   = sc$mean_p_true,
      log_skill     = sc$log_skill,
      value         = if (!is.null(vary)) vary else NA_real_
    )
  }, mc.cores = parallel::detectCores() - 1L)
  do.call(rbind, results)
}

sweep_s <- do.call(rbind, lapply(c(0.05, 0.1, 0.2, 0.4), function(s) {
  theta_s <- theta
  theta_s$P$Test <- rep(s, length(theta$P$Test))
  run_validation(n_reps = 20, T = 100, theta = theta_s, vary = s)
}))

# --- Plots -----------------------------------------------------------
library(tidyr)

sweep_long <- sweep_s |>
  pivot_longer(
    cols      = c(mode_accuracy, mean_p_true, log_skill),
    names_to  = "metric",
    values_to = "score"
  ) |>
  transform(
    test_rate = factor(value),
    metric    = factor(metric,
                       levels = c("mode_accuracy", "mean_p_true", "log_skill"),
                       labels = c("Mode accuracy", "Mean P(true ancestor)", "Log-skill"))
  )

# 1. Performance metrics vs testing rate (faceted)
ggplot(sweep_long, aes(x = test_rate, y = score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Daily testing probability", y = NULL,
       title = "Inference quality vs. testing rate")

# 2. Detected cases vs testing rate
ggplot(sweep_s, aes(x = factor(value), y = n_cases)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  labs(x = "Daily testing probability", y = "Detected cases",
       title = "Outbreak size (detected) vs. testing rate")

# 3. Performance vs detected cases — does more data help regardless of how it arises?
ggplot(sweep_long, aes(x = n_cases, y = score, colour = factor(value))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Detected cases", y = NULL, colour = "Test rate",
       title = "Inference quality vs. outbreak size")
