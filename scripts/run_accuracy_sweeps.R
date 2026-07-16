library(parallel)
library(tidyverse)
source(here::here("scripts", "_load_project.R"))

cache_dir  <- here::here("intermediate", "sweeps")
table_dir  <- here::here("results", "tables")
figure_dir <- here::here("results", "figures")
invisible(lapply(c(cache_dir, table_dir, figure_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

T_demo    <- 45
theta_demo <- local({
  th           <- theta
  th$P$Init    <- c(1, rep(0, NumBeds - 1), rep(0, NumRooms))
  th
})

# ── Trial function ────────────────────────────────────────────────────────────
# Simulate one outbreak; skip MCMC if n_obs falls outside [obs_min, obs_max].
# This conditions all conditions on the same observed outbreak size, removing
# the selection bias that arises when filtering only on a minimum.
# Note: do NOT use return(NULL) inside a tryCatch block — it exits the enclosing
# function, not just the tryCatch expression. Use NULL as the last expression instead.
run_trial_windowed <- function(theta, T, N_samples = 2000,
                                obs_min = 8, obs_max = 12) {
  od    <- simulate_outbreak(T, theta)
  n_obs <- nrow(od$ObsRec)
  if (n_obs < obs_min || n_obs > obs_max) return(NULL)  # safe here — not inside tryCatch
  Y  <- mcmc(od, N_samples = N_samples)
  sc <- ancestry_score(od$ObsRec$Anc2, Y$anc,
                       od$ObsRec$Adm, od$ObsRec$PTest)
  data.frame(
    n_obs       = n_obs,
    n_total     = nrow(od$CaseRec),
    mode_acc    = sc$mode_accuracy,
    mean_p_true = sc$mean_p_true,
    log_skill   = sc$log_skill
  )
}

# ── Sweep ─────────────────────────────────────────────────────────────────────
test_rates  <- seq(0.05, 0.975, by = 0.025)  # up to ~1; rate=1.0 causes MCMC errors
N_rep       <- 35
N_samples   <- 2000
max_att     <- 500
n_cores     <- min(detectCores() - 1L, length(test_rates))

set.seed(7823, kind = "L'Ecuyer-CMRG")

results <- do.call(rbind, mclapply(test_rates, function(rate) {
  theta_r        <- theta_demo
  theta_r$P$Test <- rep(rate, length(theta_demo$P$Test))
  valid   <- list()
  attempt <- 0
  while (length(valid) < N_rep && attempt < max_att) {
    attempt <- attempt + 1
    res <- tryCatch(
      run_trial_windowed(theta_r, T_demo, N_samples, obs_min = 8, obs_max = 12),
      error = function(e) NULL
    )
    if (!is.null(res)) valid[[length(valid) + 1]] <- cbind(test_rate = rate, res)
  }
  if (length(valid)) do.call(rbind, valid)
}, mc.cores = n_cores))

# ── Save data ─────────────────────────────────────────────────────────────────
dir.create("data", showWarnings = FALSE)
saveRDS(results, file.path(cache_dir, "sweep_testing_rate.rds"))
write.csv(results, file.path(table_dir, "sweep_testing_rate.csv"), row.names = FALSE)

# ── Plot ──────────────────────────────────────────────────────────────────────
# SE must be computed before mode_acc is overwritten in summarise
means <- results |>
  group_by(test_rate) |>
  summarise(
    se       = sd(mode_acc) / sqrt(n()),
    mode_acc = mean(mode_acc),
    .groups  = "drop"
  )

ggplot(means, aes(x = test_rate, y = mode_acc)) +
  geom_line(data = dplyr::mutate(means, mode_acc = mode_acc + se),
            aes(x = test_rate, y = mode_acc),
            color = "steelblue", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE) +
  geom_line(data = dplyr::mutate(means, mode_acc = mode_acc - se),
            aes(x = test_rate, y = mode_acc),
            color = "steelblue", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE) +
  geom_smooth(data = results,
              aes(x = test_rate, y = mode_acc),
              method = "loess", span = 0.4, se = FALSE,
              color = "steelblue", linewidth = 1,
              inherit.aes = FALSE) +
  geom_point(size = 2.5, color = "steelblue") +
  scale_x_continuous(breaks = seq(0.0, 1.0, by = 0.10)) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x        = "Testing rate (per bed per day)",
    y        = "Mode accuracy",
    title    = "Ancestry mode accuracy vs. testing rate",
    subtitle = sprintf(
      "%d reps x %d MCMC samples, window %d-%d obs, T=%d, two-ward hospital",
      N_rep, N_samples, 8L, 12L, T_demo)
  )

dir.create("figures", showWarnings = FALSE)
ggsave(file.path(figure_dir, "sweep_testing_rate.pdf"), width = 8, height = 5)
ggsave(file.path(figure_dir, "sweep_testing_rate.png"), width = 8, height = 5, dpi = 300)

# ── Mutation rate sweep (testing rate fixed at 0.10) ──────────────────────────
# Close any open ChromoteSession before mclapply to avoid fork issues on macOS.
mu_rates   <- seq(0.25, 8, by = 0.25)
N_rep_mu   <- 35
N_samp_mu  <- 2000
max_att_mu <- 500
fix_rate   <- 0.10

n_cores_mu <- min(detectCores() - 1L, length(mu_rates))
set.seed(4492, kind = "L'Ecuyer-CMRG")

results_mu <- do.call(rbind, mclapply(mu_rates, function(mu) {
  theta_r        <- theta_demo
  theta_r$P$Test <- rep(fix_rate, length(theta_demo$P$Test))
  theta_r$mu     <- mu
  valid <- list(); attempt <- 0
  while (length(valid) < N_rep_mu && attempt < max_att_mu) {
    attempt <- attempt + 1
    res <- tryCatch(
      run_trial_windowed(theta_r, T_demo, N_samp_mu, obs_min = 8, obs_max = 12),
      error = function(e) NULL
    )
    if (!is.null(res)) valid[[length(valid) + 1]] <- cbind(mu = mu, res)
  }
  if (length(valid)) do.call(rbind, valid)
}, mc.cores = n_cores_mu))

saveRDS(results_mu, file.path(cache_dir, "sweep_mutation_rate.rds"))
write.csv(results_mu, file.path(table_dir, "sweep_mutation_rate.csv"), row.names = FALSE)

means_mu <- results_mu |>
  group_by(mu) |>
  summarise(se = sd(mode_acc) / sqrt(n()), mode_acc = mean(mode_acc), .groups = "drop")

ggplot(means_mu, aes(x = mu, y = mode_acc)) +
  geom_line(data = dplyr::mutate(means_mu, mode_acc = mode_acc + se),
            aes(x = mu, y = mode_acc), color = "steelblue",
            linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = dplyr::mutate(means_mu, mode_acc = mode_acc - se),
            aes(x = mu, y = mode_acc), color = "steelblue",
            linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE) +
  geom_smooth(data = results_mu, aes(x = mu, y = mode_acc),
              method = "loess", span = 0.5, se = FALSE,
              color = "steelblue", linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 2, color = "steelblue") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Mutation rate (expected mutations per generation)", y = "Mode accuracy",
       title = "Ancestry mode accuracy vs. mutation rate",
       subtitle = sprintf(
         "%d reps x %d MCMC samples, window 8-12 obs, testing rate = %.2f, T=%d",
         N_rep_mu, N_samp_mu, fix_rate, T_demo))

ggsave(file.path(figure_dir, "sweep_mutation_rate.pdf"), width = 8, height = 5)
ggsave(file.path(figure_dir, "sweep_mutation_rate.png"), width = 8, height = 5, dpi = 300)

# ── R0 sweep (scale all three spatial betas proportionally) ──────────────────
# Varying R0 by keeping the spatial hierarchy (beta ratio) fixed and scaling
# beta[1:3] uniformly. W, LogComp, and Odds must be recomputed from the new betas.
r0_vals    <- seq(1.2, 5.0, by = 0.2)
N_rep_r0   <- 35
N_samp_r0  <- 2000
max_att_r0 <- 500
current_r0 <- r_0(theta_demo)

n_cores_r0 <- min(detectCores() - 1L, length(r0_vals))
set.seed(3871, kind = "L'Ecuyer-CMRG")

results_r0 <- do.call(rbind, mclapply(r0_vals, function(r0_target) {
  scale       <- r0_target / current_r0
  theta_r     <- theta_demo
  theta_r$beta[1:3] <- theta_demo$beta[1:3] * scale
  W           <- weights(theta_r$beta, Contact)
  theta_r$LogComp <- log1p(-W)
  theta_r$Odds    <- W / (1 - W)

  valid <- list(); attempt <- 0
  while (length(valid) < N_rep_r0 && attempt < max_att_r0) {
    attempt <- attempt + 1
    res <- tryCatch(
      run_trial_windowed(theta_r, T_demo, N_samp_r0, obs_min = 8, obs_max = 12),
      error = function(e) NULL
    )
    if (!is.null(res)) valid[[length(valid) + 1]] <- cbind(r0 = r0_target, res)
  }
  if (length(valid)) do.call(rbind, valid)
}, mc.cores = n_cores_r0))

saveRDS(results_r0, file.path(cache_dir, "sweep_r0.rds"))
write.csv(results_r0, file.path(table_dir, "sweep_r0.csv"), row.names = FALSE)

means_r0 <- results_r0 |>
  group_by(r0) |>
  summarise(se = sd(mode_acc) / sqrt(n()), mode_acc = mean(mode_acc), .groups = "drop")

ggplot() +
  geom_jitter(data = results_r0, aes(x = r0, y = mode_acc),
              width = 0.04, height = 0, alpha = 0.12, size = 0.7,
              color = "#3498db") +
  geom_ribbon(data = means_r0,
              aes(x = r0, ymin = mode_acc - se, ymax = mode_acc + se),
              fill = "#3498db", alpha = 0.18) +
  geom_smooth(data = results_r0, aes(x = r0, y = mode_acc),
              method = "loess", span = 0.5, se = FALSE,
              color = "#3498db", linewidth = 1.3) +
  geom_point(data = means_r0, aes(x = r0, y = mode_acc),
             size = 2.2, color = "#3498db") +
  geom_vline(xintercept = current_r0, linetype = "dashed",
             color = "#e07b39", linewidth = 0.7) +
  annotate("text", x = current_r0 + 0.1, y = 0.05,
           label = sprintf("baseline\nR\u2080 = %.1f", current_r0),
           hjust = 0, size = 3.2, color = "#e07b39") +
  scale_x_continuous(breaks = seq(1, 5, by = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(add = c(0.01, 0.02))) +
  labs(x = "Basic reproduction number R\u2080", y = "Ancestry mode accuracy") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#ebebeb"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "#cccccc"),
    axis.ticks       = element_line(color = "#cccccc")
  )

ggsave(file.path(figure_dir, "sweep_r0.pdf"), width = 8, height = 5)
ggsave(file.path(figure_dir, "sweep_r0.png"), width = 8, height = 5, dpi = 300)

# ── Mutation rate sweep at higher testing rate (0.40) for comparison ──────────
fix_rate_high <- 0.40
set.seed(6614, kind = "L'Ecuyer-CMRG")

results_mu_high <- do.call(rbind, mclapply(mu_rates, function(mu) {
  theta_r        <- theta_demo
  theta_r$P$Test <- rep(fix_rate_high, length(theta_demo$P$Test))
  theta_r$mu     <- mu
  valid <- list(); attempt <- 0
  while (length(valid) < N_rep_mu && attempt < max_att_mu) {
    attempt <- attempt + 1
    res <- tryCatch(
      run_trial_windowed(theta_r, T_demo, N_samp_mu, obs_min = 8, obs_max = 12),
      error = function(e) NULL
    )
    if (!is.null(res)) valid[[length(valid) + 1]] <- cbind(mu = mu, res)
  }
  if (length(valid)) do.call(rbind, valid)
}, mc.cores = n_cores_mu))

saveRDS(results_mu_high, file.path(cache_dir, "sweep_mutation_rate_high_testing.rds"))
write.csv(results_mu_high, file.path(table_dir, "sweep_mutation_rate_high_testing.csv"), row.names = FALSE)

means_mu_high <- results_mu_high |>
  group_by(mu) |>
  summarise(se = sd(mode_acc) / sqrt(n()), mode_acc = mean(mode_acc), .groups = "drop")

results_both <- bind_rows(
  mutate(results_mu,      test_rate = "0.10"),
  mutate(results_mu_high, test_rate = "0.40")
)
means_both <- bind_rows(
  mutate(means_mu,      test_rate = "0.10"),
  mutate(means_mu_high, test_rate = "0.40")
)

ggplot(means_both, aes(x = mu, y = mode_acc, color = test_rate)) +
  geom_line(data = dplyr::mutate(means_both, mode_acc = mode_acc + se),
            aes(x = mu, y = mode_acc, color = test_rate),
            linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = dplyr::mutate(means_both, mode_acc = mode_acc - se),
            aes(x = mu, y = mode_acc, color = test_rate),
            linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE) +
  geom_smooth(data = results_both,
              aes(x = mu, y = mode_acc, color = test_rate),
              method = "loess", span = 0.5, se = FALSE,
              linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 2) +
  scale_color_manual(values = c("0.10" = "steelblue", "0.40" = "#e07b39"),
                     name = "Testing rate") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Mutation rate (expected mutations per generation)",
       y = "Mode accuracy",
       title = "Ancestry mode accuracy vs. mutation rate",
       subtitle = sprintf(
         "%d reps x %d MCMC samples, window 8-12 obs, T=%d",
         N_rep_mu, N_samp_mu, T_demo))

ggsave(file.path(figure_dir, "sweep_mutation_rate_comparison.pdf"), width = 8, height = 5)
ggsave(file.path(figure_dir, "sweep_mutation_rate_comparison.png"), width = 8, height = 5, dpi = 300)
