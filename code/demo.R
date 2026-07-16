source("two_wards_data.R")
source("prep_functions.R")
source("outbreak_simulation.R")
source("mcmc.R")
source("analysis_functions.R")
source("parameters.R")
source("plot_functions.R")

T_demo     <- 45
theta_demo <- local({
  th        <- theta
  th$P$Init <- c(1, rep(0, NumBeds - 1), rep(0, NumRooms))
  th
})

# ── Simulate demo outbreak (seed 7104: 18 obs / 40 total cases) ───────────────
set.seed(7104)
OutbreakDemo <- simulate_outbreak(T_demo, theta_demo)

# ── Run MCMC ──────────────────────────────────────────────────────────────────
set.seed(2413)
Y_demo <- mcmc(OutbreakDemo, N_samples = 5000)

score <- ancestry_score(OutbreakDemo$ObsRec$Anc2, Y_demo$anc,
                        OutbreakDemo$ObsRec$Adm, OutbreakDemo$ObsRec$PTest)
print(score)

# ── Plots ─────────────────────────────────────────────────────────────────────
plot_timeline(OutbreakDemo, observed = FALSE, show_all_stays = TRUE)
plot_timeline(OutbreakDemo, observed = TRUE,  show_all_stays = TRUE)

plot_timeline_posterior(OutbreakDemo, Y_demo, min_prob = 0.2)

plot_timeline_mode(OutbreakDemo, Y_demo, color_edges = FALSE, color_nodes = FALSE)
plot_timeline_mode(OutbreakDemo, Y_demo, color_edges = TRUE,  color_nodes = FALSE)

# ── MCMC convergence diagnostics ─────────────────────────────────────────────
plot_mcmc_diagnostics(Y_demo, OutbreakDemo)
plot_mode_vs_truth(Y_demo, OutbreakDemo)

# ── Save poster figures ───────────────────────────────────────────────────────
p_input <- plot_timeline(OutbreakDemo,
                         uniform_edge_color = "#555555", edge_width_by = "fixed")
save_timeline(p_input, "input_data.png")
save_timeline(p_input, "input_data.pdf")

p_mode <- plot_timeline_mode(OutbreakDemo, Y_demo,
                             color_edges = TRUE, color_nodes = FALSE)
save_timeline(p_mode, "demo_mode.png")
save_timeline(p_mode, "demo_mode.pdf")
