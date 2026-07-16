source(here::here("scripts", "_load_project.R"))

read_positive_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  value
}

n_samples <- read_positive_integer("NOSOCOMIAL_QUICK_START_SAMPLES", 1000L)
burn_in <- read_positive_integer("NOSOCOMIAL_QUICK_START_BURN_IN", 200L)
skip_visualization <- identical(
  tolower(Sys.getenv("NOSOCOMIAL_SKIP_VISUALIZATION", "false")),
  "true"
)

theta_quick <- local({
  th <- theta
  th$P$Init <- c(1, rep(0, NumBeds - 1), rep(0, NumRooms))
  th
})

# Fixed seeds and explicit RNG settings make the core example reproducible.
set.seed(
  30,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
outbreak_quick <- simulate_outbreak(T = 45, theta = theta_quick)

set.seed(
  2413,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
fit_quick <- mcmc(
  outbreak_quick,
  N_samples = n_samples,
  burn_in = burn_in
)

score_quick <- ancestry_score(
  outbreak_quick$ObsRec$Anc2,
  fit_quick$anc,
  outbreak_quick$ObsRec$Adm,
  outbreak_quick$ObsRec$PTest
)

quick_start_summary <- list(
  total_cases = nrow(outbreak_quick$CaseRec),
  observed_cases = nrow(outbreak_quick$ObsRec),
  mcmc_samples = n_samples,
  burn_in = burn_in,
  mode_accuracy = score_quick$mode_accuracy,
  mean_p_true = score_quick$mean_p_true
)

cat(
  "Quick start complete\n",
  sprintf("  Total simulated cases: %d\n", quick_start_summary$total_cases),
  sprintf("  Observed cases: %d\n", quick_start_summary$observed_cases),
  sprintf("  MCMC samples: %d (+ %d burn-in)\n", n_samples, burn_in),
  sprintf("  Ancestry mode accuracy: %.3f\n", quick_start_summary$mode_accuracy),
  sprintf("  Mean P(true ancestor): %.3f\n", quick_start_summary$mean_p_true),
  sep = ""
)

if (!skip_visualization) {
  output_file <- Sys.getenv(
    "NOSOCOMIAL_QUICK_START_OUTPUT",
    here::here("output", "quick-start-reconstruction.html")
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  reconstruction <- plot_timeline_mode(
    outbreak_quick,
    fit_quick,
    color_edges = TRUE,
    color_nodes = FALSE
  )
  htmlwidgets::saveWidget(
    reconstruction,
    output_file,
    selfcontained = FALSE,
    title = "Nosocomial outbreak reconstruction"
  )
  cat("  Visualization: ", normalizePath(output_file), "\n", sep = "")
}
