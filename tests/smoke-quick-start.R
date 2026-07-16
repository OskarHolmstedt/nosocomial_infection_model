Sys.setenv(
  NOSOCOMIAL_QUICK_START_SAMPLES = "1000",
  NOSOCOMIAL_QUICK_START_BURN_IN = "200",
  NOSOCOMIAL_SKIP_VISUALIZATION = "true"
)

source(here::here("scripts", "quick_start.R"))

stopifnot(
  identical(quick_start_summary$total_cases, 17L),
  identical(quick_start_summary$observed_cases, 10L),
  isTRUE(all.equal(quick_start_summary$mode_accuracy, 0.4)),
  isTRUE(all.equal(quick_start_summary$mean_p_true, 0.4629, tolerance = 1e-4))
)

cat("Quick-start smoke check passed.\n")
