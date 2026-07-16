smoke_output <- here::here("tmp", "quick-start-smoke.html")
smoke_dependencies <- sub("\\.html$", "_files", smoke_output)
unlink(c(smoke_output, smoke_dependencies), recursive = TRUE, force = TRUE)

Sys.setenv(
  NOSOCOMIAL_QUICK_START_SAMPLES = "1000",
  NOSOCOMIAL_QUICK_START_BURN_IN = "200",
  NOSOCOMIAL_QUICK_START_OUTPUT = smoke_output,
  NOSOCOMIAL_SKIP_VISUALIZATION = "false"
)

source(here::here("scripts", "quick_start.R"), local = TRUE)

stopifnot(
  identical(quick_start_summary$total_cases, 17L),
  identical(quick_start_summary$observed_cases, 10L),
  isTRUE(all.equal(quick_start_summary$mode_accuracy, 0.4)),
  isTRUE(all.equal(quick_start_summary$mean_p_true, 0.4629, tolerance = 1e-4)),
  file.exists(smoke_output)
)

unlink(c(smoke_output, smoke_dependencies), recursive = TRUE, force = TRUE)
cat("Quick-start smoke check passed.\n")
