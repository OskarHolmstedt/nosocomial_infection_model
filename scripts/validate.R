source(here::here("scripts", "_load_project.R"))

required_functions <- c(
  "hospital",
  "simulate_outbreak",
  "mcmc",
  "ancestry_score",
  "plot_timeline_mode"
)
missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]
if (length(missing_functions)) {
  stop(
    "Model source check failed; missing: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}
cat("Model source check passed.\n")

source(here::here("tests", "testthat.R"))

smoke_environment <- new.env(parent = globalenv())
sys.source(
  here::here("tests", "smoke-quick-start.R"),
  envir = smoke_environment
)

skip_reports <- identical(
  tolower(Sys.getenv("NOSOCOMIAL_SKIP_REPORT_CHECKS", "false")),
  "true"
)

if (!skip_reports) {
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) {
    stop("Quarto is required for report validation.", call. = FALSE)
  }

  report_output <- tempfile(
    pattern = "report-validation-",
    tmpdir = here::here("tmp")
  )
  dir.create(report_output, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(report_output, recursive = TRUE, force = TRUE), add = TRUE)

  reports <- here::here(
    "reports",
    c("model_overview.qmd", "testing_rate_sweep.qmd")
  )
  for (report in reports) {
    report_name <- tools::file_path_sans_ext(basename(report))
    report_directory <- file.path(report_output, report_name)
    dir.create(report_directory, recursive = TRUE, showWarnings = FALSE)
    cat("Rendering ", basename(report), " without execution...\n", sep = "")
    status <- system2(
      quarto,
      c(
        "render",
        shQuote(report),
        "--no-execute",
        "--output-dir",
        shQuote(report_directory)
      )
    )
    if (!identical(status, 0L)) {
      stop("Quarto validation failed for ", basename(report), call. = FALSE)
    }
  }
  cat("Quarto report checks passed.\n")
}

cat("All validation checks passed.\n")
