testthat::test_dir(
  here::here("tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)
