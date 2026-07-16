repositories <- c(
  reconhub = "https://reconhub.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
)

if (!requireNamespace("epicontacts", quietly = TRUE)) {
  install.packages("epicontacts", repos = repositories)
}

required_packages <- c(
  "Matrix",
  "igraph",
  "epicontacts",
  "visNetwork",
  "expm",
  "here",
  "htmlwidgets",
  "tidyverse",
  "patchwork",
  "ggrepel",
  "knitr",
  "rmarkdown",
  "testthat"
)

optional_packages <- c(
  "webshot2",
  "chromote",
  "magick",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R packages after environment setup: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

missing_optional_packages <- optional_packages[
  !vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_optional_packages)) {
  warning(
    "Optional high-resolution export packages are unavailable: ",
    paste(missing_optional_packages, collapse = ", "),
    call. = FALSE
  )
}

cat("All required R dependencies are available.\n")
