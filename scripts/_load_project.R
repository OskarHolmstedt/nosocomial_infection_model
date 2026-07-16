suppressPackageStartupMessages({
  library(here)
  library(Matrix)
  library(igraph)
  library(epicontacts)
  library(visNetwork)
  library(expm)
  library(ggplot2)
})

source(here::here("R", "hospital.R"))
source(here::here("R", "simulation.R"))
source(here::here("R", "inference.R"))
source(here::here("R", "analysis.R"))
source(here::here("R", "plotting.R"))

# Default research configuration used by the demo, reports, and sweeps.
source(here::here("config", "hospital_two_wards.R"))
source(here::here("config", "parameters.R"))
