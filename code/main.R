library(Matrix)
library(igraph)
library(epicontacts)
library(visNetwork)
library(expm)
source("plot_functions.R")
source("prep_functions.R")
list2env(hospital(Wards = 3, RoomsPerWard = 4, BedsPerRoom = 3), envir = .GlobalEnv)
source("parameters.R")
source("outbreak_simulation.R")
source("mcmc.R")
source("analysis_functions.R")
#set.seed(51)
T = 50
OutbreakData = simulate_outbreak(T, theta)
plot_timeline(OutbreakData, observed = FALSE, show_all_stays = TRUE)
plot_timeline(OutbreakData, observed = TRUE, show_all_stays = TRUE)
plot_timeline(OutbreakData, show_all_stays = TRUE)
plot_timeline(OutbreakData, show_tree = FALSE, show_all_stays = TRUE)
plot_timeline(OutbreakData, show_tree = FALSE, observed = FALSE, show_all_stays = TRUE)
plot_timeline(OutbreakData, show_tree = FALSE, observed = TRUE, show_all_stays = FALSE)
plot_timeline(OutbreakData, flip = TRUE)
plot_timeline(OutbreakData, edge_width_by = "gen")

Y = mcmc(OutbreakData, N_samples = 1000)
Z = mcmc(OutbreakData, N_samples = 10000)

#posterior_mode_anc(Y$anc)
#OutbreakData$ObsRec$Anc2
x = OutbreakData$ObsRec$Anc2
#y = posterior_mode_anc(Y$anc)
#res <- (x == y) | (is.na(x) & is.na(y))
#res <- (!is.na(x) & !is.na(y) & x == y) | (is.na(x) & is.na(y))
#mean(res)

ancestry_score(x, Y$anc, OutbreakData$ObsRec$Adm, OutbreakData$ObsRec$PTest)
ancestry_score(x, Z$anc, OutbreakData$ObsRec$Adm, OutbreakData$ObsRec$PTest)
plot_timeline_posterior(OutbreakData, Y, min_prob = 0.3)
plot_timeline_posterior(OutbreakData, Z, min_prob = 0.1)
plot_timeline_mode(OutbreakData, Y)
plot_timeline_mode(OutbreakData, Z)


# --- Endemic analysis ---
# Low-R0 variant: scale all betas proportionally to R0 = 2.5,
# keeping hospital-wide spread non-zero. Start with 5 infected to
# avoid early stochastic extinction at this lower R0.
theta_low           <- theta
theta_low$beta[1:3] <- theta$beta[1:3] * (2.5 / r_0(theta))
theta_low$P$Init    <- c(rep(1, 5), rep(0, NumBeds - 5), rep(0, NumRooms))
W                   <- weights(theta_low$beta, Contact)
theta_low$LogComp   <- log1p(-W)
theta_low$Odds      <- W / (1 - W)

T_endemic      <- 500
OD_endemic     <- simulate_outbreak(T_endemic, theta)
OD_endemic_low <- simulate_outbreak(T_endemic, theta_low)

res_endemic     <- outbreak_prevalence(OutbreakData)
res_endemic_low <- outbreak_prevalence(OD_endemic_low)

plot_prevalence(res_endemic)
plot_prevalence(res_endemic_low)

library(htmlwidgets)
library(webshot2)

p <- plot_timeline_mode(OutbreakData, Y)
saveWidget(p, tmp <- tempfile(fileext = ".html"), selfcontained = TRUE)
webshot(tmp, "p9.png", vwidth = 1400, vheight = 900)
