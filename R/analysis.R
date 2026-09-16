# Basic reproduction number for the spatial SIS model.
# Sums per-contact transmission rates across spatial scales,
# weighted by number of contacts at each level, divided by recovery rate.
r_0 <- function(theta) {
  b  <- theta$beta[1:3]
  ss <- theta$spatial_sizes
  n  <- c(ss[1] - 1,      # same room, excluding self
          ss[2] - ss[1],  # same ward, different room
          ss[3] - ss[2])  # different ward
  sum(b * n) / theta$P$Dis[[1L]]
}

# Compute prevalence statistics from a simulate_outbreak() result.
#
# Returns a list with:
#   $series         data.frame of I, S, N, prevalence, I_to_S per time step
#   $mean_prevalence  time-averaged I/N
#   $mean_I_to_S      time-averaged I/S ratio
#   $sis_endemic    list with R0, endemic_prev, endemic_I_to_S from the
#                   deterministic SIS equilibrium (1 - 1/R0)
outbreak_prevalence <- function(OutbreakData) {
  theta  <- OutbreakData$theta
  State  <- OutbreakData$State
  n_beds <- theta$spatial_sizes[3]   # spatial_sizes[3] = NumBeds (total beds)

  # Restrict to bed rows only (exclude room-aggregate nodes)
  State <- State[seq_len(n_beds), , drop = FALSE]

  N <- nrow(State)
  I <- as.numeric(colSums(State))
  S <- N - I

  series <- data.frame(
    t          = seq_along(I) - 1L,
    I          = I,
    S          = S,
    N          = N,
    prevalence = I / N,
    I_to_S     = ifelse(S > 0, I / S, NA_real_)
  )

  R0           <- r_0(theta)
  endemic_prev <- max(0, 1 - 1 / R0)

  list(
    series          = series,
    mean_prevalence = mean(series$prevalence),
    mean_I_to_S     = mean(series$I_to_S, na.rm = TRUE),
    sis_endemic     = list(
      R0             = R0,
      endemic_prev   = endemic_prev,
      endemic_I_to_S = if (endemic_prev < 1) endemic_prev / (1 - endemic_prev) else Inf
    )
  )
}

# Build room-aware ground-truth ancestry for scoring
# mcmc(inference_mode = "patients_and_rooms") reconstructions, in the same
# combined node ordering the MCMC trace uses (ObsRec rows, then room_source
# rows — see mcmc_setup_rooms()). Unlike outbreak$ObsRec$Anc2 (which always
# skips rooms), this lets the true ancestor be a room, matching what the
# room-aware MCMC's own ancestor field can point to.
# use_room_tests must match the setting the MCMC was actually run with
# (selects outbreak$*$Anc3 vs. outbreak$*$Anc3_all as the source truth — see
# sim_postprocess() in simulation.R).
#
# Returns a named list with three parallel vectors, all aligned one-to-one
# with the columns of fit$anc from mcmc(inference_mode = "patients_and_rooms"):
#   $true_anc   — integer vector of true ancestor indices (NA = community import)
#   $adm_times  — lower bound on each node's infection/contamination time,
#                 using the *observable* Adm column (previous cleaning for rooms,
#                 admission day for patients). Do NOT substitute CaseRoomRec$Infc
#                 here: that is the latent ground-truth contamination time, which
#                 the MCMC does not observe and which would narrow the eligible-
#                 ancestor set relative to what the sampler actually saw.
#   $ptest_times — upper bound (positive-test day for patients/rooms with a swab,
#                  cleaning day for rooms without one — matching harmonise_room_rows).
#
# Pass $true_anc, $adm_times, $ptest_times directly to ancestry_score().
build_room_aware_truth <- function(outbreak, use_room_tests = TRUE) {
  CaseRec     <- outbreak$CaseRec
  CaseRoomRec <- outbreak$CaseRoomRec
  ObsRec      <- outbreak$ObsRec
  room_source <- if (use_room_tests) outbreak$ObsRoomRec else CaseRoomRec
  anc_col     <- if (use_room_tests) "Anc3" else "Anc3_all"

  n_pat_obs   <- nrow(ObsRec)
  n_room_case <- nrow(CaseRoomRec)

  # Position of each CaseRec/CaseRoomRec row within the *actual* combined MCMC
  # ordering (NA if that case is not part of the observed/used set).
  pat_pos      <- match(CaseRec$Id, ObsRec$Id)
  room_pos     <- if (use_room_tests) match(CaseRoomRec$Id, room_source$Id) else seq_len(n_room_case)
  combined_pos <- c(pat_pos, n_pat_obs + room_pos)

  # Raw truth (CaseRec/CaseRoomRec position space) for each row actually used
  # by the MCMC, in its own (ObsRec ++ room_source) order.
  raw_truth <- c(ObsRec[[anc_col]], room_source[[anc_col]])
  true_anc  <- combined_pos[raw_truth]

  # adm_times: use the observable Adm column (= previous cleaning time for rooms,
  # admission day for patients) — NOT CaseRoomRec$Infc (latent contamination time).
  # Using Adm makes the eligible-ancestor count consistent with what the MCMC saw.
  room_ptest <- ifelse(!is.na(room_source$PTest), room_source$PTest, room_source$Clean)
  adm_times  <- c(ObsRec$Adm,   room_source$Adm)
  ptest_times <- c(ObsRec$PTest, room_ptest)

  list(true_anc = true_anc, adm_times = adm_times, ptest_times = ptest_times)
}

# Plot simulated prevalence over time against the SIS endemic equilibrium.
plot_prevalence <- function(res) {
  library(ggplot2)

  sis_prev <- res$sis_endemic$endemic_prev
  sim_prev <- res$mean_prevalence
  T        <- max(res$series$t)

  ggplot(res$series, aes(x = t, y = prevalence)) +
    geom_line(colour = "#2c7bb6", linewidth = 0.5) +
    geom_hline(yintercept = sis_prev, colour = "#d7191c",
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = sim_prev, colour = "#1a9641",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = T * 0.7, y = sis_prev + 0.03,
             label = sprintf("SIS endemic = %.0f%%", sis_prev * 100),
             colour = "#d7191c", size = 3.5) +
    annotate("text", x = T * 0.7, y = sim_prev - 0.03,
             label = sprintf("Sim. mean = %.0f%%", sim_prev * 100),
             colour = "#1a9641", size = 3.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      x     = "Day",
      y     = "Prevalence (I / N)",
      title = sprintf("R\u2080 = %.1f \u2014 simulated prevalence vs. SIS endemic",
                      res$sis_endemic$R0)
    )
}
