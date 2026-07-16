if (!exists("hospital", mode = "function")) {
  source(here::here("R", "hospital.R"))
}

# Flat hospital: same number of beds as the global hospital, all in one room.
# Beta is calibrated so the total per-step force of infection matches — i.e.
# R0 is preserved. All spatial hierarchy is removed; S=1 for every bed pair.
#
# Usage: source this file to get theta_flat (and flat_hospital for diagnostics).
# Does NOT overwrite the global hospital variables (NumBeds, Contact, etc.).

flat_hospital <- local({
  # ---- Read actual global hospital structure ----------------------------------
  # Compute total per-step force of infection from bed 1 in the live W matrix
  W_ref      <- weights(theta_demo$beta, Contact)
  total_rate <- sum(W_ref[1L, seq_len(NumBeds)])  # bed 1 → all other beds
  n_beds     <- NumBeds

  # ---- Flat hospital layout (one room, all n_beds beds) ----------------------
  hosp          <- hospital(Wards = 1L, RoomsPerWard = 1L, BedsPerRoom = n_beds)
  BrB_flat      <- total_rate / (n_beds - 1L)
  beta          <- c(BrB = BrB_flat, BwB = 0, BhB = 0, BR = 0, RB = 0, RR = 0)

  # ---- Theta (mirrors parameters.R structure) --------------------------------
  P <- list(
    Init = c(1, rep(0, n_beds - 1), rep(0, hosp$NumRooms)),
    Imp  = rep(0L,  n_beds),
    Test = rep(0.1, n_beds),
    Dis  = rep(0.1, n_beds)
  )
  W     <- weights(beta, hosp$Contact)
  theta <- list(
    mu            = theta_demo$mu,
    beta          = beta,
    P             = P,
    lambda_import = theta_demo$lambda_import,
    pi_import     = theta_demo$pi_import,
    LogComp       = log1p(-W),
    Odds          = W / (1 - W),
    spatial_sizes = hosp$spatial_sizes,
    Positions     = hosp$Positions
  )

  list(NumBeds       = n_beds,
       NumRooms      = hosp$NumRooms,
       Contact       = hosp$Contact,
       Positions     = hosp$Positions,
       spatial_sizes = hosp$spatial_sizes,
       theta         = theta,
       BrB_flat      = BrB_flat,
       total_rate    = total_rate)
})

theta_flat <- flat_hospital$theta
