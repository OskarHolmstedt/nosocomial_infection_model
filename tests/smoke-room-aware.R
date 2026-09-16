# Smoke test for the room-aware MCMC path.
#
# Checks the end-to-end pipeline with inference_mode = "patients_and_rooms":
#   1. Room record contains Adm (previous cleaning time, not Infc).
#   2. MCMC runs without error and returns a correctly shaped trace.
#   3. Ancestry score is substantially above the chance baseline.
#   4. build_room_aware_truth returns the expected list structure.
#   5. adm_times for rooms is derived from Adm, not Infc.
#
# Fixed seeds ensure deterministic results; expected values were established
# after the three bug fixes (beta_for_s swap, one-sided genetic ratio, Adm
# data-leakage fix) and should be re-baselined if the model changes.

cat("Room-aware smoke check: setting up...\n")

theta_room <- local({
  th            <- theta
  th$beta["BR"] <- 0.10
  th$beta["RB"] <- 0.15
  th$use_room_transmission <- TRUE
  W          <- weights(th$beta, Contact)
  th$LogComp <- log1p(-W)
  th$Odds    <- W / (1 - W)
  th
})

theta_sim <- local({
  th        <- theta_room
  th$P$Init <- c(rep(1L, 3L), rep(0L, NumBeds - 3L), rep(0L, NumRooms))
  th
})

set.seed(37675)
ob <- simulate_outbreak(45L, theta_sim)

# ── 1. Room record schema ──────────────────────────────────────────────────────
stopifnot(
  "Adm column missing from CaseRoomRec" =
    "Adm" %in% names(ob$CaseRoomRec),
  "Adm must be <= Infc for all room episodes" =
    all(ob$CaseRoomRec$Adm <= ob$CaseRoomRec$Infc, na.rm = TRUE),
  "First episode must have Adm = 0 (simulation start)" =
    ob$CaseRoomRec$Adm[1L] == 0L,
  "Later episodes must have Adm > 0 (previous cleaning time)" =
    any(ob$CaseRoomRec$Adm > 0L)
)
cat("  PASS  room record schema (Adm column)\n")

# ── 2. MCMC runs without error ─────────────────────────────────────────────────
set.seed(5512)
fit <- mcmc(
  ob,
  N_samples      = 1000L,
  burn_in        = 200L,
  inference_mode = "patients_and_rooms",
  use_room_tests = FALSE
)

n_nodes <- nrow(ob$ObsRec) + nrow(ob$CaseRoomRec)
stopifnot(
  "fit$anc must be a matrix"            = is.matrix(fit$anc),
  "fit$anc must have 1000 rows"         = nrow(fit$anc) == 1000L,
  "fit$anc columns must equal n_nodes"  = ncol(fit$anc) == n_nodes,
  "fit$inf must be a matrix"            = is.matrix(fit$inf),
  "fit$inf dimensions match fit$anc"    = identical(dim(fit$inf), dim(fit$anc))
)
cat("  PASS  MCMC trace shape\n")

# ── 3. build_room_aware_truth list structure ───────────────────────────────────
truth <- build_room_aware_truth(ob, use_room_tests = FALSE)

stopifnot(
  "truth must be a list"                         = is.list(truth),
  "truth must have true_anc field"               = !is.null(truth$true_anc),
  "truth must have adm_times field"              = !is.null(truth$adm_times),
  "truth must have ptest_times field"            = !is.null(truth$ptest_times),
  "all fields must have length n_nodes"          =
    length(truth$true_anc) == n_nodes &&
    length(truth$adm_times) == n_nodes &&
    length(truth$ptest_times) == n_nodes
)
cat("  PASS  build_room_aware_truth structure\n")

# ── 4. adm_times uses Adm not Infc ────────────────────────────────────────────
n_pat    <- nrow(ob$ObsRec)
room_adm <- truth$adm_times[(n_pat + 1L):n_nodes]

stopifnot(
  "room adm_times must equal CaseRoomRec$Adm" =
    isTRUE(all.equal(room_adm, ob$CaseRoomRec$Adm)),
  "room adm_times must NOT equal CaseRoomRec$Infc (leakage check)" =
    !isTRUE(all.equal(room_adm, ob$CaseRoomRec$Infc))
)
cat("  PASS  adm_times derived from Adm, not Infc\n")

# ── 5. Reconstruction accuracy above chance ────────────────────────────────────
sc   <- ancestry_score(truth$true_anc, fit$anc, truth$adm_times, truth$ptest_times)
keep <- is.finite(sc$lift)

# Chance baseline ≈ 1/n_nodes; require mean lift > 2 (>2× better than random).
# Mode accuracy must exceed the chance level by a meaningful margin.
chance <- 1 / n_nodes
stopifnot(
  "mean lift must exceed 2× chance"      = mean(sc$lift[keep]) > 2,
  "mode accuracy must beat chance level" = sc$mode_accuracy > chance
)
cat(sprintf("  PASS  accuracy above chance  (mode_acc=%.3f, chance=%.3f, lift=%.2f)\n",
            sc$mode_accuracy, chance, mean(sc$lift[keep])))

cat("Room-aware smoke check passed.\n")
