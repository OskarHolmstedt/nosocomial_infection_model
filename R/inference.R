## Run MCMC outbreak reconstruction.
# inference_mode: "patients"            — infer over detected patients only (default)
#                 "patients_and_rooms"  — include room contamination episodes as
#                                         potential ancestors; requires CaseRoomRec
#                                         in outbreak_data
# use_room_tests: when TRUE and inference_mode = "patients_and_rooms", room positive
#                 swabs constrain room infection times (uses theta$P$RoomTest).
#                 When FALSE, room infection times are fully latent.
mcmc <- function(outbreak_data, N_samples, burn_in = 0,
                 use_genetics   = TRUE,
                 inference_mode = c("patients", "patients_and_rooms"),
                 use_room_tests = TRUE) {
  inference_mode <- match.arg(inference_mode)
  mcmc_item <- mcmc_setup(outbreak_data, N_samples, burn_in,
                          use_genetics, inference_mode, use_room_tests)
  mcmc_item <- mcmc_loop(mcmc_item)
  mcmc_cleanup(mcmc_item)
}
mcmc_setup <- function(outbreak_data, N_samples, burn_in,
                       use_genetics   = TRUE,
                       inference_mode = "patients",
                       use_room_tests = TRUE) {
  theta <- outbreak_data$theta

  # P$Test and P$Dis are per-bed vectors in the simulation; MCMC needs scalars
  theta$P$Test  <- theta$P$Test[[1L]]
  theta$P$Dis   <- theta$P$Dis[[1L]]
  theta$P$Clean    <- theta$P$Clean[[1L]]
  theta$P$RoomTest <- theta$P$RoomTest[[1L]]

  if (inference_mode == "patients_and_rooms") {
    return(mcmc_setup_rooms(outbreak_data, theta, N_samples, burn_in,
                            use_genetics, use_room_tests))
  }

  ## Patient-only mode (original behaviour) ─────────────────────────────────
  Rec    <- outbreak_data$ObsRec
  D      <- outbreak_data$ObsDist

  Rec    <- init(Rec)
  S      <- create_s_matrix(Rec)
  C      <- create_c_matrix(Rec)
  t_root <- min(Rec$Adm) - 1L
  Tms    <- create_transition_matrices(theta, kap_max = max(Rec$PTest) - min(Rec$Adm) + 1L)

  n_cases <- nrow(Rec)
  Trace   <- list(inf = matrix(NA_integer_, N_samples, n_cases),
                  anc = matrix(NA_integer_, N_samples, n_cases),
                  kap = matrix(NA_integer_, N_samples, n_cases))
  list(Rec            = Rec,
       theta          = theta,
       Dist           = D,
       S              = S,
       C              = C,
       t_root         = t_root,
       Tms            = Tms,
       Trace          = Trace,
       N_samples      = N_samples,
       burn_in        = burn_in,
       use_genetics   = use_genetics,
       inference_mode = inference_mode,
       use_room_tests = use_room_tests)
}
mcmc_loop <- function(mcmc_item) {
  list2env(mcmc_item, envir = environment())

  for (n in seq_len(burn_in + N_samples)) {
    Rec = propose_infection_times(Rec, theta, S, t_root)
    Rec = propose_ancestors(Rec, Dist, theta, S, C, Tms$TmPowF, t_root, use_genetics)
    Rec = propose_kappa(Rec, Dist, theta, S, C, Tms$TmPowF, t_root, use_genetics)
    Rec = propose_subtree_shift(Rec, theta, S, t_root)
    if (n > burn_in) {
      i <- n - burn_in
      Trace$inf[i, ] <- Rec$inf
      Trace$anc[i, ] <- Rec$anc
      Trace$kap[i, ] <- Rec$kap
    }
  }
  mcmc_item$Rec   <- Rec
  mcmc_item$Trace <- Trace
  mcmc_item
}
mcmc_cleanup <- function(mcmc_item) {
  mcmc_item$Trace
}
## Room-aware MCMC setup ────────────────────────────────────────────────────
# Builds a combined patient + room record and extended S/C matrices.
mcmc_setup_rooms <- function(outbreak_data, theta, N_samples, burn_in,
                             use_genetics, use_room_tests) {
  pat_rec <- outbreak_data$ObsRec
  D_pat   <- outbreak_data$ObsDist

  # Choose room source: observed (positive swab) or all contaminated episodes
  room_source <- if (use_room_tests) outbreak_data$ObsRoomRec
                 else                outbreak_data$CaseRoomRec

  # Shared schema for combined record: keep only the columns MCMC needs,
  # plus node_type / has_ptest to distinguish patients from rooms in proposals.
  mcmc_cols <- c("Id","Bed","Room","Ward","Adm","Dis","Infc","NTest","PTest","Anc","Mut")
  pat_rows  <- data.frame(pat_rec[, mcmc_cols, drop = FALSE],
                           node_type = "patient",
                           has_ptest = TRUE,
                           stringsAsFactors = FALSE)

  # Combine patients and rooms into one record with shared schema
  Rec <- if (!is.null(room_source) && nrow(room_source) > 0)
    rbind(pat_rows, harmonise_room_rows(room_source, nrow(pat_rec)))
  else
    pat_rows

  n_pat   <- nrow(pat_rec)
  n_rooms <- nrow(Rec) - n_pat

  # Extend genetic distance matrix: rooms have Inf distance to all others
  # (no genetic information from environmental samples)
  n_total <- nrow(Rec)
  D <- matrix(Inf, n_total, n_total)
  D[seq_len(n_pat), seq_len(n_pat)] <- D_pat

  # Initialise latent variables on the combined record.
  # init_inf uses PTest as upper bound; rooms without a real swab have PTest set
  # to Clean (Dis) by harmonise_room_rows, so this works for all rows.
  Rec    <- init(Rec)
  S      <- create_s_matrix_extended(Rec)
  C      <- create_c_matrix_extended(Rec)
  t_root <- min(Rec$Adm, na.rm = TRUE) - 1L
  Tms    <- create_transition_matrices_extended(
              theta,
              kap_max = max(Rec$PTest, na.rm = TRUE) - min(Rec$Adm, na.rm = TRUE) + 1L)

  Trace <- list(inf = matrix(NA_integer_, N_samples, n_total),
                anc = matrix(NA_integer_, N_samples, n_total),
                kap = matrix(NA_integer_, N_samples, n_total))

  list(Rec            = Rec,
       theta          = theta,
       Dist           = D,
       S              = S,
       C              = C,
       t_root         = t_root,
       Tms            = Tms,
       Trace          = Trace,
       N_samples      = N_samples,
       burn_in        = burn_in,
       use_genetics   = use_genetics,
       inference_mode = "patients_and_rooms",
       use_room_tests = use_room_tests,
       n_pat          = n_pat,
       n_rooms        = n_rooms)
}

## Map room contamination episode columns to the patient-record schema.
# Rooms have no Bed.
#   Adm  = room_rec$Adm   — time room last became susceptible (previous cleaning,
#                            or 0 for the first episode). This is the observable
#                            lower bound on contamination time; using Infc here
#                            would leak the latent ground-truth to the sampler.
#   Dis  = room_rec$Clean — current decontamination time (observable: staff clean).
#   Infc = NA             — contamination time is latent; the sampler infers it
#                            in [Adm, min(PTest, Dis)].
# node_type = "room" is used in proposals to apply correct rates.
# has_ptest = TRUE means there is a real observed swab; FALSE means PTest is set
# to Dis as a bounding proxy (and the test-rate likelihood is suppressed).
harmonise_room_rows <- function(room_rec, n_pat_rows) {
  n      <- nrow(room_rec)
  has_pt <- !is.na(room_rec$PTest)
  PTest  <- ifelse(has_pt, room_rec$PTest, room_rec$Clean)
  data.frame(
    Id        = as.character(n_pat_rows + seq_len(n)),
    Bed       = NA_character_,
    Room      = room_rec$Room,
    Ward      = room_rec$Ward,
    Adm       = room_rec$Adm,    # previous cleaning time -- observable lower bound
    Dis       = room_rec$Clean,  # current cleaning time  -- observable upper bound
    Infc      = NA_integer_,     # latent: inferred by the sampler
    NTest     = room_rec$NTest,
    PTest     = PTest,
    Anc       = NA_character_,
    Mut       = NA_integer_,
    node_type = "room",
    has_ptest = has_pt,
    stringsAsFactors = FALSE
  )
}

#########
# Setup #
#########
create_s_matrix <- function(Rec) {
  same_room <- outer(Rec$Room, Rec$Room, "==")
  same_ward <- outer(Rec$Ward, Rec$Ward, "==")
  3L - same_ward - same_room
}
create_c_matrix <- function(Rec) {
  pmax(outer(Rec$Dis, Rec$Dis, pmin) -
         outer(Rec$Adm, Rec$Adm, pmax) + 1L, 0L)
}
create_transition_matrices <- function(theta, kap_max = 10L) {
  b  <- theta$beta[1:3]
  Tm <- outer(1:3, 1:3, function(i, j) b[pmax(i, j)]) %*% diag(theta$spatial_sizes)

  # Cache all powers Tm^1 ... Tm^kap_max as a 3x3xkap_max tensor
  TmPow        <- array(0, dim = c(3L, 3L, kap_max))
  TmPow[,,1L]  <- Tm
  for (k in seq_len(kap_max - 1L) + 1L)
    TmPow[,,k] <- TmPow[,,k - 1L] %*% Tm

  # Lookup table: TmPowF[k, s] = (Tm^k %*% f)[s], replaces T %^% (kap-1) %*% f at runtime
  TmPowF <- t(apply(TmPow, 3L, function(M) drop(M %*% b)))

  list(TmPow = TmPow, TmPowF = TmPowF)
}

## Extended spatial matrix for patient + room nodes.
# S[i, j] encodes the directional spatial relationship from node j to node i:
#   1: patient <- patient, same room  (beta_BrB)
#   2: patient <- patient, same ward  (beta_BwB)
#   3: patient <- patient, hospital   (beta_BhB)
#   4: patient <- room,    same room  (beta_RB)
#   5: room    <- patient, same room  (beta_BR)
#   6: room    <- room,    same ward  (beta_RR)
create_s_matrix_extended <- function(Rec) {
  n         <- nrow(Rec)
  is_room   <- !is.na(Rec$node_type) & Rec$node_type == "room"
  is_pat    <- !is_room

  same_room <- outer(Rec$Room, Rec$Room, "==")
  same_ward <- outer(Rec$Ward, Rec$Ward, "==")

  S <- matrix(0L, n, n)

  # patient <- patient
  pp <- outer(which(is_pat), which(is_pat))
  S[pp] <- (3L - same_ward[pp] - same_room[pp])

  # patient <- room (same room only; cross-room contamination not modelled)
  pr_i <- which(is_pat)
  pr_j <- which(is_room)
  if (length(pr_i) > 0 && length(pr_j) > 0) {
    for (j in pr_j)
      S[pr_i[same_room[pr_i, j]], j] <- 4L
  }

  # room <- patient (same room only)
  rp_i <- which(is_room)
  rp_j <- which(is_pat)
  if (length(rp_i) > 0 && length(rp_j) > 0) {
    for (i in rp_i)
      S[i, rp_j[same_room[i, rp_j]]] <- 5L
  }

  # room <- room (same ward, excluding self)
  rr_i <- which(is_room)
  if (length(rr_i) > 1) {
    for (i in rr_i)
      S[i, rr_i[same_ward[i, rr_i] & rr_i != i]] <- 6L
  }

  S
}

## Extended contact matrix: days of overlap between node i and node j.
# For room nodes, Adm = contamination start and Dis = cleaning time.
create_c_matrix_extended <- function(Rec) {
  pmax(outer(Rec$Dis, Rec$Dis, pmin) -
         outer(Rec$Adm, Rec$Adm, pmax) + 1L, 0L)
}

## Extended 6×6 transition matrix for patient + room chains.
# Spatial codes 1-3 = patient-patient (uses existing 3x3 block).
# Codes 4-6 (room-involved) use direct beta scaling; the one-step transition
# into codes 4/5/6 from the 3x3 patient block is via off-diagonal extensions.
create_transition_matrices_extended <- function(theta, kap_max = 10L) {
  b   <- theta$beta  # length 6: BrB, BwB, BhB, BR, RB, RR
  sz  <- theta$spatial_sizes

  # 6x6 transition matrix (rows/cols = spatial codes 1..6)
  # Entry [i,j]: expected transmission from scale j to scale i in one step
  # Patient-patient block (1:3, 1:3): same as original
  Tm6 <- matrix(0, 6L, 6L)
  Tm6[1:3, 1:3] <- outer(1:3, 1:3, function(i, j) b[pmax(i, j)]) %*% diag(sz)
  Tm6[5L, 1:3]  <- b[4L]   # patient (any scale) -> room, same room
  Tm6[4L, 5L]   <- b[5L]   # room -> patient, same room
  Tm6[6L, 5L]   <- b[6L]   # room -> room (within ward)
  Tm6[6L, 6L]   <- b[6L]   # room -> room self (excluded in contact matrix but keep rate)

  # Codes 4 and 6 are "a patient" / "a room" that arrived via a room-mediated /
  # room-sourced edge respectively. As ONGOING SOURCES for further hidden hops,
  # they must behave identically to codes 1 and 5 (same underlying node type —
  # the specific edge that produced them doesn't change their own future
  # transmission rates). Without this, codes 4 and 6 are unreachable as
  # intermediates for kap > 1 chains (column 4 was otherwise all-zero, and
  # column 6 fed nothing into row 4), so any reconstructed ancestry with more
  # than one hidden generation, ending in a room-mediated edge, always got
  # probability exactly 0 regardless of how well timing/genetics supported it.
  Tm6[1:3, 4L] <- Tm6[1:3, 1L]  # room-infected patient -> onward patient-patient spread
  Tm6[5L,  4L] <- Tm6[5L,  1L]  # room-infected patient -> can still contaminate a room
  Tm6[4L,  6L] <- Tm6[4L,  5L]  # room-sourced room -> can still infect a patient (RB)
  # (row 6 col 6 already covers room-sourced room -> room; row 5 col 6 stays 0 —
  # "room contaminated by a patient" (row 5's own definition) doesn't apply when
  # the source is itself a room)

  TmPow        <- array(0, dim = c(6L, 6L, kap_max))
  TmPow[,,1L]  <- Tm6
  for (k in seq_len(kap_max - 1L) + 1L)
    TmPow[,,k] <- TmPow[,,k - 1L] %*% Tm6

  TmPowF <- t(apply(TmPow, 3L, function(M) drop(M %*% b)))

  list(TmPow = TmPow, TmPowF = TmPowF)
}

#########################
# Initial distributions #
#########################
init <- function(Rec) {
  Rec$inf <- init_inf(Rec)
  Rec$anc <- init_anc(Rec)
  Rec$kap <- init_kap(Rec)
  Rec
}
## Initial infection times
# Sampled uniformly between admission and discharge date
init_inf <- function(Rec) {
  as.integer(runif(nrow(Rec),
                   Rec$Adm,
                   Rec$PTest + 1))
}
## Initial ancestors
# Sample uniformly from cases with earlier infection time; NA = community case
init_anc <- function(Rec) {
  order_inf <- order(Rec$inf) # sorted positions → original indices
  n_eligible <- rank(Rec$inf, ties.method = "min") - 1L # how many eligible ancestors each case has
  sampled <- as.integer(runif(nrow(Rec), 0, n_eligible)) # vectorised uniform draw over {0,...,k-1}
  ifelse(n_eligible == 0L, NA_integer_, order_inf[sampled + 1L])
}
## Initial number of generations
# Sample uniformly between 1 and generation time
init_kap <- function(Rec) {
  anc_inf <- ifelse(is.na(Rec$anc),
                    0,
                    Rec$inf[Rec$anc])
  as.integer(runif(nrow(Rec), 0, Rec$inf - anc_inf)) + 1L
}


##########################
# Proposal distributions #
##########################
generate_inf_proposals <- function(Rec, theta) {
  n      <- nrow(Rec)
  delta  <- pmax(1L, as.integer(round((Rec$PTest - Rec$Adm) * 0.25)))
  NewInf <- reflect(Rec$inf + as.integer(round(runif(n, -1, 1) * delta)), Rec$Adm, Rec$PTest)
  Accept <- log(runif(n)) - ll_testing_inf_move(NewInf, Rec, theta)
  list(NewInf = NewInf, Accept = Accept)
}

## Proposal function for infection times
propose_infection_times <- function(Rec, theta, S, t_root) {
  n   <- nrow(Rec)
  idx <- seq_len(n)

  children_list <- get_children(Rec$anc)
  is_room       <- node_is_room(Rec)

  # Hoist: these only depend on fixed quantities (theta, S, Rec$anc, Rec$kap)
  # beta_anc only used for hospital-acquired cases (non-NA anc); NA-safe indexing via replace
  anc_idx      <- replace(Rec$anc, is.na(Rec$anc), 1L)
  beta_anc     <- beta_for_s(S[cbind(idx, anc_idx)], theta$beta)  # imported slots are never read
  kap_ch_list  <- lapply(idx, \(i) Rec$kap[children_list[[as.character(i)]]])
  # S[child, i] (not S[i, child]): S[a,b] encodes the relationship *from* b *to*
  # a, and for a child c of i, the edge runs from i to c, so the code we want is
  # S[c, i]. This block of S is symmetric for patient-patient pairs (which is
  # why it was invisible there), but not for edges touching a room — e.g. a
  # child that is a room needs the patient-> room rate (BR), not the room->
  # patient rate (RB) that S[i, child] would silently substitute.
  beta_ch_list <- lapply(idx, \(i) beta_for_s(S[children_list[[as.character(i)]], i], theta$beta))

  p      <- generate_inf_proposals(Rec, theta)
  NewInf <- p$NewInf
  Accept <- p$Accept

  for (i in sample.int(n)) {
    if (is.na(Rec$anc[i])) {
      own_ta <- own_tb <- own_kap <- own_beta <- NULL
      # Rooms cannot be community imports; community prior only applies to patients
      ImportLL <- if (is_room[i]) 0
                  else (NewInf[i] - Rec$inf[i]) * log(1 - theta$lambda_import)
    } else {
      anc_time <- Rec$inf[Rec$anc[i]]
      own_ta   <- Rec$inf[i] - anc_time
      own_tb   <- NewInf[i]  - anc_time
      own_kap  <- Rec$kap[i]
      own_beta <- beta_anc[i]
      ImportLL <- 0
    }

    children <- children_list[[as.character(i)]]

    TimeLL <- ll_time_inf_move_local(
      c(own_ta,   Rec$inf[children] - Rec$inf[i]),
      c(own_tb,   Rec$inf[children] - NewInf[i]),
      c(own_kap,  kap_ch_list[[i]]),
      c(own_beta, beta_ch_list[[i]])
    )

    if (Accept[i] < TimeLL + ImportLL)
      Rec$inf[i] <- NewInf[i]
  }
  Rec
}

propose_ancestors <- function(Rec, D, theta, S, C, TmPowF, t_root, use_genetics = TRUE) {
  n        <- nrow(Rec)
  idx      <- seq_len(n)
  is_room  <- node_is_room(Rec)

  # --- Propose new ancestors --------------------------------------------------
  # Eligible ancestors: cases infected at least kap steps earlier.
  # NA (community) is included as one of n_elig+1 candidates (index 0) so that
  # hospital ↔ community moves are possible for patient nodes.
  # Room nodes cannot be community cases: NA is excluded from their proposal.
  order_inf <- order(Rec$inf)
  n_elig    <- findInterval(Rec$inf - Rec$kap, Rec$inf[order_inf])
  # For rooms: sample only from eligible hospital ancestors (n_elig candidates, no NA).
  # For patients: n_elig + 1 candidates (including NA = community).
  n_cand  <- ifelse(is_room, n_elig, n_elig + 1L)
  sampled <- as.integer(runif(n, 0, pmax(n_cand, 1L)))  # pmax avoids runif(n, 0, 0)
  # Rooms: sampled in 0..(n_elig-1) → always pick a real ancestor
  # Patients: sampled == 0 → NA (community)
  NewAnc  <- ifelse(!is_room & sampled == 0L,
                    NA_integer_,
                    order_inf[pmax(ifelse(is_room, sampled + 1L, sampled), 1L)])
  # Rooms with no eligible ancestors: keep current ancestor (no move)
  NewAnc  <- ifelse(is_room & n_elig == 0L, Rec$anc, NewAnc)

  # Classify transition type for each case
  to_real   <- is.na(Rec$anc) & !is.na(NewAnc)   # community → hospital
  to_import <- !is.na(Rec$anc) & is.na(NewAnc)   # hospital → community
  real_real <- !is.na(Rec$anc) & !is.na(NewAnc)  # hospital → hospital

  # Dummy-safe integer indices for array lookups (NA positions replaced with 1;
  # the resulting values are never read for those cases)
  anc_idx    <- replace(Rec$anc, is.na(Rec$anc), 1L)
  newanc_idx <- replace(NewAnc,  is.na(NewAnc),  1L)
  # Safe beta extraction: S=0 (no valid spatial link) maps to beta=0, not numeric(0)
  old_beta   <- beta_for_s(S[cbind(idx, anc_idx)],    theta$beta)
  new_beta   <- beta_for_s(S[cbind(idx, newanc_idx)], theta$beta)

  # --- Generation time log-likelihood ratio -----------------------------------
  TimeLL <- numeric(n)
  # hospital → hospital: NB ratio (generation times both under NB model)
  TimeLL[real_real] <- ll_time_anc_move_local(
    Rec$inf[real_real] - Rec$inf[anc_idx[real_real]],
    Rec$inf[real_real] - Rec$inf[newanc_idx[real_real]],
    Rec$kap[real_real], old_beta[real_real], new_beta[real_real]
  )
  # community → hospital: NB(new ancestor) − Geo(community prior)
  TimeLL[to_real] <-
    ll_time_nb(Rec$inf[to_real] - Rec$inf[newanc_idx[to_real]],
               Rec$kap[to_real], new_beta[to_real]) -
    ll_time_geo(Rec$inf[to_real], Rec$Adm[to_real], theta$lambda_import)
  # hospital → community: Geo(community prior) − NB(old ancestor)
  TimeLL[to_import] <-
    ll_time_geo(Rec$inf[to_import], Rec$Adm[to_import], theta$lambda_import) -
    ll_time_nb(Rec$inf[to_import] - Rec$inf[anc_idx[to_import]],
               Rec$kap[to_import], old_beta[to_import])

  # --- Genetic log-likelihood ratio -------------------------------------------
  # For NA↔real moves there is no genetic model for community cases (GenLL=0).
  # For real→real moves we split by whether each ancestor is a room (D = Inf):
  #   patient → patient: full two-sided log-ratio (usual case)
  #   room    → patient: one-sided +ll_genetic_single (gain genetic support)
  #   patient → room:    one-sided -ll_genetic_single (lose genetic support)
  #   room    → room:    0 (no genetic information on either side)
  # Zeroing the whole term for room-involved moves (the previous behaviour via
  # ll_genetic_anc_move_local returning NaN→0) discards the one-sided signal.
  GenLL <- numeric(n)
  if (use_genetics) {
    d_old <- D[cbind(idx, anc_idx)]
    d_new <- D[cbind(idx, newanc_idx)]

    pp <- real_real & is.finite(d_old) & is.finite(d_new)   # patient → patient
    rp <- real_real & !is.finite(d_old) & is.finite(d_new)  # room    → patient
    pr <- real_real & is.finite(d_old) & !is.finite(d_new)  # patient → room

    if (any(pp))
      GenLL[pp] <- ll_genetic_anc_move_local(d_old[pp], d_new[pp],
                                              Rec$kap[pp], theta$mu)
    if (any(rp))
      GenLL[rp] <-  ll_genetic_single(d_new[rp], Rec$kap[rp], theta$mu)
    if (any(pr))
      GenLL[pr] <- -ll_genetic_single(d_old[pr], Rec$kap[pr], theta$mu)
  }

  # --- Ancestor selection log-likelihood ratio --------------------------------
  # hospital → hospital: qlogis scores are comparable because Z cancels.
  # NA↔real: qlogis(anc_prob) cannot be compared against pi without the
  # normalisation constant Z = sum_k anc_prob(k); use the prior term only.
  log_or <- log(theta$pi_import / (1 - theta$pi_import))
  AncLL  <- numeric(n)
  for (i in sample.int(n)) {
    if (!real_real[i] && !to_real[i] && !to_import[i]) next
    p_old     <- if (real_real[i]) qlogis(anc_prob(Rec$kap[i], S[i, anc_idx[i]],    C[i, anc_idx[i]],    old_beta[i], TmPowF)) else 0
    p_new     <- if (real_real[i]) qlogis(anc_prob(Rec$kap[i], S[i, newanc_idx[i]], C[i, newanc_idx[i]], new_beta[i], TmPowF)) else 0
    prior_adj <- if (!is_room[i] && to_real[i]) -log_or else if (!is_room[i] && to_import[i]) log_or else 0
    diff      <- p_new - p_old
    AncLL[i]  <- prior_adj + if (is.nan(diff)) 0 else diff
  }

  # --- MH acceptance (vectorised) ---------------------------------------------
  Rec$anc <- ifelse(log(runif(n)) < TimeLL + GenLL + AncLL, NewAnc, Rec$anc)
  Rec
}

propose_kappa <- function(Rec, D, theta, S, C, TmPowF, t_root, use_genetics = TRUE) {
  n      <- nrow(Rec)
  delta  <- 3L
  NewKap <- Rec$kap + sample(-delta:delta, n, replace = TRUE)
  rates  <- node_rates(Rec, theta)

  for (i in sample.int(n)) {
    anc      <- Rec$anc[i]
    anc_tinf <- if (is.na(anc)) t_root else Rec$inf[anc]
    NewKapTmp <- reflect(NewKap[i], 1L, Rec$inf[i] - anc_tinf)

    CaseLL <- ll_detection_kap_move_local(Rec$kap[i], NewKapTmp, theta,
                                          test = rates$test[i],
                                          dis  = rates$dis[i])

    if (is.na(anc)) {
      TimeLL <- 0
      GenLL  <- 0
      AncLL  <- 0
    } else {
      beta   <- beta_for_s(S[i, anc], theta$beta)
      TimeLL <- ll_time_kap_move_local(Rec$kap[i], NewKapTmp,
                                             Rec$inf[i] - Rec$inf[anc], beta)
      GenLL  <- if (use_genetics) ll_genetic_kap_move_local(Rec$kap[i], NewKapTmp, theta$mu, D[i, anc]) else 0
      AncLL  <- ll_ancestry_kap_move_local(Rec$kap[i], NewKapTmp,
                                            S[i, anc], C[i, anc], beta, TmPowF)
    }

    if (log(runif(1)) < CaseLL + TimeLL + GenLL + AncLL)
      Rec$kap[i] <- NewKapTmp
  }
  Rec
}


propose_subtree_shift <- function(Rec, theta, S, t_root) {
  n             <- nrow(Rec)
  idx           <- which(!is.na(Rec$anc))
  children_list <- split(idx, Rec$anc[idx])

  # BFS: collect all descendants of root, inclusive
  get_subtree <- function(root) {
    result <- root
    front  <- root
    repeat {
      ch <- unlist(children_list[as.character(front)], use.names = FALSE)
      if (length(ch) == 0L) break
      result <- c(result, ch)
      front  <- ch
    }
    result
  }

  # n proposals: one randomly chosen subtree root each time
  for (k in seq_len(n)) {
    i       <- sample.int(n, 1L)
    subtree <- get_subtree(i)

    # Tightest window constraint across all subtree members
    max_pos <- min(Rec$PTest[subtree] - Rec$inf[subtree])  # max forward shift
    max_neg <- min(Rec$inf[subtree]   - Rec$Adm[subtree])  # max backward shift

    # Ancestor constraint: inf[i] + delta - anc_inf >= kap[i]
    anc_i   <- Rec$anc[i]
    anc_inf <- if (is.na(anc_i)) t_root else Rec$inf[anc_i]
    max_neg <- min(max_neg, Rec$inf[i] - anc_inf - Rec$kap[i])

    if (max_neg + max_pos == 0L) next  # no movement possible

    # Sample delta uniformly on {-max_neg, ..., max_pos}
    # Proposal ratio = 1 because max_neg + max_pos is conserved under the shift
    delta <- sample.int(max_neg + max_pos + 1L, 1L) - (max_neg + 1L)
    if (delta == 0L) next

    new_inf          <- Rec$inf
    new_inf[subtree] <- Rec$inf[subtree] + delta

    # Log acceptance ratio:
    # test-LL: each subtree member shifts by delta — use per-node test rates
    rates    <- node_rates(Rec[subtree, , drop = FALSE], theta)
    TestLL   <- -delta * sum(log(1 - rates$test))
    # time-LL: only the boundary bond i -> anc[i] changes
    beta_anc <- if (is.na(anc_i)) theta$beta[3L] else beta_for_s(S[i, anc_i], theta$beta)
    TimeLL   <- ll_time_inf_move_local(
      Rec$inf[i] - anc_inf,
      new_inf[i]    - anc_inf,
      Rec$kap[i],
      beta_anc
    )

    # Geometric prior anchoring imported cases toward t_root (patients only)
    ImportLL <- if (is.na(anc_i) && !node_is_room(Rec)[i])
                  delta * log(1 - theta$lambda_import)
                else 0

    if (log(runif(1L)) < TestLL + TimeLL + ImportLL) Rec$inf <- new_inf
  }
  Rec
}

########################
# Likelihood functions #
########################
## Detection likelihood helpers — node-type aware ────────────────────────────

# Returns per-node test and discharge/clean rate vectors.
# For patient-only mode, returns scalar values as before.
# For rooms-and-patients mode, uses node_type column if present.
node_rates <- function(Rec, theta) {
  if ("node_type" %in% names(Rec)) {
    is_room <- !is.na(Rec$node_type) & Rec$node_type == "room"
    # Suppress test-rate likelihood for rooms that have no actual positive swab
    # (has_ptest = FALSE means PTest was set to Dis as a proxy, not a real obs).
    has_pt  <- if ("has_ptest" %in% names(Rec)) Rec$has_ptest else TRUE
    room_test_r <- ifelse(has_pt, theta$P$RoomTest, 0)
    test_r  <- ifelse(is_room, room_test_r, theta$P$Test)
    dis_r   <- ifelse(is_room, theta$P$Clean,       theta$P$Dis)
  } else {
    test_r <- rep(theta$P$Test, nrow(Rec))
    dis_r  <- rep(theta$P$Dis,  nrow(Rec))
  }
  list(test = test_r, dis = dis_r)
}

ll_detection_kap_move_local <- function(ka, kb, theta, test = NULL, dis = NULL) {
  if (is.null(test)) test <- theta$P$Test[1L]
  if (is.null(dis))  dis  <- theta$P$Dis[1L]
  p_dec <- test * (1 - dis) / (test + dis - test * dis)
  (kb - ka) * log(1 - p_dec)
}
ll_testing_inf_move <- function(t, Rec, theta) {
  rates <- node_rates(Rec, theta)
  (Rec[, "inf"] - t) * log(1 - rates$test)
}
ll_time_inf_move_local <- function(ta, tb, kappa, beta) {
  l = sum(lgamma(tb) - lgamma(ta) +
          lgamma(ta - kappa + 1) - lgamma(tb - kappa + 1) +
          (tb - ta) * log(1 - beta))
  if (is.nan(l)) -Inf else l
}
ll_time_anc_move_local <- function(ta, tb, kap, beta_a, beta_b) {
  l <- lgamma(tb) - lgamma(ta) +
       lgamma(ta - kap + 1) - lgamma(tb - kap + 1) +
       kap * log(beta_b / beta_a) +
       (tb - kap) * log(1 - beta_b) - (ta - kap) * log(1 - beta_a)
  ifelse(is.nan(l), -Inf, l)
}
ll_time_kap_move_local <- function(ka, kb, t, beta) {
  # This is exactly ll_time_nb(t, kb, beta) - ll_time_nb(t, ka, beta): the
  # sign()-wrapped form this replaced miscomputed that difference (it flipped
  # the sign of the lgamma(ka)-lgamma(kb) combinatorial-penalty term whenever
  # kb > ka, and zeroed a term entirely for the single most common move,
  # kb == ka + 1) so it rewarded rather than penalised inflating kap.
  l <- ll_time_nb(t, kb, beta) - ll_time_nb(t, ka, beta)
  if (is.nan(l)) 0 else l
}
ll_genetic_anc_move_local <- function(Da, Db, kap, mu) {
  # Inf distance means cross-component (no genetic path) — treat as no information
  l <- (Db - Da) * log(kap * mu) + lgamma(Da + 1) - lgamma(Db + 1)
  ifelse(!is.finite(l), 0, l)
}
ll_genetic_kap_move_local <- function(ka, kb, mu, D) {
  if (ka == kb || !is.finite(D)) return(0)
  D * (log(kb) - log(ka)) - mu * (kb - ka)
}
# Probability that a case with generation kap was infected by a specific ancestor,
# given spatial relationship s and contact duration c.
# Returns 0 for impossible links (s=0 means no valid spatial relationship, beta=0).
anc_prob <- function(kap, s, c, beta, TmPowF) {
  if (s == 0L || beta == 0) return(0)
  if (kap == 1L) 1 - (1 - beta)^c
  else TmPowF[min(kap - 1L, nrow(TmPowF)), min(s, ncol(TmPowF))]
}

# Absolute log-likelihood of NB(kap, beta) generation time t.
ll_time_nb <- function(t, kap, beta) {
  lgamma(t) - lgamma(kap) - lgamma(t - kap + 1) + kap * log(beta) + (t - kap) * log(1 - beta)
}

# Absolute log-likelihood of Geo(lambda) infection time for a community case.
ll_time_geo <- function(inf, adm, lambda) {
  (inf - adm) * log(1 - lambda) + log(lambda)
}

# Absolute Poisson log-likelihood for genomic distance D given kap mutations per step.
ll_genetic_single <- function(D, kap, mu) {
  D * log(kap * mu) - kap * mu - lgamma(D + 1)
}

ll_ancestry_kap_move_local <- function(ka, kb, s, c, beta, TmPowF) {
  l <- qlogis(anc_prob(kb, s, c, beta, TmPowF)) - qlogis(anc_prob(ka, s, c, beta, TmPowF))
  if (is.nan(l)) 0 else l
}
ll_ancestry_anc_move_local <- function(kap, s_old, s_new, c_old, c_new, beta_old, beta_new, TmPowF) {
  l <- qlogis(anc_prob(kap, s_new, c_new, beta_new, TmPowF)) - qlogis(anc_prob(kap, s_old, c_old, beta_old, TmPowF))
  if (is.nan(l)) 0 else l
}
####################
# Post-processing  #
####################

## For each case, compute three summaries of reconstruction quality against ground truth:
##   p_true       — per-case posterior probability assigned to the true ancestor
##   log_score    — mean log(p_true), a proper scoring rule (higher = better; 0 is perfect)
##   mean_p_true  — mean posterior probability of the true ancestor across cases
##   mode_accuracy — fraction of cases where the posterior mode matches truth
ancestry_score <- function(true_anc, anc_trace, adm_times, ptest_times) {
  n <- length(true_anc)
  stopifnot(ncol(anc_trace) == n, length(adm_times) == n, length(ptest_times) == n)

  p_true <- vapply(seq_len(n), function(i) {
    ta <- true_anc[i]
    tr <- anc_trace[, i]
    if (is.na(ta)) mean(is.na(tr)) else mean(!is.na(tr) & tr == ta)
  }, numeric(1))

  # Eligible ancestors from the observed-data perspective:
  # any case j ≠ i admitted before i's positive test (Adm[j] < PTest[i]),
  # plus community (NA). Infection times are latent so are not used here.
  # For room-aware scoring, pass adm_times = CaseRoomRec$Adm (the previous
  # cleaning time — observable lower bound), NOT CaseRoomRec$Infc (the latent
  # contamination time). Using Infc narrows the eligible set relative to what
  # the MCMC actually saw and inflates the lift metric. Use the $adm_times
  # field returned by build_room_aware_truth() to get the right values.
  n_elig   <- vapply(seq_len(n), \(i) sum(adm_times[-i] < ptest_times[i]) + 1L, integer(1))
  p_random <- 1 / n_elig
  lift     <- p_true / p_random   # > 1 means better than random

  mode_anc     <- posterior_mode_anc(anc_trace)
  mode_correct <- ifelse(is.na(true_anc),
                         is.na(mode_anc),
                         !is.na(mode_anc) & mode_anc == true_anc)

  list(
    p_true        = p_true,
    n_elig        = n_elig,
    lift          = lift,                                       # p_true / p_random
    log_skill     = mean(log(pmax(lift, 1e-10))),               # mean log-lift (0 = random)
    log_score     = mean(log(pmax(p_true, 1e-10))),             # raw log score
    mean_p_true   = mean(p_true),
    mode_accuracy = mean(mode_correct)
  )
}

## For each case, return the most frequently sampled ancestor across MCMC iterations.
## NA (community/imported case) is treated as a valid ancestor state and counted.
## Ties are broken by first occurrence in the frequency table (i.e. lowest index).
posterior_mode_anc <- function(anc_trace) {
  apply(anc_trace, 2, function(x) {
    tab  <- table(x, useNA = "always")
    best <- names(tab)[which.max(tab)]
    if (is.na(best)) NA_integer_ else as.integer(best)
  })
}

## For each case, summarise reconstruction quality of posterior infection times
## against ground truth. Returns per-case vectors and aggregate scalars.
##
## true_inf   — integer vector of true infection times (ObsRec$Infc)
## inf_trace  — N_samples × n_cases integer matrix (Y$inf)
## prob       — credible interval levels to report coverage for
##
## Per-case outputs:
##   post_mean / post_mode / post_median — point summaries of the posterior
##   err_mean / err_mode                 — signed error (positive = overestimate)
##   ci_lower / ci_upper                 — equal-tailed CI bounds for each level in prob
##   ci_width                            — CI width per case, per level
##
## Aggregate outputs:
##   mae_mean / mae_mode / mae_median — mean absolute error of each point summary
##   bias_mean                        — mean signed error of posterior mean
##   coverage                         — named vector: fraction of cases where truth
##                                      falls inside each credible interval
##   mean_ci_width                    — mean CI width per level (narrower = sharper)
inf_score <- function(true_inf, inf_trace, prob = c(0.50, 0.90, 0.95)) {
  n <- length(true_inf)
  stopifnot(ncol(inf_trace) == n)

  # ── Per-case point summaries ────────────────────────────────────────────────
  post_mean   <- colMeans(inf_trace)
  post_median <- apply(inf_trace, 2, median)
  post_mode   <- apply(inf_trace, 2, function(x) {
    tab <- table(x)
    as.integer(names(tab)[which.max(tab)])
  })

  err_mean   <- post_mean   - true_inf
  err_mode   <- post_mode   - true_inf
  err_median <- post_median - true_inf

  # ── Per-case credible intervals ─────────────────────────────────────────────
  # ci_lower[[k]] and ci_upper[[k]] are length-n vectors for prob[k]
  lo_probs <- (1 - prob) / 2
  hi_probs <- 1 - lo_probs

  ci_lower <- lapply(lo_probs, function(p) apply(inf_trace, 2, quantile, probs = p))
  ci_upper <- lapply(hi_probs, function(p) apply(inf_trace, 2, quantile, probs = p))
  ci_width <- Map(`-`, ci_upper, ci_lower)
  names(ci_lower) <- names(ci_upper) <- names(ci_width) <- paste0("ci_", prob * 100)

  # ── Aggregate: coverage + mean CI width ─────────────────────────────────────
  coverage <- vapply(seq_along(prob), function(k) {
    mean(true_inf >= ci_lower[[k]] & true_inf <= ci_upper[[k]])
  }, numeric(1))
  names(coverage) <- paste0("cov_", prob * 100)

  mean_ci_width <- vapply(ci_width, mean, numeric(1))

  list(
    # Per-case
    post_mean   = post_mean,
    post_median = post_median,
    post_mode   = post_mode,
    err_mean    = err_mean,
    err_mode    = err_mode,
    err_median  = err_median,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    ci_width    = ci_width,
    # Aggregate
    mae_mean      = mean(abs(err_mean)),
    mae_mode      = mean(abs(err_mode)),
    mae_median    = mean(abs(err_median)),
    bias_mean     = mean(err_mean),
    coverage      = coverage,
    mean_ci_width = mean_ci_width
  )
}

## Summarise reconstruction quality of posterior kappa against ground truth,
## restricted to hospital-acquired cases (non-NA true ancestor) because kappa
## for community imports is measured from the artificial t_root boundary.
##
## true_kap   — integer vector of true generation counts (ObsRec$Gen)
## kap_trace  — N_samples × n_cases integer matrix (Y$kap)
## true_anc   — integer/NA vector of true ancestors (ObsRec$Anc2);
##              if provided, only non-NA cases are scored
## prob       — credible interval levels for coverage reporting
##
## Per-case outputs (indexed to hospital cases only):
##   case_idx                  — original column indices of the scored cases
##   post_mean/median/mode     — point summaries
##   err_mean/median/mode      — signed errors (positive = overestimate)
##   ci_lower / ci_upper       — equal-tailed CI bounds per level
##   ci_width                  — CI width per case per level
##
## Aggregate outputs:
##   n_hospital    — number of cases scored
##   mae_mean / mae_median / mae_mode
##   bias_mean
##   coverage      — named vector of empirical coverage per CI level
##   mean_ci_width — mean CI width per level
kap_score <- function(true_kap, kap_trace, true_anc = NULL,
                      prob = c(0.50, 0.90, 0.95)) {
  n <- length(true_kap)
  stopifnot(ncol(kap_trace) == n)

  idx <- if (!is.null(true_anc)) which(!is.na(true_anc)) else seq_len(n)
  nh  <- length(idx)

  true_k <- true_kap[idx]
  tr     <- kap_trace[, idx, drop = FALSE]

  # ── Point summaries ─────────────────────────────────────────────────────────
  post_mean   <- colMeans(tr)
  post_median <- apply(tr, 2, median)
  post_mode   <- apply(tr, 2, function(x) {
    tab <- table(x)
    as.integer(names(tab)[which.max(tab)])
  })

  err_mean   <- post_mean   - true_k
  err_median <- post_median - true_k
  err_mode   <- post_mode   - true_k

  # ── Credible intervals ───────────────────────────────────────────────────────
  lo <- (1 - prob) / 2
  hi <- 1 - lo
  ci_lower <- lapply(lo, function(p) apply(tr, 2, quantile, probs = p))
  ci_upper <- lapply(hi, function(p) apply(tr, 2, quantile, probs = p))
  ci_width <- Map(`-`, ci_upper, ci_lower)
  nm <- paste0("ci_", prob * 100)
  names(ci_lower) <- names(ci_upper) <- names(ci_width) <- nm

  # ── Coverage ─────────────────────────────────────────────────────────────────
  coverage <- vapply(seq_along(prob), function(k)
    mean(true_k >= ci_lower[[k]] & true_k <= ci_upper[[k]]), numeric(1))
  names(coverage) <- paste0("cov_", prob * 100)

  mean_ci_width <- vapply(ci_width, mean, numeric(1))

  list(
    n_hospital    = nh,
    case_idx      = idx,
    post_mean     = post_mean,
    post_median   = post_median,
    post_mode     = post_mode,
    err_mean      = err_mean,
    err_median    = err_median,
    err_mode      = err_mode,
    ci_lower      = ci_lower,
    ci_upper      = ci_upper,
    ci_width      = ci_width,
    mae_mean      = mean(abs(err_mean)),
    mae_median    = mean(abs(err_median)),
    mae_mode      = mean(abs(err_mode)),
    bias_mean     = mean(err_mean),
    coverage      = coverage,
    mean_ci_width = mean_ci_width
  )
}

####################
# Helper functions #
####################
get_children <- function(anc) {
  idx <- which(!is.na(anc))
  split(idx, anc[idx])
}

## Returns a logical vector: TRUE for room nodes, FALSE for patient nodes.
# Works for both patient-only records (no node_type column) and combined records.
node_is_room <- function(Rec) {
  if ("node_type" %in% names(Rec))
    !is.na(Rec$node_type) & Rec$node_type == "room"
  else
    rep(FALSE, nrow(Rec))
}

## Safe extraction of beta values for a vector of spatial codes.
# S = 0 means no valid spatial link (different room, no contamination path).
# These are mapped to beta = 0 rather than causing a 0-index error.
#
# S codes and their correct beta indices:
#   1 → BrB [1]  patient ← patient, same room
#   2 → BwB [2]  patient ← patient, same ward
#   3 → BhB [3]  patient ← patient, hospital-wide
#   4 → RB  [5]  patient ← room    (room-to-bed rate, NOT bed-to-room)
#   5 → BR  [4]  room    ← patient (bed-to-room rate, NOT room-to-bed)
#   6 → RR  [6]  room    ← room
# Codes 4 and 5 are intentionally crossed: the direct index s → beta[s]
# swaps BR and RB, which is wrong once room betas are non-zero.
beta_for_s <- function(s_vec, beta_vec) {
  beta_idx         <- c(1L, 2L, 3L, 5L, 4L, 6L)  # S-code → beta position
  result           <- beta_vec[beta_idx[pmax(s_vec, 1L)]]
  result[s_vec == 0L] <- 0
  result
}
reflect <- function(x, a, b) {
  out <- x
  same <- a == b
  out[same] <- a[same]
  diff <- !same
  L <- b[diff] - a[diff]
  y <- (x[diff] - a[diff]) %% (2 * L)
  out[diff] <- a[diff] + ifelse(y <= L, y, 2 * L - y)
  out
}