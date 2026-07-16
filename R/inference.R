mcmc <- function(outbreak_data, N_samples, burn_in = 0, use_genetics = TRUE) {
  mcmc_item <- mcmc_setup(outbreak_data, N_samples, burn_in, use_genetics)
  mcmc_item <- mcmc_loop(mcmc_item)
  mcmc_cleanup(mcmc_item)
}
mcmc_setup <- function(outbreak_data, N_samples, burn_in, use_genetics = TRUE) {
  ## Extraction from outbreak_data
  Rec   <- outbreak_data$ObsRec
  D     <- outbreak_data$ObsDist
  theta <- outbreak_data$theta

  # P$Test and P$Dis are per-bed vectors in the simulation; MCMC needs scalars
  theta$P$Test <- theta$P$Test[[1L]]
  theta$P$Dis  <- theta$P$Dis[[1L]]

  ## Fix rec
  Rec = init(Rec)
  S = create_s_matrix(Rec)
  C = create_c_matrix(Rec)
  t_root = min(Rec$Adm) - 1L

  ## Create transition matrices
  Tms <- create_transition_matrices(theta, kap_max = max(Rec$PTest) - min(Rec$Adm) + 1L)

  ## Create trace object
  n_cases <- nrow(Rec)
  Trace <- list(inf = matrix(NA_integer_, N_samples, n_cases),
                anc = matrix(NA_integer_, N_samples, n_cases),
                kap = matrix(NA_integer_, N_samples, n_cases))
  list(Rec          = Rec,
       theta        = theta,
       Dist         = D,
       S            = S,
       C            = C,
       t_root       = t_root,
       Tms          = Tms,
       Trace        = Trace,
       N_samples    = N_samples,
       burn_in      = burn_in,
       use_genetics = use_genetics)
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

  # Hoist: these only depend on fixed quantities (theta, S, Rec$anc, Rec$kap)
  # beta_anc only used for hospital-acquired cases (non-NA anc); NA-safe indexing via replace
  anc_idx      <- replace(Rec$anc, is.na(Rec$anc), 1L)
  beta_anc     <- theta$beta[S[cbind(idx, anc_idx)]]   # imported slots are never read
  kap_ch_list  <- lapply(idx, \(i) Rec$kap[children_list[[as.character(i)]]])
  beta_ch_list <- lapply(idx, \(i) theta$beta[S[i, children_list[[as.character(i)]]]])

  p      <- generate_inf_proposals(Rec, theta)
  NewInf <- p$NewInf
  Accept <- p$Accept

  for (i in sample.int(n)) {
    if (is.na(Rec$anc[i])) {
      own_ta <- own_tb <- own_kap <- own_beta <- NULL
      ImportLL <- (NewInf[i] - Rec$inf[i]) * log(1 - theta$lambda_import)
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
  n   <- nrow(Rec)
  idx <- seq_len(n)

  # --- Propose new ancestors --------------------------------------------------
  # Eligible ancestors: cases infected at least kap steps earlier.
  # NA (community) is included as one of n_elig+1 candidates (index 0) so that
  # hospital ↔ community moves are possible. Proposal is symmetric → ratio cancels.
  order_inf <- order(Rec$inf)
  n_elig    <- findInterval(Rec$inf - Rec$kap, Rec$inf[order_inf])
  sampled   <- as.integer(runif(n, 0, n_elig + 1L))
  NewAnc    <- ifelse(sampled == 0L, NA_integer_, order_inf[pmax(sampled, 1L)])

  # Classify transition type for each case
  to_real   <- is.na(Rec$anc) & !is.na(NewAnc)   # community → hospital
  to_import <- !is.na(Rec$anc) & is.na(NewAnc)   # hospital → community
  real_real <- !is.na(Rec$anc) & !is.na(NewAnc)  # hospital → hospital

  # Dummy-safe integer indices for array lookups (NA positions replaced with 1;
  # the resulting values are never read for those cases)
  anc_idx    <- replace(Rec$anc, is.na(Rec$anc), 1L)
  newanc_idx <- replace(NewAnc,  is.na(NewAnc),  1L)
  old_beta   <- theta$beta[S[cbind(idx, anc_idx)]]
  new_beta   <- theta$beta[S[cbind(idx, newanc_idx)]]

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
  # Only hospital → hospital moves: the kap*mu term cancels in the ratio,
  # keeping this well-behaved. For NA↔real moves there is no genetic model
  # for community cases, so the term is 0 (stays initialised at zero).
  GenLL <- numeric(n)
  if (use_genetics)
    GenLL[real_real] <- ll_genetic_anc_move_local(
      D[cbind(idx[real_real], anc_idx[real_real])],
      D[cbind(idx[real_real], newanc_idx[real_real])],
      Rec$kap[real_real], theta$mu
    )

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
    prior_adj <- if (to_real[i]) -log_or else if (to_import[i]) log_or else 0
    AncLL[i]  <- prior_adj + p_new - p_old
  }

  # --- MH acceptance (vectorised) ---------------------------------------------
  Rec$anc <- ifelse(log(runif(n)) < TimeLL + GenLL + AncLL, NewAnc, Rec$anc)
  Rec
}

propose_kappa <- function(Rec, D, theta, S, C, TmPowF, t_root, use_genetics = TRUE) {
  n      <- nrow(Rec)
  delta  <- 3L
  NewKap <- Rec$kap + sample(-delta:delta, n, replace = TRUE)

  for (i in sample.int(n)) {
    anc      <- Rec$anc[i]
    anc_tinf <- if (is.na(anc)) t_root else Rec$inf[anc]
    NewKapTmp <- reflect(NewKap[i], 1L, Rec$inf[i] - anc_tinf)

    CaseLL <- ll_detection_kap_move_local(Rec$kap[i], NewKapTmp, theta)

    if (is.na(anc)) {
      TimeLL <- 0
      GenLL  <- 0
      AncLL  <- 0
    } else {
      beta   <- theta$beta[S[i, anc]]
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
    # test-LL: each subtree member shifts by delta (proposal always feasible)
    TestLL   <- -delta * length(subtree) * log(1 - theta$P$Test)
    # time-LL: only the boundary bond i -> anc[i] changes
    beta_anc <- if (is.na(anc_i)) theta$beta[3L] else theta$beta[S[i, anc_i]]
    TimeLL   <- ll_time_inf_move_local(
      Rec$inf[i] - anc_inf,
      new_inf[i]    - anc_inf,
      Rec$kap[i],
      beta_anc
    )

    # Geometric prior anchoring imported cases toward t_root
    ImportLL <- if (is.na(anc_i)) delta * log(1 - theta$lambda_import) else 0

    if (log(runif(1L)) < TestLL + TimeLL + ImportLL) Rec$inf <- new_inf
  }
  Rec
}

########################
# Likelihood functions #
########################
ll_detection_kap_move_local <- function(ka, kb, theta) {
  # Use scalar rates; if beds ever have heterogeneous rates this needs to be bed-specific
  test  <- theta$P$Test[1L]
  dis   <- theta$P$Dis[1L]
  p_dec <- test * (1 - dis) / (test + dis - test * dis)
  (kb - ka) * log(1 - p_dec)
}
ll_testing_inf_move <- function(t, Rec, theta) {
  (Rec[, "inf"] - t) * log(1 - theta$P$Test)
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
  sign(ka - 1 - kb) * (lgamma(ka) - lgamma(kb)) +
    sign(t - ka - (t - kb + 1)) * (lgamma(t - ka + 1) - lgamma(t - kb + 1)) +
    (kb - ka) * (log(beta) - log(1 - beta))
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
anc_prob <- function(kap, s, c, beta, TmPowF) {
  if (kap == 1L) 1 - (1 - beta)^c
  else TmPowF[min(kap - 1L, nrow(TmPowF)), s]
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