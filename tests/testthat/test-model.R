library(testthat)
source(here::here("R", "simulation.R"))
source(here::here("R", "inference.R"))

# Helper: build a minimal Cases data frame for uplift_ancestry tests.
# Row names are set to `ids` to match what sim_postprocess produces.
make_cases <- function(ids, ancs, infc_times, ptests, muts) {
  data.frame(
    Id    = ids,
    Bed   = ids,
    Room  = rep("R1", length(ids)),
    Ward  = rep("W1", length(ids)),
    Adm   = 0L,
    Dis   = 10L,
    Infc  = infc_times,
    NTest = NA_integer_,
    PTest = ptests,
    Anc   = ancs,
    Mut   = muts,
    row.names = ids
  )
}

# ── posterior_mode_anc ────────────────────────────────────────────────────────

test_that("posterior_mode_anc: returns modal ancestor for each case", {
  # Case 1 mostly has ancestor 2; case 2 mostly has ancestor 1
  anc <- matrix(c(2L, 2L, 2L, 1L,
                  1L, 1L, 2L, 1L), nrow = 4, ncol = 2)
  result <- posterior_mode_anc(anc)
  expect_equal(result, c(2L, 1L))
})

test_that("posterior_mode_anc: NA (community case) counts as a valid state", {
  anc <- matrix(c(NA, NA, NA, 1L), nrow = 4, ncol = 1)
  result <- posterior_mode_anc(anc)
  expect_true(is.na(result[1]))
})

test_that("posterior_mode_anc: NA wins when it is the plurality", {
  anc <- matrix(c(NA, NA, NA, 1L, 2L), nrow = 5, ncol = 1)
  result <- posterior_mode_anc(anc)
  expect_true(is.na(result[1]))
})

test_that("posterior_mode_anc: returns integer vector of length n_cases", {
  anc <- matrix(c(1L, 2L, 1L, 2L, 3L, 3L), nrow = 3, ncol = 2)
  result <- posterior_mode_anc(anc)
  expect_type(result, "integer")
  expect_length(result, 2)
})

# ── prob_of_inf ────────────────────────────────────────────────────────────────

test_that("prob_of_inf: no infected beds gives zero probabilities", {
  # rowSums of a zero-column matrix = 0, so 1 - exp(0) = 0 for all susceptibles
  LogComp <- matrix(log(0.5), nrow = 3, ncol = 3)
  Idx <- list(S = 1:3, I = integer(0))
  expect_equal(prob_of_inf(LogComp, Idx), rep(0, 3))
})

test_that("prob_of_inf: no susceptible beds gives empty result", {
  LogComp <- matrix(log(0.5), nrow = 3, ncol = 3)
  Idx <- list(S = integer(0), I = 1:3)
  expect_equal(prob_of_inf(LogComp, Idx), numeric(0))
})

test_that("prob_of_inf: single susceptible, single infected with known weight", {
  # w = 0.5  →  LogComp = log(1 - 0.5)  →  prob = 1 - exp(log(0.5)) = 0.5
  LogComp <- matrix(log(0.5), nrow = 2, ncol = 2)
  Idx <- list(S = 1, I = 2)
  expect_equal(prob_of_inf(LogComp, Idx), 0.5)
})

test_that("prob_of_inf: two infecteds combine as independent competing hazards", {
  # w1=0.3, w2=0.4  →  prob = 1 - (1-0.3)*(1-0.4) = 0.58
  LogComp <- matrix(0, nrow = 3, ncol = 3)
  LogComp[1, 2] <- log(1 - 0.3)
  LogComp[1, 3] <- log(1 - 0.4)
  Idx <- list(S = 1, I = c(2, 3))
  expect_equal(prob_of_inf(LogComp, Idx), 1 - 0.7 * 0.6, tolerance = 1e-10)
})

test_that("prob_of_inf: w=1 gives certain infection (prob = 1)", {
  LogComp <- matrix(-Inf, nrow = 2, ncol = 2)
  Idx <- list(S = 1, I = 2)
  expect_equal(prob_of_inf(LogComp, Idx), 1)
})

# ── sample_ancestors ───────────────────────────────────────────────────────────

test_that("sample_ancestors: empty NewInfIdx returns character(0)", {
  Odds <- matrix(1, nrow = 3, ncol = 3)
  Idx  <- list(S = c(1, 3), I = 2)
  result <- sample_ancestors(Idx, integer(0), c("a", "b", "c"), Odds)
  expect_equal(result, character(0))
})

test_that("sample_ancestors: single infected is always chosen as ancestor", {
  # Only one infected bed; it must be sampled regardless of odds value
  Odds <- matrix(1, nrow = 3, ncol = 3)
  Idx  <- list(S = c(1, 3), I = 2)
  Id   <- c("a", "b", "c")
  expect_equal(sample_ancestors(Idx, 1L, Id, Odds), "b")
})

test_that("sample_ancestors: returns character vector with one entry per new infection", {
  Odds <- matrix(1, nrow = 4, ncol = 4)
  Idx  <- list(S = c(1, 2), I = c(3, 4))
  Id   <- c("a", "b", "c", "d")
  result <- sample_ancestors(Idx, c(1L, 2L), Id, Odds)
  expect_type(result, "character")
  expect_length(result, 2)
  expect_true(all(result %in% c("c", "d")))
})

test_that("sample_ancestors: extreme odds consistently select the dominant infector", {
  set.seed(8823)
  Odds <- matrix(0, nrow = 3, ncol = 3)
  Odds[1, 2] <- 1e9   # bed 2 overwhelmingly likely
  Odds[1, 3] <- 1e-9
  Idx <- list(S = 1, I = c(2, 3))
  Id  <- c("a", "b", "c")
  results <- replicate(200, sample_ancestors(Idx, 1L, Id, Odds))
  expect_true(mean(results == "b") > 0.99)
})

# ── create_record ──────────────────────────────────────────────────────────────

test_that("create_record: output has T * N rows and all expected columns", {
  N <- 3; T <- 5
  Positions <- data.frame(
    Room = paste0("R", 1:N), Ward = rep("W1", N),
    row.names = as.character(1:N)
  )
  State <- matrix(0, nrow = N, ncol = T + 1)
  Rec <- create_record(Positions, State, as.character(1:N), T)
  expect_equal(nrow(Rec), T * N)
  expect_named(Rec, c("Id","Bed","Room","Ward","Adm","Dis","Infc","NTest","PTest","Anc","Mut"))
})

test_that("create_record: initial patients have Id, Bed, Room, Ward, Adm set", {
  N <- 2; T <- 3
  Positions <- data.frame(
    Room = c("R1", "R2"), Ward = c("W1", "W1"),
    row.names = c("1", "2")
  )
  State <- matrix(0, nrow = N, ncol = T + 1)
  Rec <- create_record(Positions, State, c("1", "2"), T)
  expect_equal(Rec["1", "Id"],   "1")
  expect_equal(Rec["1", "Bed"],  "1")
  expect_equal(Rec["1", "Room"], "R1")
  expect_equal(Rec["1", "Ward"], "W1")
  expect_equal(Rec["1", "Adm"],  0L)
  expect_equal(Rec["2", "Room"], "R2")
})

test_that("create_record: initially infected beds get Infc=0 and Anc='0'", {
  N <- 2; T <- 3
  Positions <- data.frame(
    Room = c("R1", "R2"), Ward = c("W1", "W1"),
    row.names = c("1", "2")
  )
  State <- matrix(0, nrow = N, ncol = T + 1)
  State[1, 1] <- 1  # bed "1" infected at t=0
  Rec <- create_record(Positions, State, c("1", "2"), T)
  expect_equal(Rec["1", "Infc"], 0L)
  expect_equal(Rec["1", "Anc"],  "0")
  expect_true(is.na(Rec["2", "Infc"]))
  expect_equal(Rec["2", "Anc"],  "")  # susceptible, no ancestor
})

test_that("create_record: future patient slots are empty/NA", {
  N <- 2; T <- 3
  Positions <- data.frame(
    Room = c("R1", "R2"), Ward = c("W1", "W1"),
    row.names = c("1", "2")
  )
  State <- matrix(0, nrow = N, ncol = T + 1)
  Rec <- create_record(Positions, State, c("1", "2"), T)
  # Rows "3"–"6" are future admission slots, should be uninitialised
  expect_equal(Rec["3", "Id"],  "")
  expect_true(is.na(Rec["3", "Adm"]))
  expect_true(is.na(Rec["3", "Infc"]))
})

# ── uplift_ancestry ────────────────────────────────────────────────────────────

test_that("uplift_ancestry: root case gets Anc2='0', Mut2=own Mut, Gen=1", {
  Cases  <- make_cases("1", "0", 0L, 1L, 3L)
  result <- uplift_ancestry(Cases)
  expect_equal(result["1", "Anc2"], "0")
  expect_equal(result["1", "Mut2"], 3L)
  expect_equal(result["1", "Gen"],  1L)
})

test_that("uplift_ancestry: tested direct parent → Anc2=parent, Mut2=own Mut, Gen=1", {
  Cases  <- make_cases(c("1","2"), c("0","1"), c(0L,1L), c(1L,2L), c(2L,3L))
  result <- uplift_ancestry(Cases)
  expect_equal(result["2", "Anc2"], "1")
  expect_equal(result["2", "Mut2"], 3L)
  expect_equal(result["2", "Gen"],  1L)
})

test_that("uplift_ancestry: single unobserved parent is skipped, mutations accumulated", {
  # Tree: 1(root,tested) → 2(untested) → 3(tested)
  Cases <- make_cases(
    c("1","2","3"), c("0","1","2"),
    c(0L, 1L, 2L), c(1L, NA_integer_, 3L), c(1L, 2L, 3L)
  )
  result <- uplift_ancestry(Cases)
  # Case 2: parent 1 is tested → direct link, no chain
  expect_equal(result["2", "Anc2"], "1")
  expect_equal(result["2", "Mut2"], 2L)
  expect_equal(result["2", "Gen"],  1L)
  # Case 3: parent 2 is untested → skip to 1, accumulate mut[3] + mut[2]
  expect_equal(result["3", "Anc2"], "1")
  expect_equal(result["3", "Mut2"], 3L + 2L)
  expect_equal(result["3", "Gen"],  2L)
})

test_that("uplift_ancestry: long unobserved chain accumulates mutations and increments Gen", {
  # Tree: 1(tested) → 2(untested) → 3(untested) → 4(tested)
  Cases <- make_cases(
    c("1","2","3","4"), c("0","1","2","3"),
    c(0L, 1L, 2L, 3L), c(1L, NA, NA, 4L), c(1L, 2L, 3L, 4L)
  )
  result <- uplift_ancestry(Cases)
  expect_equal(result["4", "Anc2"], "1")
  expect_equal(result["4", "Mut2"], 4L + 3L + 2L)  # own + all intermediate
  expect_equal(result["4", "Gen"],  3L)
})

test_that("uplift_ancestry: output row count equals input row count", {
  Cases <- make_cases(c("1","2","3"), c("0","1","2"), c(0L,1L,2L),
                      c(1L, NA, 3L), c(1L, 2L, 3L))
  expect_equal(nrow(uplift_ancestry(Cases)), 3L)
})

test_that("uplift_ancestry: output is sorted by Infc ascending", {
  # Input given in reverse infection-time order
  Cases <- make_cases(c("3","2","1"), c("2","1","0"), c(2L,1L,0L),
                      c(3L, NA, 1L), c(3L, 2L, 1L))
  result <- uplift_ancestry(Cases)
  expect_equal(result$Infc, c(0L, 1L, 2L))
})

test_that("uplift_ancestry: Anc2, Mut2, and Gen columns are added to output", {
  Cases  <- make_cases("1", "0", 0L, 1L, 1L)
  result <- uplift_ancestry(Cases)
  expect_true(all(c("Anc2", "Mut2", "Gen") %in% names(result)))
})

# ── create_room_record ────────────────────────────────────────────────────────

make_room_positions <- function(num_beds, num_rooms) {
  all_names <- c(as.character(seq_len(num_beds)),
                 paste0("r", seq_len(num_rooms)))
  data.frame(
    Room = c(paste0("R", rep(seq_len(num_rooms), each = num_beds %/% num_rooms)),
             paste0("R", seq_len(num_rooms))),
    Ward = rep("W1", num_beds + num_rooms),
    row.names = all_names,
    stringsAsFactors = FALSE
  )
}

test_that("create_room_record: output has correct columns including Adm", {
  Positions <- make_room_positions(4, 2)
  State     <- matrix(0, nrow = 6, ncol = 6)
  RoomId    <- c("r1_ep1", "r2_ep1")
  Rec       <- create_room_record(Positions, State, RoomId, 4, 2, 5)
  expect_named(Rec, c("Id","Room","Ward","Adm","Infc","Clean","NTest","PTest","Anc"))
})

test_that("create_room_record: initial episodes have Adm = 0", {
  Positions <- make_room_positions(4, 2)
  State     <- matrix(0, nrow = 6, ncol = 6)
  RoomId    <- c("r1_ep1", "r2_ep1")
  Rec       <- create_room_record(Positions, State, RoomId, 4, 2, 5)
  expect_equal(Rec["r1_ep1", "Adm"], 0L)
  expect_equal(Rec["r2_ep1", "Adm"], 0L)
})

test_that("create_room_record: initial episode rows have Id/Room/Ward set", {
  Positions <- make_room_positions(4, 2)
  State     <- matrix(0, nrow = 6, ncol = 6)
  RoomId    <- c("r1_ep1", "r2_ep1")
  Rec       <- create_room_record(Positions, State, RoomId, 4, 2, 5)
  expect_equal(Rec["r1_ep1", "Id"],   "r1_ep1")
  expect_equal(Rec["r1_ep1", "Room"], "R1")
  expect_equal(Rec["r2_ep1", "Room"], "R2")
})

test_that("create_room_record: initially contaminated rooms get Infc=0 and Anc='0'", {
  Positions <- make_room_positions(4, 2)
  State     <- matrix(0, nrow = 6, ncol = 6)
  State[5, 1] <- 1  # room index 1 (global row 5 = NumBeds + 1)
  RoomId    <- c("r1_ep1", "r2_ep1")
  Rec       <- create_room_record(Positions, State, RoomId, 4, 2, 5)
  expect_equal(Rec["r1_ep1", "Infc"], 0L)
  expect_equal(Rec["r1_ep1", "Anc"],  "0")
  expect_true(is.na(Rec["r2_ep1", "Infc"]))
})

# ── uplift_ancestry with room intermediaries ──────────────────────────────────

make_room_cases <- function(ids, ancs, infc_times, ptests) {
  data.frame(
    Id    = ids,
    Room  = rep("R1", length(ids)),
    Ward  = rep("W1", length(ids)),
    Infc  = infc_times,
    Clean = NA_integer_,
    NTest = NA_integer_,
    PTest = ptests,
    Anc   = ancs,
    stringsAsFactors = FALSE,
    row.names = ids
  )
}

test_that("uplift_ancestry: unobserved room intermediary is skipped", {
  # Chain: patient 1 (tested) -> room r1_ep1 (not tested) -> patient 2 (tested)
  cases <- make_cases(c("1","2"), c("0","r1_ep1"), c(0L, 2L),
                      c(1L, 3L), c(1L, 2L))
  room  <- make_room_cases("r1_ep1", "1", 1L, NA_integer_)
  result <- uplift_ancestry(cases, room_cases = room)
  # Patient 2's Anc2 should skip the room and point to patient 1
  expect_equal(result[result$Id == "2", "Anc2"], "1")
  expect_equal(result[result$Id == "2", "Gen"],  2L)
})

test_that("uplift_ancestry: row count unchanged when room_cases supplied", {
  cases  <- make_cases(c("1","2"), c("0","r1_ep1"), c(0L, 2L),
                       c(1L, 3L), c(1L, 2L))
  room   <- make_room_cases("r1_ep1", "1", 1L, NA_integer_)
  result <- uplift_ancestry(cases, room_cases = room)
  expect_equal(nrow(result), 2L)
})

test_that("uplift_ancestry: NULL room_cases reproduces original behaviour", {
  cases  <- make_cases(c("1","2","3"), c("0","1","2"),
                       c(0L, 1L, 2L), c(1L, NA, 3L), c(1L, 2L, 3L))
  expect_equal(uplift_ancestry(cases, room_cases = NULL),
               uplift_ancestry(cases))
})

# ── node_is_room ──────────────────────────────────────────────────────────────

test_that("node_is_room: FALSE for all rows when node_type column absent", {
  Rec <- data.frame(Id = c("1","2"), Room = c("R1","R1"), Ward = c("W1","W1"))
  expect_equal(node_is_room(Rec), c(FALSE, FALSE))
})

test_that("node_is_room: correctly identifies room rows", {
  Rec <- data.frame(
    Id        = c("1","2","3"),
    node_type = c(NA_character_, "room", NA_character_),
    stringsAsFactors = FALSE
  )
  expect_equal(node_is_room(Rec), c(FALSE, TRUE, FALSE))
})

# ── beta_for_s ────────────────────────────────────────────────────────────────
# S-codes 4 and 5 use crossed beta indices:
#   4 = patient <- room  → should return beta[5] = RB (room-to-bed rate)
#   5 = room <- patient  → should return beta[4] = BR (bed-to-room rate)
# The direct index s → beta[s] swap was the original bug.

test_that("beta_for_s: S=0 returns 0 (no spatial link)", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)
  expect_equal(beta_for_s(0L, beta), 0)
})

test_that("beta_for_s: S=1,2,3 map to BrB, BwB, BhB directly", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)
  expect_equal(beta_for_s(1L, beta), 0.08)
  expect_equal(beta_for_s(2L, beta), 0.005)
  expect_equal(beta_for_s(3L, beta), 0.0005)
})

test_that("beta_for_s: S=4 (patient<-room) returns RB = beta[5], not BR = beta[4]", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)  # BR=0.10, RB=0.15
  expect_equal(beta_for_s(4L, beta), 0.15)  # RB
  expect_false(isTRUE(all.equal(beta_for_s(4L, beta), 0.10)))  # not BR
})

test_that("beta_for_s: S=5 (room<-patient) returns BR = beta[4], not RB = beta[5]", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)
  expect_equal(beta_for_s(5L, beta), 0.10)  # BR
  expect_false(isTRUE(all.equal(beta_for_s(5L, beta), 0.15)))  # not RB
})

test_that("beta_for_s: S=6 (room<-room) returns RR = beta[6]", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)
  expect_equal(beta_for_s(6L, beta), 0.02)
})

test_that("beta_for_s: vectorised over mixed S codes", {
  beta <- c(0.08, 0.005, 0.0005, 0.10, 0.15, 0.02)
  result <- beta_for_s(c(0L, 1L, 4L, 5L), beta)
  expect_equal(result, c(0, 0.08, 0.15, 0.10))
})

# ── harmonise_room_rows ───────────────────────────────────────────────────────

make_room_source <- function(adm, infc, clean, ptest = NA_integer_) {
  data.frame(
    Room  = "R1",
    Ward  = "W1",
    Adm   = adm,
    Infc  = infc,
    Clean = clean,
    NTest = NA_integer_,
    PTest = ptest,
    stringsAsFactors = FALSE
  )
}

test_that("harmonise_room_rows: Adm comes from room_rec$Adm, not $Infc", {
  rs  <- make_room_source(adm = 3L, infc = 7L, clean = 20L)
  out <- harmonise_room_rows(rs, n_pat_rows = 5L)
  expect_equal(out$Adm,  3L)   # previous cleaning time
  expect_false(out$Adm == 7L)  # must NOT be Infc
})

test_that("harmonise_room_rows: Infc is NA (latent, not passed to MCMC)", {
  rs  <- make_room_source(adm = 3L, infc = 7L, clean = 20L)
  out <- harmonise_room_rows(rs, n_pat_rows = 5L)
  expect_true(is.na(out$Infc))
})

test_that("harmonise_room_rows: Dis equals Clean", {
  rs  <- make_room_source(adm = 0L, infc = 5L, clean = 15L)
  out <- harmonise_room_rows(rs, n_pat_rows = 0L)
  expect_equal(out$Dis, 15L)
})

test_that("harmonise_room_rows: PTest used when swab present, else Clean", {
  rs_swab    <- make_room_source(adm = 0L, infc = 2L, clean = 10L, ptest = 6L)
  rs_no_swab <- make_room_source(adm = 0L, infc = 2L, clean = 10L)
  expect_equal(harmonise_room_rows(rs_swab,    0L)$PTest, 6L)
  expect_equal(harmonise_room_rows(rs_no_swab, 0L)$PTest, 10L)
  expect_true( harmonise_room_rows(rs_swab,    0L)$has_ptest)
  expect_false(harmonise_room_rows(rs_no_swab, 0L)$has_ptest)
})

# ── build_room_aware_truth ────────────────────────────────────────────────────

test_that("build_room_aware_truth: returns a named list with three fields", {
  source(here::here("R", "analysis.R"))
  source(here::here("R", "hospital.R"))
  source(here::here("config", "hospital_two_wards.R"))
  source(here::here("config", "parameters.R"))

  theta_room <- local({
    th <- theta; th$beta["BR"] <- 0.10; th$beta["RB"] <- 0.15
    th$use_room_transmission <- TRUE
    W <- weights(th$beta, Contact)
    th$LogComp <- log1p(-W); th$Odds <- W / (1 - W); th
  })
  th_sim <- local({
    th <- theta_room
    th$P$Init <- c(rep(1, 3), rep(0, NumBeds - 3), rep(0, NumRooms)); th
  })
  set.seed(37675)
  ob <- simulate_outbreak(45, th_sim)

  truth <- build_room_aware_truth(ob, use_room_tests = FALSE)
  expect_type(truth, "list")
  expect_named(truth, c("true_anc", "adm_times", "ptest_times"))
})

test_that("build_room_aware_truth: adm_times for rooms equals CaseRoomRec$Adm", {
  source(here::here("R", "analysis.R"))
  source(here::here("R", "hospital.R"))
  source(here::here("config", "hospital_two_wards.R"))
  source(here::here("config", "parameters.R"))

  theta_room <- local({
    th <- theta; th$beta["BR"] <- 0.10; th$beta["RB"] <- 0.15
    th$use_room_transmission <- TRUE
    W <- weights(th$beta, Contact)
    th$LogComp <- log1p(-W); th$Odds <- W / (1 - W); th
  })
  th_sim <- local({
    th <- theta_room
    th$P$Init <- c(rep(1, 3), rep(0, NumBeds - 3), rep(0, NumRooms)); th
  })
  set.seed(37675)
  ob <- simulate_outbreak(45, th_sim)

  truth    <- build_room_aware_truth(ob, use_room_tests = FALSE)
  n_pat    <- nrow(ob$ObsRec)
  room_adm <- truth$adm_times[(n_pat + 1):length(truth$adm_times)]

  expect_equal(room_adm, ob$CaseRoomRec$Adm)
  expect_true(all(room_adm <= ob$CaseRoomRec$Infc, na.rm = TRUE))
})
