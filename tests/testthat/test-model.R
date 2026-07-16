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
