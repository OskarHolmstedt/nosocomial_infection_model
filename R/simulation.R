## Function for simulating a nosocomial infectious disease outbreak
# T         - number of days to simulate
# theta     - model parameters (beta, P, LogComp, Odds, mu, Positions)
simulate_outbreak = function(T, theta) {
  sim = sim_setup(T, theta)
  sim = sim_loop(sim, T, theta)
  sim_postprocess(sim, T, theta)
}
sim_setup <- function(T, theta) {
  N     = nrow(theta$Positions)
  State = Matrix(0,
                 nrow = N,
                 ncol = T + 1,
                 dimnames = list(rownames(theta$Positions), NULL))
  State[, 1] = rbern(N, theta$P$Init)
  Id  = rownames(theta$Positions)
  Rec = create_record(theta$Positions, State, Id, T)
  list(State = State, Id = Id, Rec = Rec)
}

sim_loop <- function(sim, T, theta) {
  State   = sim$State
  Id      = sim$Id
  Rec     = sim$Rec
  NumBeds = length(theta$P$Dis)

  for (t in 1:T) {
    col_t = as.vector(State[, t])
    Idx   = get_idx(col_t)

    # Build next state: infected carry forward, susceptibles may become infected
    P_Inf           = prob_of_inf(theta$LogComp, Idx)
    next_col        = col_t
    next_col[Idx$S] = rbern(length(Idx$S), P_Inf)

    # Record new infections and sample ancestors
    NewInfIdx           = which(next_col > col_t)
    NewIds              = Id[NewInfIdx]
    Rec[NewIds, "Infc"] = t
    Rec[NewIds, "Anc"]  = sample_ancestors(Idx, NewInfIdx, Id, theta$Odds)

    # Test
    TestIdx = which(rbern(NumBeds, theta$P$Test) == 1)
    TestRes = col_t[TestIdx]
    Rec[Id[TestIdx[TestRes == 0]], "NTest"] = t
    PosIds  = Id[TestIdx[TestRes == 1]]
    Rec[PosIds[is.na(Rec[PosIds, "PTest"])], "PTest"] = t

    # Discharge: overwrite admitted beds in next_col, then assign whole column
    DisIdx           = which(rbern(NumBeds, theta$P$Dis) == 1)
    NextIds          = as.numeric(Id[DisIdx]) + NumBeds
    Rec[Id[DisIdx], "Dis"]  = t
    next_col[DisIdx] = rbern(length(DisIdx), theta$P$Imp[DisIdx])
    Id[DisIdx]       = NextIds
    Pos = theta$Positions[DisIdx, , drop = FALSE]
    Rec[NextIds, c("Id","Bed","Room","Ward","Adm")] = list(
      NextIds, rownames(Pos), Pos$Room, Pos$Ward, t + 1
    )
    State[, t+1] = next_col
  }

  list(State = State, Id = Id, Rec = Rec)
}

sim_postprocess <- function(sim, T, theta) {
  State = sim$State
  Id    = sim$Id
  Rec   = sim$Rec

  Rec[Id, "Dis"] = T
  Rec = Rec[Rec$Id != "", ]

  # Sample mutations for all cases
  CaseRows = !is.na(Rec$Infc)
  Rec[CaseRows, "Mut"] = rpois(sum(CaseRows), theta$mu)

  CaseRec   = uplift_ancestry(Rec[CaseRows, ])
  ObsRec    = CaseRec[!is.na(CaseRec$PTest), ]

  rownames(Rec)     <- NULL
  rownames(CaseRec) <- NULL
  rownames(ObsRec)  <- NULL

  # Convert ancestry columns to integer row indices (NA = community/root case).
  # Must happen before Id is overwritten, since match() looks up against old Ids.
  CaseRec$Anc  <- match(CaseRec$Anc,  CaseRec$Id)  # "0" not in Id  →  NA
  CaseRec$Anc2 <- match(CaseRec$Anc2, CaseRec$Id)
  ObsRec$Anc2  <- match(ObsRec$Anc2,  ObsRec$Id)

  # Replace character patient IDs with sequential row positions.
  # Character so visNetwork/epicontacts IDs stay consistent with auxiliary nodes.
  Rec$Id     <- as.character(seq_len(nrow(Rec)))
  CaseRec$Id <- as.character(seq_len(nrow(CaseRec)))
  ObsRec$Id  <- as.character(seq_len(nrow(ObsRec)))

  FullTree = make_tree(CaseRec, anc_col = "Anc",  mut_col = "Mut")
  ObsTree  = make_tree(ObsRec,  anc_col = "Anc2", mut_col = "Mut2")

  list(T        = T,
       theta    = theta,
       State    = State,
       FullRec  = Rec,
       CaseRec  = CaseRec,
       ObsRec   = ObsRec,
       FullTree = FullTree,
       ObsTree  = ObsTree,
       FullDist = distances(FullTree),
       ObsDist  = distances(ObsTree))
}
#=========================================================
# Internal functions
#=========================================================
rbern <- function(n, p) rbinom(n, 1L, p)

make_tree <- function(record, anc_col, mut_col) {
  has_anc      <- !is.na(record[[anc_col]])
  n            <- nrow(record)
  tree         <- make_empty_graph(n, directed = TRUE)
  V(tree)$name <- as.character(record$Id)  # enables character-ID edge lists
  if (any(has_anc)) {
    # Use integer vertex positions for add_edges (Anc is already an integer row index)
    from_v         <- record[[anc_col]][has_anc]
    to_v           <- which(has_anc)
    tree           <- add_edges(tree, c(rbind(from_v, to_v)))
    E(tree)$weight <- record[[mut_col]][has_anc]
  }
  tree
}

## Index functions
# Get both infected and susceptible indices
get_idx = function(State) list(S = which(State == 0), I = which(State == 1))

## Infection probabilities from pre-computed log(1 - w) matrix
prob_of_inf <- function(LogComp, Idx) {
  1 - exp(rowSums(LogComp[Idx$S, Idx$I, drop = FALSE]))
}

## Sample ancestors using pre-computed odds matrix w/(1-w)
# Competing hazards: P(j sole infector) ∝ w_j * prod_{k≠j}(1-w_k) ∝ w_j/(1-w_j)
sample_ancestors <- function(Idx, NewInfIdx, Id, Odds) {
  if (length(NewInfIdx) == 0) return(character(0))
  W_odds = Odds[NewInfIdx, Idx$I, drop = FALSE]
  vapply(seq_len(nrow(W_odds)),
         function(i) sample(Id[Idx$I], 1L, prob = W_odds[i, ]),
         character(1L))
}

## Create data frame for keeping track of patients hospital records
create_record <- function(Positions, State, Id, T) {
  Size = T * nrow(State)
  Rec = data.frame(
    Id    = character(Size),
    Bed   = character(Size),
    Room  = character(Size),
    Ward  = character(Size),
    Adm   = NA_integer_,
    Dis   = NA_integer_,
    Infc  = NA_integer_,
    NTest = NA_integer_,
    PTest = NA_integer_,
    Anc   = character(Size),
    Mut   = NA_integer_
  )
  Rec[Id, c("Id", "Bed")]   = Id
  Rec[Id, c("Room","Ward")] = Positions[Id, ]
  Rec[Id, "Adm"]            = 0L

  InitInf = Id[which(State[, 1] > 0)]
  Rec[InitInf, "Infc"] = 0L
  Rec[InitInf, "Anc"]  = "0"
  Rec
}

uplift_ancestry <- function(Cases) {
  # Sort by infection time: guarantees parents are processed before children
  Cases <- Cases[order(Cases$Infc), ]
  ids   <- Cases$Id

  # Named vectors for O(1) lookup (faster than repeated data frame indexing)
  anc    <- setNames(Cases$Anc,           ids)
  mut    <- setNames(Cases$Mut,           ids)
  tested <- setNames(!is.na(Cases$PTest), ids)

  anc2 <- setNames(vector("character", length(ids)), ids)
  mut2 <- setNames(vector("integer",   length(ids)), ids)
  gen  <- setNames(vector("integer",   length(ids)), ids)

  for (i in ids) {
    p <- anc[[i]]
    if (p == "0" || tested[[p]]) {
      # Direct ancestor is root or observed: no chain to walk
      anc2[[i]] <- p
      mut2[[i]] <- mut[[i]]
      gen[[i]]  <- 1L
    } else {
      # Parent already resolved — reuse its result (memoization)
      anc2[[i]] <- anc2[[p]]
      mut2[[i]] <- mut[[i]] + mut2[[p]]
      gen[[i]]  <- gen[[p]] + 1L
    }
  }

  Cases$Anc2 <- anc2
  Cases$Mut2 <- mut2
  Cases$Gen  <- gen
  Cases
}
