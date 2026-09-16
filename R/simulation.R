## Function for simulating a nosocomial infectious disease outbreak
# T         - number of days to simulate
# theta     - model parameters (beta, P, LogComp, Odds, mu, Positions)
simulate_outbreak = function(T, theta) {
  sim = sim_setup(T, theta)
  sim = sim_loop(sim, T, theta)
  sim_postprocess(sim, T, theta)
}
sim_setup <- function(T, theta) {
  N        <- nrow(theta$Positions)
  NumBeds  <- length(theta$P$Dis)
  NumRooms <- length(theta$P$Clean)

  State <- Matrix(0,
                  nrow = N,
                  ncol = T + 1,
                  dimnames = list(rownames(theta$Positions), NULL))
  State[, 1] <- rbern(N, theta$P$Init)

  # Bed episode IDs: bed name strings ("1", "2", ...), shifted upward on each discharge
  Id <- rownames(theta$Positions)[seq_len(NumBeds)]

  # Room episode IDs: "<roomname>_ep1", rotated on each decontamination
  RoomNames   <- rownames(theta$Positions)[NumBeds + seq_len(NumRooms)]
  RoomId      <- paste0(RoomNames, "_ep1")
  RoomEpCount <- rep(1L, NumRooms)

  Rec     <- create_patient_record(theta$Positions, State, Id, NumBeds, T)
  RoomRec <- create_room_record(theta$Positions, State, RoomId, NumBeds, NumRooms, T)

  list(State = State, Id = Id, RoomId = RoomId, RoomEpCount = RoomEpCount,
       Rec = Rec, RoomRec = RoomRec)
}

sim_loop <- function(sim, T, theta) {
  State       <- sim$State
  Id          <- sim$Id
  RoomId      <- sim$RoomId
  RoomEpCount <- sim$RoomEpCount
  Rec         <- sim$Rec
  RoomRec     <- sim$RoomRec
  NumBeds     <- length(theta$P$Dis)
  NumRooms    <- length(theta$P$Clean)

  for (t in 1:T) {
    col_t <- as.vector(State[, t])
    Idx   <- get_idx(col_t)
    AllId <- c(Id, RoomId)  # full node ID vector: beds then rooms

    # Build next state: infected carry forward, susceptibles may become infected
    P_Inf           <- prob_of_inf(theta$LogComp, Idx)
    next_col        <- col_t
    next_col[Idx$S] <- rbern(length(Idx$S), P_Inf)

    # Record new infections — split into patient vs room events
    NewInfIdx     <- which(next_col > col_t)
    NewPatInfIdx  <- NewInfIdx[NewInfIdx <= NumBeds]
    NewRoomInfIdx <- NewInfIdx[NewInfIdx >  NumBeds]

    if (length(NewPatInfIdx) > 0) {
      NewIds              <- Id[NewPatInfIdx]
      Rec[NewIds, "Infc"] <- t
      Rec[NewIds, "Anc"]  <- sample_ancestors(Idx, NewPatInfIdx, AllId, theta$Odds)
    }

    # Patient testing (beds only)
    TestIdx <- which(rbern(NumBeds, theta$P$Test) == 1)
    TestRes <- col_t[TestIdx]
    Rec[Id[TestIdx[TestRes == 0]], "NTest"] <- t
    PosIds  <- Id[TestIdx[TestRes == 1]]
    Rec[PosIds[is.na(Rec[PosIds, "PTest"])], "PTest"] <- t

    # Patient discharge and re-admission
    DisIdx  <- which(rbern(NumBeds, theta$P$Dis) == 1)
    NextIds <- as.numeric(Id[DisIdx]) + NumBeds
    Rec[Id[DisIdx], "Dis"]  <- t
    next_col[DisIdx] <- rbern(length(DisIdx), theta$P$Imp[DisIdx])
    Id[DisIdx]       <- NextIds
    Pos <- theta$Positions[DisIdx, , drop = FALSE]
    Rec[NextIds, c("Id","Bed","Room","Ward","Adm")] <- list(
      NextIds, rownames(Pos), Pos$Room, Pos$Ward, t + 1
    )

    # Room dynamics — only active when use_room_transmission = TRUE.
    # Gated here so that disabling rooms makes the simulation bit-for-bit
    # identical to the original (no extra RNG calls are made).
    if (isTRUE(theta$use_room_transmission)) {

      # Record newly infected rooms
      if (length(NewRoomInfIdx) > 0) {
        NewRoomLocal                <- NewRoomInfIdx - NumBeds
        NewRoomIds                  <- RoomId[NewRoomLocal]
        RoomRec[NewRoomIds, "Infc"] <- t
        RoomRec[NewRoomIds, "Anc"]  <- sample_ancestors(Idx, NewRoomInfIdx, AllId, theta$Odds)
      }

      # Room testing (environmental swabs) — must run before decontamination
      # rotates episode IDs below, so that a swab taken while a room is still
      # contaminated is attributed to that (still-open) episode, not to the
      # fresh episode the room rotates into if it is also cleaned today.
      RoomTestLocalIdx <- which(rbern(NumRooms, theta$P$RoomTest) == 1L)
      if (length(RoomTestLocalIdx) > 0) {
        RoomTestRes <- col_t[NumBeds + RoomTestLocalIdx]
        NegRoomIds  <- RoomId[RoomTestLocalIdx[RoomTestRes == 0L]]
        RoomRec[NegRoomIds[is.na(RoomRec[NegRoomIds, "NTest"])], "NTest"] <- t
        PosRoomIds  <- RoomId[RoomTestLocalIdx[RoomTestRes == 1L]]
        RoomRec[PosRoomIds[is.na(RoomRec[PosRoomIds, "PTest"])], "PTest"] <- t
      }

      # Room decontamination (cleaning) — only contaminated rooms are eligible
      RoomStateLocal <- col_t[NumBeds + seq_len(NumRooms)]
      ContamLocalIdx <- which(RoomStateLocal == 1L)
      if (length(ContamLocalIdx) > 0) {
        CleanLocalIdx <- ContamLocalIdx[
          rbern(length(ContamLocalIdx), theta$P$Clean[ContamLocalIdx]) == 1L
        ]
        if (length(CleanLocalIdx) > 0) {
          CleanRoomIds <- RoomId[CleanLocalIdx]
          RoomRec[CleanRoomIds, "Clean"] <- t
          next_col[NumBeds + CleanLocalIdx] <- 0L

          # Rotate episode IDs for cleaned rooms (ready for re-contamination).
          # Adm of the new episode = t: the room is susceptible from this day.
          RoomEpCount[CleanLocalIdx] <- RoomEpCount[CleanLocalIdx] + 1L
          RoomNames_now <- rownames(theta$Positions)[NumBeds + CleanLocalIdx]
          NewEpIds      <- paste0(RoomNames_now, "_ep", RoomEpCount[CleanLocalIdx])
          RoomRec[NewEpIds, c("Id","Room","Ward","Adm")] <- list(
            NewEpIds,
            theta$Positions$Room[NumBeds + CleanLocalIdx],
            theta$Positions$Ward[NumBeds + CleanLocalIdx],
            t
          )
          RoomId[CleanLocalIdx] <- NewEpIds
        }
      }

    } # end use_room_transmission

    State[, t + 1] <- next_col
  }

  list(State = State, Id = Id, RoomId = RoomId, Rec = Rec, RoomRec = RoomRec)
}

sim_postprocess <- function(sim, T, theta) {
  State   <- sim$State
  Id      <- sim$Id
  RoomId  <- sim$RoomId
  Rec     <- sim$Rec
  RoomRec <- sim$RoomRec

  # Close open patient spells still admitted at end of simulation
  Rec[Id, "Dis"] <- T
  Rec <- Rec[Rec$Id != "", ]

  # Close open room episodes still contaminated at end of simulation
  RoomRec[RoomId, "Clean"] <- T
  RoomRec <- RoomRec[RoomRec$Id != "", ]

  # Sample mutations for all patient cases
  CaseRows <- !is.na(Rec$Infc)
  Rec[CaseRows, "Mut"] <- rpois(sum(CaseRows), theta$mu)

  # Room contamination records (episodes with an infection event)
  CaseRoomRec <- RoomRec[!is.na(RoomRec$Infc), , drop = FALSE]

  # Patient ancestry: rooms treated as unobserved intermediaries so Anc2 always
  # resolves to a patient (or "0"). Raw room ancestry is preserved in RoomRec.
  CaseRec <- uplift_ancestry(Rec[CaseRows, ], room_cases = CaseRoomRec)
  ObsRec  <- CaseRec[!is.na(CaseRec$PTest), ]

  # Room-aware ground truth for scoring mcmc(inference_mode = "patients_and_rooms")
  # reconstructions. Unlike Anc2, this lets a resolved ancestor be a room episode,
  # matching what the room-aware MCMC's own ancestor field can point to. Two
  # variants mirror the two `use_room_tests` settings:
  #   Anc3     — only swab-positive rooms count as observed (use_room_tests = TRUE)
  #   Anc3_all — every contaminated room episode counts as observed (use_room_tests = FALSE)
  # Both are integer indices into c(CaseRec$Id, CaseRoomRec$Id) — patients first
  # (CaseRec's own order), then rooms (CaseRoomRec's own order) — a stable
  # reference frame independent of any later filtering to ObsRec/ObsRoomRec.
  reindex_by_id <- function(target_ids, source_ids, values) values[match(target_ids, source_ids)]

  swab_room_ids <- CaseRoomRec$Id[!is.na(CaseRoomRec$PTest)]
  up_swab <- uplift_ancestry(Rec[CaseRows, ], room_cases = CaseRoomRec,
                             room_observed_ids = swab_room_ids, return_rooms = TRUE)
  up_all  <- uplift_ancestry(Rec[CaseRows, ], room_cases = CaseRoomRec,
                             room_observed_ids = CaseRoomRec$Id, return_rooms = TRUE)

  combined_ids <- c(CaseRec$Id, CaseRoomRec$Id)  # still raw strings at this point

  CaseRec$Anc3     <- match(reindex_by_id(CaseRec$Id, up_swab$patients$Id, up_swab$patients$Anc2), combined_ids)
  CaseRec$Mut3     <- reindex_by_id(CaseRec$Id, up_swab$patients$Id, up_swab$patients$Mut2)
  CaseRec$Anc3_all <- match(reindex_by_id(CaseRec$Id, up_all$patients$Id,  up_all$patients$Anc2),  combined_ids)
  CaseRec$Mut3_all <- reindex_by_id(CaseRec$Id, up_all$patients$Id,  up_all$patients$Mut2)

  CaseRoomRec$Anc3     <- match(reindex_by_id(CaseRoomRec$Id, up_swab$rooms$Id, up_swab$rooms$Anc2), combined_ids)
  CaseRoomRec$Mut3     <- reindex_by_id(CaseRoomRec$Id, up_swab$rooms$Id, up_swab$rooms$Mut2)
  CaseRoomRec$Anc3_all <- match(reindex_by_id(CaseRoomRec$Id, up_all$rooms$Id,  up_all$rooms$Anc2),  combined_ids)
  CaseRoomRec$Mut3_all <- reindex_by_id(CaseRoomRec$Id, up_all$rooms$Id,  up_all$rooms$Mut2)

  ObsRoomRec <- CaseRoomRec[!is.na(CaseRoomRec$PTest), , drop = FALSE]
  ObsRec$Anc3     <- CaseRec$Anc3[match(ObsRec$Id, CaseRec$Id)]
  ObsRec$Anc3_all <- CaseRec$Anc3_all[match(ObsRec$Id, CaseRec$Id)]

  rownames(Rec)         <- NULL
  rownames(CaseRec)     <- NULL
  rownames(ObsRec)      <- NULL
  rownames(CaseRoomRec) <- NULL
  rownames(ObsRoomRec)  <- NULL

  # Convert ancestry columns to integer row indices (NA = community/root case).
  # Anc for patients with room ancestors → NA (room not in patient Id set),
  # which is consistent with treating the room link as latent. Anc2 is always
  # a patient index because uplift_ancestry resolves through rooms.
  CaseRec$Anc  <- match(CaseRec$Anc,  CaseRec$Id)
  CaseRec$Anc2 <- match(CaseRec$Anc2, CaseRec$Id)
  ObsRec$Anc2  <- match(ObsRec$Anc2,  ObsRec$Id)

  # Replace character patient IDs with sequential row positions.
  Rec$Id         <- as.character(seq_len(nrow(Rec)))
  CaseRec$Id     <- as.character(seq_len(nrow(CaseRec)))
  ObsRec$Id      <- as.character(seq_len(nrow(ObsRec)))
  CaseRoomRec$Id <- as.character(seq_len(nrow(CaseRoomRec)))
  ObsRoomRec$Id  <- as.character(seq_len(nrow(ObsRoomRec)))

  FullTree <- make_tree(CaseRec, anc_col = "Anc",  mut_col = "Mut")
  ObsTree  <- make_tree(ObsRec,  anc_col = "Anc2", mut_col = "Mut2")

  list(T            = T,
       theta        = theta,
       State        = State,
       FullRec      = Rec,
       CaseRec      = CaseRec,
       ObsRec       = ObsRec,
       RoomRec      = RoomRec,
       CaseRoomRec  = CaseRoomRec,
       ObsRoomRec   = ObsRoomRec,
       FullTree     = FullTree,
       ObsTree      = ObsTree,
       FullDist     = distances(FullTree),
       ObsDist      = distances(ObsTree))
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

## Patient-only record: thin wrapper over create_record for beds 1:NumBeds only
create_patient_record <- function(Positions, State, Id, NumBeds, T) {
  create_record(
    Positions[seq_len(NumBeds), , drop = FALSE],
    State[seq_len(NumBeds),     , drop = FALSE],
    Id, T
  )
}

## Room contamination episode record
# Each contamination episode per room gets one row. Pre-allocates T * NumRooms rows
# (generous upper bound: a room cleaned and re-contaminated every day).
create_room_record <- function(Positions, State, RoomId, NumBeds, NumRooms, T) {
  Size <- T * NumRooms
  Rec  <- data.frame(
    Id    = character(Size),
    Room  = character(Size),
    Ward  = character(Size),
    Adm   = 0L,           # time room last became susceptible (0 = simulation start)
    Infc  = NA_integer_,
    Clean = NA_integer_,
    NTest = NA_integer_,
    PTest = NA_integer_,
    Anc   = character(Size),
    stringsAsFactors = FALSE
  )
  RoomPos <- Positions[NumBeds + seq_len(NumRooms), , drop = FALSE]
  Rec[RoomId, c("Id","Room","Ward","Adm")] <- list(
    RoomId, RoomPos$Room, RoomPos$Ward, 0L
  )

  # Handle rooms initialised as contaminated (Adm stays 0 = simulation start)
  InitContam <- which(State[NumBeds + seq_len(NumRooms), 1] > 0)
  if (length(InitContam) > 0) {
    InitIds              <- RoomId[InitContam]
    Rec[InitIds, "Infc"] <- 0L
    Rec[InitIds, "Anc"]  <- "0"
  }
  Rec
}

## Walk the ancestry chain, skipping unobserved intermediaries (including rooms).
# room_cases: optional data frame of room contamination episodes.
# room_observed_ids: room episode Ids to treat as valid stopping points (like a
# tested patient); default none, so every room is latent and chains pass through
# them to the patient source — this is what produces Anc2. Passing specific room
# ids instead builds a room-aware ground truth (used for Anc3, see sim_postprocess)
# for scoring "patients_and_rooms" MCMC reconstructions, whose own ancestor field
# can legitimately point at a room.
# return_rooms: if TRUE, also resolve room_cases itself and return
# list(patients = Cases, rooms = room_cases) instead of just Cases.
uplift_ancestry <- function(Cases, room_cases = NULL, room_observed_ids = character(0),
                            return_rooms = FALSE) {
  # Sort by infection time: guarantees parents are processed before children
  Cases <- Cases[order(Cases$Infc), ]
  ids   <- Cases$Id

  # Build combined lookup across patients and rooms
  anc    <- setNames(Cases$Anc,           ids)
  mut    <- setNames(Cases$Mut,           ids)
  # Patients are observed if they have a positive test; rooms are observed only
  # if explicitly listed in room_observed_ids (default: none, i.e. all latent).
  tested <- setNames(!is.na(Cases$PTest), ids)

  has_rooms <- !is.null(room_cases) && nrow(room_cases) > 0
  if (has_rooms) {
    rids   <- room_cases$Id
    anc    <- c(anc,  setNames(room_cases$Anc,                        rids))
    mut    <- c(mut,  setNames(rep(NA_integer_, length(rids)),         rids))
    tested <- c(tested, setNames(rids %in% room_observed_ids,          rids))
  }

  # Resolve ancestry for all nodes (patients + rooms) in infection-time order.
  # Use order() not sort() to preserve stable ordering on tied infection times,
  # matching the behaviour of Cases[order(Cases$Infc), ] in the original code.
  all_ids  <- names(anc)
  all_infc <- c(Cases$Infc,
                if (has_rooms) room_cases$Infc else integer(0L))
  names(all_infc) <- all_ids
  ordered_ids <- all_ids[order(all_infc)]

  all_anc2 <- setNames(vector("character", length(all_ids)), all_ids)
  all_mut2 <- setNames(vector("integer",   length(all_ids)), all_ids)
  all_gen  <- setNames(vector("integer",   length(all_ids)), all_ids)

  for (i in ordered_ids) {
    p <- anc[[i]]
    if (p == "0" || isTRUE(tested[p])) {
      # Root or observed ancestor — stop walking
      all_anc2[[i]] <- p
      all_mut2[[i]] <- if (!is.na(mut[[i]])) mut[[i]] else 0L
      all_gen[[i]]  <- 1L
    } else {
      # Unobserved (or room) intermediary — reuse its resolved result (memoization)
      all_anc2[[i]] <- all_anc2[[p]]
      all_mut2[[i]] <- if (!is.na(mut[[i]])) mut[[i]] + all_mut2[[p]] else all_mut2[[p]]
      all_gen[[i]]  <- all_gen[[p]] + 1L
    }
  }

  Cases$Anc2 <- all_anc2[ids]
  Cases$Mut2 <- all_mut2[ids]
  Cases$Gen  <- all_gen[ids]

  if (return_rooms) {
    if (has_rooms) {
      room_cases$Anc2 <- all_anc2[rids]
      room_cases$Mut2 <- all_mut2[rids]
      room_cases$Gen  <- all_gen[rids]
    }
    return(list(patients = Cases, rooms = room_cases))
  }
  Cases
}
