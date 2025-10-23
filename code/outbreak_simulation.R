outbreak_simulation = function(NumBeds,
                               NumRooms,
                               T,
                               P_InitialInfection,
                               P_OutsideInfection,
                               P_Test,
                               beta_B,
                               beta_BR_RB,
                               beta_R,
                               Contact_B,
                               Contact_BR,
                               Contact_R) {
  # Setup  before simulation
  RoomNames <- letters[1:NumRooms]
  BedNames <- as.character(1:NumBeds)
  
  colnames(Contact_B[[1]]) = BedNames
  rownames(Contact_B[[1]]) = BedNames
  colnames(Contact_B[[2]]) = BedNames
  rownames(Contact_B[[2]]) = BedNames
  colnames(Contact_B[[3]]) = BedNames
  rownames(Contact_B[[3]]) = BedNames
  colnames(Contact_BR) = BedNames
  rownames(Contact_BR) = RoomNames
  colnames(Contact_R) = RoomNames
  rownames(Contact_R) = RoomNames

  Weights = weights(beta_B, beta_R, beta_BR_RB, Contact_B, Contact_BR, Contact_R)
  NextId = NumBeds + 1
  State_B = Matrix(0, nrow = NumBeds, ncol = T + 1)
  TestRes = Matrix(NA, nrow = NumBeds, ncol = T + 1)
  
  State_B[1,1] = 1
  
  State_R = Matrix(0, nrow = NumRooms, ncol = T + 1)
  Id = matrix(nrow = NumBeds, ncol = T + 1)
  Id[,1] = 1:NumBeds
  LOS = round(rgamma(NumBeds, 20, 2))
  I = c()
  TP = c(NumBeds)
  TI = c(sum(State_B))
  Ancestor_B = c()
  Ancestor_B[State_B[, 1] == 1] = 0
  Ancestor_R = character(NumRooms)

  for (t in 1:T) {
    # Extract sets of indices
    B = get_idx(State_B[, t])
    R = get_idx(State_R[, t])
    
    # Infect people
    State_B[, t + 1] =
      infect_patients(beta_B,
                      beta_BR_RB,
                      Contact_B,
                      Contact_BR,
                      B,
                      State_B[, t],
                      State_R[, t])
    
    # Infect rooms
    State_R[, t + 1] =
      infect_rooms(beta_R,
                   beta_BR_RB,
                   Contact_R,
                   Contact_BR,
                   R,
                   State_R[, t],
                   State_B[, t])
    
    # Choose ancestors for patients
    Ancestor_B = sample_ancestors(State_B, t, B, Ancestor_B, Id, Weights$B, R, RoomNames)
    
    # Choose ancestors for rooms
    NewInfectedRooms = new_inf_idx(State_R, t, R$I)
    if (length(NewInfectedRooms) > 0) {
      Ancestor_R[NewInfectedRooms] = apply(Weights$R[NewInfectedRooms, c(B$I, R$I), drop = FALSE],
                                          1,
                                          function(w) sample(c(Id[B$I, t],
                                                               RoomNames[R$I]),
                                                             size = 1,
                                                             prob = w))
    }
    
    # Testing
    Res = test_patients(NumBeds, P_Test, t, State_B, TestRes, LOS)
    LOS = Res$LOS
    TestRes = Res$TestRes
    ## Testing patients
    ## Testing rooms
    
    # Discharge people
    LOS = LOS - 1
    DischargeIdx = which(LOS < 0)
    NumDischarge = length(DischargeIdx)
    TP[t + 1] = TP[t] + NumDischarge
    State_B[DischargeIdx, t + 1] = rbinom(NumDischarge, 1, P_OutsideInfection)
    LOS[DischargeIdx] = round(rgamma(NumDischarge, 20, 2))
    Id[, t + 1] = Id[, t]
    Id[DischargeIdx, t+1] = NextId:(NextId + NumDischarge - 1)
    NextId = NextId + NumDischarge
  }

  CaseIds = which(!is.na(Ancestor_B))
  CaseAncestors = Ancestor_B[CaseIds]
  RoomIds = which(Ancestor_R != "")
  RoomMatrix = matrix(c(Ancestor_R[RoomIds], RoomNames[RoomIds]), length(RoomIds),2)
  CaseMatrix = matrix(c(CaseAncestors, CaseIds), length(CaseIds),2)
  CaseGraph = graph_from_edgelist(rbind(RoomMatrix, CaseMatrix))
  E(CaseGraph)$weight = rbinom(length(E(CaseGraph)), 1000, 1/1000)
  Dist = distances(CaseGraph)
  list(State_B, State_R, CaseGraph)
}
#=========================================================
# Internal functions
#=========================================================
## Index functions
# Get both infected and susceptible indices
get_idx = function(State) { # State vector
  list(S = sus_idx(State),
       I = inf_idx(State))
}

# Indices of infected beds
inf_idx <- function(State) { # State vector
  which(State == 1)
}
# Indices of susceptible beds
sus_idx <- function(State) { # State vector
  which(State == 0)
}
# Indices of newly infected beds
new_inf_idx <- function(State, # State matrix
                        t, # Current time
                        I) { # Infected indices at current time
  setdiff(inf_idx(State[, t + 1]), I)
}
## Infection probabilities
# Patient probabilities
probability_of_infection <- function(beta,
                                     beta_PR_RP,
                                     Contact_P,
                                     Contact_PR,
                                     S,
                                     State_P,
                                     State_R) {
  1 - (1 - beta[1])^(Contact_P[[1]][S,] %*% State_P) * # Transmission from patients sharing room
      (1 - beta[2])^(Contact_P[[2]][S,] %*% State_P) * # Transmission from patients sharing ward
      (1 - beta[3])^(Contact_P[[3]][S,] %*% State_P) * # Transmission from patients sharing hospital
      (1 - beta_PR_RP[1])^(t(Contact_PR[,S]) %*% State_R) # Transmission from rooms
}
# Room probabilities
probability_of_infection_room <- function(beta_R,
                                          beta_PR_RP,
                                          Contact_R,
                                          Contact_PR,
                                          SR, State_R,
                                          State_P) {
  1 - (1-beta_R)^(Contact_R[SR, ] %*% State_R) * # Transmission from room to room
      (1-beta_PR_RP[2])^(Contact_PR[SR, ] %*% State_P) # Transmission from patient to room
      
}
update_State <- function(State,
                             S,
                             P_Infection) {
  State[S] = rbinom(length(S), 1, as.array(P_Infection))
  State
}


# Generating new infections
infect_patients <- function(beta_P,
                            beta_PR_RP,
                            Contact_P,
                            Contact_PR,
                            P,
                            State_P,
                            State_R) {
  # Infect people
  P_Infection = probability_of_infection(beta_P,
                                         beta_PR_RP,
                                         Contact_P,
                                         Contact_PR,
                                         P$S,
                                         State_P,
                                         State_R)
  update_State(State_P,
                   P$S,
                   P_Infection)
}
infect_rooms <- function(beta_R,
                         beta_PR_RP,
                         Contact_R,
                         Contact_PR,
                         R,
                         State_R,
                         State_P) {
  P_Infection = probability_of_infection_room(beta_R,
                                              beta_PR_RP,
                                              Contact_R,
                                              Contact_PR,
                                              R$S,
                                              State_R,
                                              State_P)
  update_State(State_R,
                   R$S,
                   P_Infection)
}

sample_ancestors <- function(State_P, t, P, Ancestor_P, Id, Weights_P, R, RoomNames) {
  NewInfIdx = new_inf_idx(State_P, t, P$I)
  if (length(NewInfIdx) > 0) {
    Ancestor_P[Id[NewInfIdx, t]] = apply(Weights_P[NewInfIdx, c(P$I, R$I), drop = FALSE],
                                                    1,
                                                    function(w) sample(c(Id[P$I, t],
                                                                         RoomNames[R$I]),
                                                                       size = 1,
                                                                       prob = w))
  }
  Ancestor_P
}

test_patients <- function(NumberOfBeds,P_Test, t, State_P, TestRes, LOS) {
  TestIdx = which(rbinom(NumberOfBeds, 1, P_Test) == 1)
  TestRes[TestIdx, t] = State_P[TestIdx, t]
  LOS[which(TestRes[, t] == 1)] = 0
  list(TestRes = TestRes, LOS = LOS)
}
discharge_admission <- function(LOS, TP, State_B, Id, P_OutsideInfection, t, NextId) {
  LOS = LOS - 1
  DischargeIdx = which(LOS < 0)
  NumDischarge = length(DischargeIdx)
  TP[t + 1] = TP[t] + NumDischarge
  State_B[DischargeIdx, t + 1] = rbinom(NumDischarge, 1, P_OutsideInfection)
  LOS[DischargeIdx] = round(rgamma(NumDischarge, 20, 2))
  Id[, t + 1] = Id[, t]
  Id[DischargeIdx, t+1] = NextId:(NextId + NumDischarge - 1)
  NextId = NextId + NumDischarge
}
weights <- function(beta_B, beta_R, beta_BR_RB, Contact_B, Contact_BR, Contact_R) {
  list(B = weights_b(beta_B, Contact_B, beta_BR_RB, Contact_BR),
       R = weigths_r(beta_PR_RP, beta_R, Contact_BR, Contact_R))
}
weights_b <- function(beta_B, Contact_B, beta_BR_RB, Contact_BR) {
  cbind(beta_B[1] * Contact_B[[1]] + # Ancestor patient sharing room
          beta_B[2] * Contact_B[[2]] + # Ancestor patient sharing ward
          beta_B[3] * Contact_B[[3]],  # Ancestor patient sharing hospital
        beta_BR_RB[1] * t(Contact_BR)) # Ancestor room
}
weigths_r <- function(beta_BR_RB, beta_R, Contact_BR, Contact_R) {
  cbind(beta_BR_RB[2] * Contact_BR, # Ancestor patient in room
        beta_R * Contact_R) # Ancestor another room
}