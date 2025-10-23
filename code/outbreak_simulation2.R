## Function for simulating a nosocomial infectious disease outbreak
# NumBeds - number of beds in the hospital
# NumRooms - number of rooms in the hospital
# T - number of days to simulate
# P - vector of probabilities
# Beta - vector of transmission rates
# Contact - list of contact matrices
outbreak_simulation = function(NumBeds, NumRooms, T, P, beta, Contact) {
  # Setup  before simulation
  RoomNames <- letters[1:NumRooms]
  BedNames <- 1:NumBeds
  Weights = weights(beta, Contact)

  NextId = NumBeds + 1
  State = Matrix(0, nrow = NumBeds+NumRooms, ncol = T + 1)
  TestRes = Matrix(NA, nrow = NumBeds+NumRooms, ncol = T + 1)
  
  #State[1,1] = 1
  State[1:NumBeds ,1] = rbinom(NumBeds, 1, P$Init)
  
  Id = matrix(nrow = NumBeds+NumRooms, ncol = T + 1)
  Id[,1] = c(BedNames, RoomNames)
  LOS = round(rgamma(NumBeds, 20, 2))
  Size = 100
  A = data.frame(Id = character(Size),
                 Bed = character(Size),
                 Room = character(Size),
                 Ward = character(Size),
                 Adm = integer(Size),
                 Dis = integer(Size),
                 Infc = integer(Size),
                 NTest = integer(Size),
                 PTest = integer(Size),
                 Anc = character(Size))
  A[1:NumBeds, "Id"] = 1:NumBeds
  A[1:NumBeds, "Bed"] = 1:NumBeds
  A[1:NumBeds, "Room"] = c("a","a","a","a","b","b","b","b","c","c","c","c")
  A[1:NumBeds, "Ward"] = "A"
  A[1:NumBeds,"Adm"] = 0
  A[which(State[,1] > 0), "Anc"] = 0
  A[RoomNames, "Id"] = RoomNames
  A[RoomNames, "Anc"] = ""
  print(A)
  for (t in 1:T) {
    # Extract sets of indices
    Idx = get_idx(State[,t])

    # Spread infection
    State[, t + 1] = infect(State[,t], Weights, Idx)
    # Find the indices of new infectees
    NewInfIdx = new_inf_idx(State, t, Idx$I)
    # Record time of infection
    A[Id[NewInfIdx, t], "Infc"] = t
    # Sample ancestors
    A = sample_ancestors(Idx, NewInfIdx, Id[, t], Weights, A)
    # Testing
    Res = test(NumBeds, P$Test, t, State[, t], TestRes, LOS, A, Id[,t])
    TestRes = Res$TestRes
    LOS = Res$LOS
    A = Res$A
    
    # Discharge people
    LOS = LOS - 1
    Id[, t + 1] = Id[, t]
    DisIdx = which(LOS < 0)
    if (length(DisIdx) > 0) {
      NumDis = length(DisIdx)
      State[DisIdx, t + 1] = rbinom(NumDis, 1, P$Out)
      LOS[DisIdx] = round(rgamma(NumDis, 20, 2))
      A[Id[DisIdx,t], "Dis"] = t 
      Id[DisIdx, t + 1] = NextId:(NextId + NumDis - 1)
      A[NextId:(NextId + NumDis - 1), "Id"] = NextId:(NextId + NumDis - 1)
      A[NextId:(NextId + NumDis - 1), "Bed"] = -11
      A[NextId:(NextId + NumDis - 1), "Room"] = "z"
      A[NextId:(NextId + NumDis - 1), "Ward"] = "Z"
      A[NextId:(NextId + NumDis - 1), "Adm"] = t+1
      A[NextId:(NextId + NumDis - 1), "Dis"] = -1
      A[NextId:(NextId + NumDis - 1), "Infc"] = -1
      A[NextId:(NextId + NumDis - 1), "NTest"] = -1
      A[NextId:(NextId + NumDis - 1), "PTest"] = -1
      A[NextId:(NextId + NumDis - 1), "Anc"] = ""
      NextId = NextId + NumDis
    }
  }
  #print(A)
  Tree = graph_from_edgelist(as.matrix(A[A$Anc != "", c("Anc", "Id")]))
  E(Tree)$weight = rbinom(length(E(Tree)), 1000, 1/1000)
  list(State = State, Tree = Tree, Id = Id, Test = TestRes, A = A)
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
prob_of_inf <- function(Weights, Idx) {
  1 - apply(1 - Weights[Idx$S,Idx$I, drop = FALSE], 1, prod)
}

# Generating new infections
infect <- function(State, Weights, Idx) {
  P_Inf = prob_of_inf(Weights, Idx)
  State[Idx$S] = rbinom(length(Idx$S), 1, as.array(P_Inf))
  State
}


sample_ancestors <- function(Idx, NewInfIdx, Id, Weights, A) {
  if (length(NewInfIdx) > 0) {
    A[Id[NewInfIdx],"Anc"] = apply(Weights[NewInfIdx, Idx$I, drop = FALSE],
                                     1,
                                     function(w) sample(Id[Idx$I], size = 1, prob = w))
    }
  A
}

test <- function(NumBeds, P_Test, t, State, TestRes, LOS,A,Id) {
  TestIdx = which(rbinom(NumBeds, 1, P_Test) == 1)
  TestRes[TestIdx, t] = State[TestIdx]
  LOS[which(TestRes[, t] == 1)] = 0
  A[Id[which(TestRes[, t] == 0)], "NTest"] = t
  A[Id[which(TestRes[, t] == 1)], "PTest"] = t
  list(TestRes = TestRes, LOS = LOS, A = A)
}

# Weights for transmission dynamics
weights <- function(beta, Contact) {
  Reduce('+', Map(function(b, C) b * C, beta, Contact))
}