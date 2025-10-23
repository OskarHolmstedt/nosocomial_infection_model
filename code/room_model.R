library(Matrix)
library(igraph)
NumberOfBeds = 12
NumberOfRooms = 3
RoomNames <- letters[1:NumberOfRooms]
BedNames <- as.character(1:NumberOfBeds)
qp = as_adjacency_matrix(make_full_graph(4))
zp = matrix(0,4,4)
op = matrix(1,4,4)
ContactR = cbind(rbind(qp,zp,zp),rbind(zp,qp,zp),rbind(zp,zp,qp))
colnames(ContactR) = BedNames
rownames(ContactR) = BedNames
ContactW = cbind(rbind(zp,op,op),rbind(op,zp,op),rbind(op,op,zp))
colnames(ContactW) = BedNames
rownames(ContactW) = BedNames
ContactH = Matrix(1, NumberOfBeds, NumberOfBeds) - diag(NumberOfBeds) - ContactR - ContactW
colnames(ContactH) = BedNames
rownames(ContactH) = BedNames
ContactPR = t(Matrix(c(1,1,1,1,0,0,0,0,0,0,0,0,
                       0,0,0,0,1,1,1,1,0,0,0,0,
                       0,0,0,0,0,0,0,0,1,1,1,1),
                     NumberOfBeds, NumberOfRooms))
colnames(ContactPR) = BedNames
rownames(ContactPR) = RoomNames
ContactRR = Matrix(0, NumberOfRooms, NumberOfRooms)
colnames(ContactRR) = RoomNames
rownames(ContactRR) = RoomNames

P_OutsideInfection = 0.00
P_Test = 1/NumberOfBeds
betaPRP = 0.08
betaPWP = 0.02
betaPHP = 0.00
betaPR = 0.01
betaRP = 0.01
betaRR = 0
WeightsP = cbind(betaR * ContactR + # Ancestor patient sharing room
                 betaW * ContactW + # Ancestor patient sharing ward
                 betaH * ContactH,  # Ancestor patient sharing hospital
                 betaPR * t(ContactPR)) # Ancestor room
WeightsR = cbind(betaRP * ContactPR, # Ancestor patient in room
                 betaRR * ContactRR) # Ancestor another room
NumberOfDays = 50
P_InitialInfection = 1 / NumberOfBeds
NextId = NumberOfBeds + 1
Condition = Matrix(0, nrow = NumberOfBeds, ncol = NumberOfDays+1)
TestResults = Matrix(NA, nrow = NumberOfBeds, ncol = NumberOfDays+1)

Condition[,1] = c(1,0,0,0, 0,0,0,0, 0,0,0,0)

ConditionR = Matrix(0, nrow = NumberOfRooms, ncol = NumberOfDays+1)
Id = matrix(nrow = NumberOfBeds, ncol = NumberOfDays+1)
Id[,1] = 1:NumberOfBeds
LengthOfStay = round(rgamma(NumberOfBeds, 20, 2))
I = c()
TP = c(NumberOfBeds)
TI = c(sum(Condition))
Ancestor = c()
Ancestor[Condition[,1] == 1] = 0
AncestorR = character(NumberOfRooms)
NextCase = sum(Condition[,1])

for (Day in 1:NumberOfDays) {
  InfectedIndices = which(Condition[, Day] == 1)
  SusceptibleIndices = which(Condition[, Day] == 0)
  InfectedRooms = which(ConditionR[,Day] == 1)
  SusceptibleRooms = which(ConditionR[,Day] == 0)

  # Infect people
  P_Infection = 1 - (1 - betaR)^(ContactR[SusceptibleIndices,] %*% Condition[, Day]) * # Transmission from patients sharing room
                    (1 - betaW)^(ContactW[SusceptibleIndices,] %*% Condition[, Day]) * # Transmission from patients sharing ward
                    (1 - betaH)^(ContactH[SusceptibleIndices,] %*% Condition[, Day]) * # Transmission from patients sharing ward
                    (1 - betaRP)^(t(ContactPR[,SusceptibleIndices]) %*% ConditionR[, Day]) # Transmission from rooms
  Condition[InfectedIndices, Day + 1] = Condition[InfectedIndices, Day]
  Condition[SusceptibleIndices, Day + 1] = rbinom(length(SusceptibleIndices), 1, as.array(P_Infection))

  # Infect rooms
  P_InfectionR = 1 - (1-betaPR)^(ContactPR[SusceptibleRooms, ]%*%Condition[,Day]) *
                     (1-betaRR)^(ContactRR[SusceptibleRooms, ]%*%ConditionR[,Day])
  ConditionR[InfectedRooms, Day + 1] = ConditionR[InfectedRooms, Day]
  ConditionR[SusceptibleRooms, Day + 1] = rbinom(length(SusceptibleRooms), 1, as.array(P_InfectionR))
  
  # Choose ancestors for patients
  NewInfectedIndices = setdiff(which(Condition[, Day + 1] == 1), InfectedIndices)
  if (length(NewInfectedIndices) > 0) {
    Ancestor[Id[NewInfectedIndices, Day]] = apply(WeightsP[NewInfectedIndices, c(InfectedIndices, InfectedRooms), drop = FALSE],
                                                  1,
                                                  function(w) sample(c(Id[InfectedIndices, Day],
                                                                       RoomNames[InfectedRooms]),
                                                                     size = 1,
                                                                     prob = w))
  }
  # Choose ancestors for rooms
  NewInfectedRooms = setdiff(which(ConditionR[, Day + 1] == 1), InfectedRooms)
  print(NewInfectedRooms)
  if (length(NewInfectedRooms) > 0) {
    AncestorR[NewInfectedRooms] = apply(WeightsR[NewInfectedRooms, c(InfectedIndices, InfectedRooms), drop = FALSE],
                                        1,
                                        function(w) sample(c(Id[InfectedIndices, Day],
                                                             RoomNames[InfectedRooms]),
                                                             size = 1,
                                                             prob = w))
  }
  
  # Testing
  TestIndices = which(rbinom(NumberOfBeds, 1, P_Test) == 1)
  TestResults[TestIndices, Day] = Condition[TestIndices, Day]
  LengthOfStay[which(TestResults[,Day] == 1)] = 0
  ## Testing patients
  ## Testing rooms
  
  # Discharge people
  LengthOfStay = LengthOfStay - 1
  DischargeIndices = which(LengthOfStay < 0)
  NumberOfDischarges = length(DischargeIndices)
  TP[Day + 1] = TP[Day] + NumberOfDischarges
  Condition[DischargeIndices, Day + 1] = rbinom(NumberOfDischarges, 1, P_OutsideInfection)
  LengthOfStay[DischargeIndices] = round(rgamma(NumberOfDischarges, 20, 2))
  Id[, Day + 1] = Id[, Day]
  Id[DischargeIndices, Day+1] = NextId:(NextId + NumberOfDischarges - 1)
  NextId = NextId + NumberOfDischarges
}
image(rbind(Condition, 2*ConditionR), xlab = "Days", ylab = "Patients/Rooms")
print(ConditionR)
CaseIds = which(!is.na(Ancestor))
CaseAncestors = Ancestor[CaseIds]
RoomIds = which(AncestorR != "")
RoomMatrix = matrix(c(AncestorR[RoomIds], RoomNames[RoomIds]), length(RoomIds),2)
CaseMatrix = matrix(c(CaseAncestors, CaseIds), length(CaseIds),2)
CaseGraph = graph_from_edgelist(rbind(RoomMatrix, CaseMatrix))
E(CaseGraph)$weight = rbinom(length(E(CaseGraph)), 1000, 1/1000)
plot(CaseGraph, layout = layout_as_tree, edge.label = E(CaseGraph)$weight)
Dist = distances(CaseGraph)
