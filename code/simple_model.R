P_OutsideInfection = 0.00
beta = 0.0055
NumberOfBeds = 40
NumberOfDays = 50
P_InitialInfection = 2 / NumberOfBeds
NextId = NumberOfBeds + 1
Condition = matrix(nrow = NumberOfBeds, ncol = NumberOfDays+1)
Condition[,1] = rbinom(NumberOfBeds, 1, P_InitialInfection)
Id = matrix(nrow = NumberOfBeds, ncol = NumberOfDays+1)
Id[,1] = 1:NumberOfBeds
LengthOfStay = round(rgamma(NumberOfBeds, 20, 2))
I = c()
TP = c(NumberOfBeds)
TI = c(sum(Condition))
Ancestor = c()
Ancestor[Condition[,1] == 1] = 0
NextCase = sum(Condition[,1])

for (Day in 1:NumberOfDays) {
  InfectedIndices = which(Condition[, Day] == 1)
  print("InfectedIndices")
  print(InfectedIndices)
  SusceptibleIndices = which(Condition[, Day] == 0)
  print("SusceptibleIndices")
  print(SusceptibleIndices)
  # Infect people
  NumberOfInfected = length(InfectedIndices)
  print("NumberOfInfected")
  print(NumberOfInfected)
  NumberOfSusceptible = length(SusceptibleIndices)
  print("NumberOfSusceptible")
  print(NumberOfSusceptible)
  P_Infection = 1 - (1-beta)^NumberOfInfected
  print("P_Infection")
  print(P_Infection)
  Condition[, Day + 1] = Condition[, Day]
  print("Condition0")
  print(Condition[, Day + 1])
  Condition[SusceptibleIndices, Day + 1] = rbinom(NumberOfSusceptible, 1, P_Infection)
  print("Condition1")
  print(Condition[, Day + 1])
  
  # Choose ancestors
  NewInfectedIndices = setdiff(which(Condition[, Day + 1] == 1), InfectedIndices)
  print("NewInfectedIndices")
  print(NewInfectedIndices)
  NumberOfNewInfected = length(NewInfectedIndices)
  print("NumberOfNewInfected")
  print(NumberOfNewInfected)
  print("Infected IDs")
  print(Id[InfectedIndices, Day])
  Ancestor[Id[NewInfectedIndices, Day]] = resample(Id[InfectedIndices, Day], NumberOfNewInfected, TRUE)
  print("Ancestor")
  print(Ancestor)
  
  # Infect rooms
  
  # Testing
  
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
  
  # Update statistics
  #I[day] = NumberOfInfected
  #TI[day+1] = TI[day] + sum(condition[susceptibleIndices, day+1])
}
image(Condition)
CaseIds = which(!is.na(Ancestor))
CaseAncestors = Ancestor[CaseIds]
CaseMatrix = matrix(c(CaseAncestors, CaseIds), length(CaseIds),2)
CaseGraph = graph_from_edgelist(apply(CaseMatrix[CaseMatrix[,1] != 0,], c(1,2), as.character))
plot(CaseGraph, layout = layout_as_tree)
CaseGraph = graph_from_edgelist(apply(CaseMatrix, c(1,2), as.character))
plot(CaseGraph, layout = layout_as_tree)

