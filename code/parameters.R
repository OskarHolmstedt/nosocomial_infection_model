# Parameters

P_OutsideInfection = 0.00
P_Test = 1/NumBeds
beta_BrB = 0.08
beta_BwB = 0.02
beta_BhB = 0.01
beta_B = c(beta_BrB, beta_BwB, beta_BhB)
beta_BR = 0.01
beta_RB = 0.01
beta_BR_RB = c(beta_BR, beta_RB)
beta_RR = 0
beta_R = c(beta_RR)
beta = c(beta_B, beta_BR_RB, beta_R)
NumDays = 50
P_InitialInfection = 1 / NumBeds
P = list(Init = P_InitialInfection, Out = P_OutsideInfection, Test = P_Test)