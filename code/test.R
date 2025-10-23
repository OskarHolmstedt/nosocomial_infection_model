library(Matrix)
library(igraph)
source("small_ward_data.R")
#source("large_ward_data.R")
#source("two_wards_data.R")
source("parameters.R")
source("outbreak_simulation2.R")

RoomNames <- letters[1:NumRooms]
BedNames <- as.character(1:NumBeds)
beta = c(beta_B, beta_BR_RB, beta_R)
Ctct = test_contact_matrices(NumBeds,
                      NumRooms,
                      BedNames,
                      RoomNames,
                      Contact,
                      ContactPR,
                      ContactRR)
Wgs = weights(beta, Ctct)
Td = 5
State = Matrix(0, nrow = NumBeds+NumRooms, ncol = Td + 1)
State[1,1] = 1
Idx = get_idx(State[,1])
probability_of_infection(beta, Ctct, Idx$S, State[,1])

test_contact_matrices <- function(NumBeds,
                                  NumRooms,
                                  BedNames,
                                  RoomNames,
                                  Contact,
                                  ContactPR,
                                  ContactRR) {
  contact_matrices(NumBeds,
                   NumRooms,
                   BedNames,
                   RoomNames,
                   Contact,
                   ContactPR,
                   ContactRR)
}