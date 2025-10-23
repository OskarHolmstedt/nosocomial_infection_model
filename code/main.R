library(Matrix)
library(igraph)
source("prep_functions.R")
source("small_ward_data.R")
#source("large_ward_data.R")
#source("two_wards_data.R")
source("parameters.R")
source("outbreak_simulation2.R")

#set.seed(9918)
res = outbreak_simulation(NumBeds,
                          NumRooms,
                          NumDays,
                          P,
                          beta,
                          Contact)
State = res$State
Tree = res$Tree
image(State, xlab = "Days", ylab = "Patients/Rooms")

plot(Tree, layout = layout_as_tree, edge.label = E(Tree)$weight)
Dist = distances(Tree)

