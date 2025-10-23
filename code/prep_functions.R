contact_matrices <- function(NumBeds,
                             NumRooms,
                             BedNames,
                             RoomNames,
                             Contact_B,
                             Contact_BR,
                             Contact_R) {
  Contact_B[[1]] = rbind(cbind(Contact_B[[1]],
                               matrix(0, nrow = NumBeds, ncol = NumRooms)),
                         matrix(0, nrow = NumRooms, ncol = NumBeds+NumRooms))
  Contact_B[[2]] = rbind(cbind(Contact_B[[2]],
                               matrix(0, nrow = NumBeds, ncol = NumRooms)),
                         matrix(0, nrow = NumRooms, ncol = NumBeds+NumRooms))
  Contact_B[[3]] = rbind(cbind(Contact_B[[3]],
                               matrix(0, nrow = NumBeds, ncol = NumRooms)),
                         matrix(0, nrow = NumRooms, ncol = NumBeds+NumRooms))
  Contact_B[[4]] = rbind(cbind(matrix(0, nrow = NumBeds, ncol = NumBeds), t(Contact_BR)),
                         matrix(0, nrow = NumRooms, ncol = NumBeds+NumRooms))
  Contact_B[[5]] = cbind(rbind(matrix(0, nrow = NumBeds, ncol = NumBeds), Contact_BR),
                         matrix(0, nrow = NumBeds+NumRooms, ncol = NumRooms))
  Contact_B[[6]] = cbind(matrix(0, nrow = NumBeds+NumRooms, ncol = NumBeds),
                         rbind(matrix(0, nrow = NumBeds, ncol = NumRooms),
                               Contact_R))
  Contact = Contact_B
  colnames(Contact[[1]]) = c(BedNames, RoomNames)
  rownames(Contact[[1]]) = c(BedNames, RoomNames)
  colnames(Contact[[2]]) = c(BedNames, RoomNames)
  rownames(Contact[[2]]) = c(BedNames, RoomNames)
  colnames(Contact[[3]]) = c(BedNames, RoomNames)
  rownames(Contact[[3]]) = c(BedNames, RoomNames)
  colnames(Contact[[4]]) = c(BedNames, RoomNames)
  rownames(Contact[[4]]) = c(BedNames, RoomNames)
  colnames(Contact[[5]]) = c(BedNames, RoomNames)
  rownames(Contact[[5]]) = c(BedNames, RoomNames)
  colnames(Contact[[6]]) = c(BedNames, RoomNames)
  rownames(Contact[[6]]) = c(BedNames, RoomNames)
  Contact
}
