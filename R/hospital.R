# Weights for transmission dynamics
weights <- function(beta, Contact) {
  Reduce('+', Map(function(b, C) b * C, beta, Contact))
}

# Build a uniform hospital layout: Wards wards, RoomsPerWard rooms each, BedsPerRoom beds each
hospital <- function(Wards, RoomsPerWard, BedsPerRoom) {
  NumBeds    <- Wards * RoomsPerWard * BedsPerRoom
  NumRooms   <- Wards * RoomsPerWard
  BedNames   <- seq_len(NumBeds)
  RoomNames  <- paste0("r", seq_len(NumRooms))
  WardNames  <- paste0("W", seq_len(Wards))

  # Bed-bed contact matrices at each spatial scale
  J_B  <- matrix(1, BedsPerRoom, BedsPerRoom)
  J_RB <- matrix(1, RoomsPerWard * BedsPerRoom, RoomsPerWard * BedsPerRoom)

  ContactR <- bdiag(replicate(NumRooms, J_B,  simplify = FALSE)) - Diagonal(NumBeds)
  ContactW <- bdiag(replicate(Wards,    J_RB, simplify = FALSE)) -
              bdiag(replicate(NumRooms, J_B,  simplify = FALSE))
  ContactH <- Matrix(1, NumBeds, NumBeds, sparse = TRUE) -
              Diagonal(NumBeds) - ContactR - ContactW

  ContactBR <- kronecker(Diagonal(NumRooms), matrix(1, 1, BedsPerRoom))

  # Embed blocks into the full (beds + rooms) space
  Z_BB <- Matrix(0, NumBeds,  NumBeds,  sparse = TRUE)
  Z_BR <- Matrix(0, NumBeds,  NumRooms, sparse = TRUE)
  Z_RB <- Matrix(0, NumRooms, NumBeds,  sparse = TRUE)
  Z_RR <- Matrix(0, NumRooms, NumRooms, sparse = TRUE)

  extend <- function(BB = Z_BB, BR = Z_BR, RB = Z_RB, RR = Z_RR)
    rbind(cbind(BB, BR), cbind(RB, RR))

  node_names <- c(as.character(BedNames), RoomNames)
  Contact <- lapply(
    list(
      extend(BB = ContactR),      # bed-bed, same room
      extend(BB = ContactW),      # bed-bed, same ward
      extend(BB = ContactH),      # bed-bed, hospital-wide
      extend(BR = t(ContactBR)),  # bed -> room
      extend(RB = ContactBR),     # room -> bed
      extend(RR = Z_RR)           # room -> room (zero)
    ),
    function(M) { dimnames(M) <- list(node_names, node_names); M }
  )

  # Positions: room and ward label for each bed and room node
  WardPositions <- c(rep(WardNames, each = RoomsPerWard * BedsPerRoom),
                     rep(WardNames, each = RoomsPerWard))
  RoomPositions <- c(RoomNames[rep(seq_len(NumRooms), each = BedsPerRoom)],
                     RoomNames)
  Positions <- data.frame(
    Room = RoomPositions,
    Ward = WardPositions,
    row.names        = c(as.character(BedNames), RoomNames),
    stringsAsFactors = FALSE
  )

  # Approximate number of patients at each spatial scale (room, ward, hospital)
  spatial_sizes <- c(BedsPerRoom,
                     RoomsPerWard * BedsPerRoom,
                     NumBeds)

  list(NumBeds       = NumBeds,
       NumRooms      = NumRooms,
       BedNames      = BedNames,
       RoomNames     = RoomNames,
       Contact       = Contact,
       Positions     = Positions,
       spatial_sizes = spatial_sizes)
}
