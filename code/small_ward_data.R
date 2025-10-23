# Small ward
# Three rooms, 4 beds each
NumBeds = 12
NumRooms = 3

qp = as_adjacency_matrix(make_full_graph(4))
zp = matrix(0,4,4)
op = matrix(1,4,4)
ContactR = cbind(rbind(qp,zp,zp),rbind(zp,qp,zp),rbind(zp,zp,qp))

ContactW = cbind(rbind(zp,op,op),rbind(op,zp,op),rbind(op,op,zp))

ContactH = Matrix(1, NumBeds, NumBeds, sparse = TRUE) - diag(NumBeds) - ContactR - ContactW
ContactB = list(ContactR, ContactW, ContactH)

ContactBR = t(Matrix(c(1,1,1,1,0,0,0,0,0,0,0,0,
                       0,0,0,0,1,1,1,1,0,0,0,0,
                       0,0,0,0,0,0,0,0,1,1,1,1),
                     NumBeds, NumRooms))

ContactRR = Matrix(0, NumRooms, NumRooms)

RoomNames <- letters[1:NumRooms]
BedNames <- as.character(1:NumBeds)
Contact = contact_matrices(NumBeds, NumRooms, BedNames, RoomNames, ContactB, ContactBR, ContactRR)
