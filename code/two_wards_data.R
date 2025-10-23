# Two small wards
# Two wards, Three rooms, 2 beds each
NumberOfBeds = 12
NumberOfRooms = 6

qp = as_adjacency_matrix(make_full_graph(2))
zp = matrix(0,2,2)
op = matrix(1,2,2)
ContactR = cbind(rbind(qp,zp,zp,zp,zp,zp),
                 rbind(zp,qp,zp,zp,zp,zp),
                 rbind(zp,zp,qp,zp,zp,zp),
                 rbind(zp,zp,zp,qp,zp,zp),
                 rbind(zp,zp,zp,zp,qp,zp),
                 rbind(zp,zp,zp,zp,zp,qp))

ContactW = cbind(rbind(zp,op,op,zp,zp,zp),
                 rbind(op,zp,op,zp,zp,zp),
                 rbind(op,op,zp,zp,zp,zp),
                 rbind(zp,zp,zp,zp,op,op),
                 rbind(zp,zp,zp,op,zp,op),
                 rbind(zp,zp,zp,op,op,zp))

ContactH = Matrix(1, NumberOfBeds, NumberOfBeds, sparse = TRUE) - diag(NumberOfBeds) - ContactR - ContactW
Contact = list(ContactR, ContactW, ContactH)

ContactPR = t(Matrix(c(1,1,0,0,0,0,0,0,0,0,0,0,
                       0,0,1,1,0,0,0,0,0,0,0,0,
                       0,0,0,0,1,1,0,0,0,0,0,0,
                       0,0,0,0,0,0,1,1,0,0,0,0,
                       0,0,0,0,0,0,0,0,1,1,0,0,
                       0,0,0,0,0,0,0,0,0,0,1,1),
                     NumberOfBeds, NumberOfRooms))

ContactRR = Matrix(0, NumberOfRooms, NumberOfRooms)