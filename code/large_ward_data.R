# Large ward
# 6 rooms, 3 beds each
NumberOfBeds = 18
NumberOfRooms = 6

qp = as_adjacency_matrix(make_full_graph(3))
zp = matrix(0,3,3)
op = matrix(1,3,3)
ContactR = cbind(rbind(qp,zp,zp,zp,zp,zp),
                 rbind(zp,qp,zp,zp,zp,zp),
                 rbind(zp,zp,qp,zp,zp,zp),
                 rbind(zp,zp,zp,qp,zp,zp),
                 rbind(zp,zp,zp,zp,qp,zp),
                 rbind(zp,zp,zp,zp,zp,qp))
ContactW = cbind(rbind(zp,op,op,op,op,op),
                 rbind(op,zp,op,op,op,op),
                 rbind(op,op,zp,op,op,op),
                 rbind(op,op,op,zp,op,op),
                 rbind(op,op,op,op,zp,op),
                 rbind(op,op,op,op,op,zp))

ContactH = Matrix(1, NumberOfBeds, NumberOfBeds, sparse = TRUE) - diag(NumberOfBeds) - ContactR - ContactW
Contact = list(ContactR, ContactW, ContactH)

ContactPR = t(Matrix(c(1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
                       0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
                       0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0,
                       0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0,
                       0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0,
                       0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1),
                     NumberOfBeds, NumberOfRooms))

ContactRR = Matrix(0, NumberOfRooms, NumberOfRooms)