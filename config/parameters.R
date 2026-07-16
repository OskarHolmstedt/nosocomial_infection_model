# Model parameters
theta <- local({
  mu = 2
  beta <- c(
    BrB = 0.08,   # bed to bed, same room
    BwB = 0.005,   # bed to bed, same ward
    BhB = 0.0005,   # bed to bed, hospital-wide
    BR  = 0.00,   # bed to room
    RB  = 0.00,   # room to bed
    RR  = 0.00    # room to room
  )
  P <- list(
    Init = c(1, rep(0, NumBeds - 1), rep(0, NumRooms)),
    Imp  = rep(0L,   NumBeds),
    Test = rep(0.1,  NumBeds),
    Dis  = rep(0.1, NumBeds)
  )
  W <- weights(beta, Contact)
  list(
    mu            = mu,
    beta          = beta,
    P             = P,
    # Geometric prior rate for imported (community) cases: P(inf) ∝ (1-lambda_import)^(inf-t_root).
    # Needs to be >= P$Test to overcome the testing-likelihood pull toward later times.
    lambda_import = 0.3,
    pi_import     = 0.05,   # prior probability that a case is community-imported
    LogComp       = log1p(-W),
    Odds          = W / (1 - W),
    spatial_sizes = spatial_sizes,
    Positions     = Positions
  )
})
