mcmc <- function(Data, N, ) {
  for (n in 1:N) {
    #propose new parameters
    #compute likelihood
    #accept or reject
  }
}

log_likelihood <- function() {
  genetic_log_likelihood() +
    testing_log_likelihood() +
    contact_log_likelihood()
}
genetic_log_likelihood() {
  # Poisson?
}
testing_log_likelihood() {
  # Geometric?
}
contact_log_likelihood() {
  # What to do?
}
