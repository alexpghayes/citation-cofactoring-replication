model4 <- function(n, k, sort_nodes, rep = 1, diagonal = 0.8, within_between_ratio = 3, off_diagonal = 0.01) {

  B <- matrix(off_diagonal, nrow = k, ncol = k)
  B_between <- (diagonal / within_between_ratio - (k - 2) * off_diagonal)
  B[1, ] <- B_between
  B[k, 2] <- B_between
  diag(B) <- diagonal

  log_debug(
    "Within community mass is {sum(diag(B))} and between community mass is {sum(B) - sum(diag(B))}."
  )

  theta_in <- 1 + rexp(n, 1/8)
  theta_out <- 1 + rexp(n, 1/8)

  m <- directed_dcsbm(
    theta_in = theta_in / sum(theta_in),
    theta_out = theta_out / sum(theta_out),
    B = B,
    expected_density = 0.15,
    sort_nodes = sort_nodes
  )
  m$rep <- rep
  m
}


model5 <- function(n, k, sort_nodes, rep = 1, diagonal = 0.8, within_between_ratio = 3, off_diagonal = 0.01) {

  B <- matrix(off_diagonal, nrow = k, ncol = k)
  B_between <- (diagonal / within_between_ratio - (k - 2) * off_diagonal)
  B[1, ] <- B_between
  B[k, 2] <- B_between
  diag(B) <- diagonal
  B <- (t(B) + B) / 2

  log_debug(
    "Within community mass is {sum(diag(B))} and between community mass is {sum(B) - sum(diag(B))}."
  )

  theta_in <- 1 + rexp(n, 1/8)
  theta_out <- theta_in

  m <- directed_dcsbm(
    theta_in = theta_in / sum(theta_in),
    theta_out = theta_out / sum(theta_out),
    B = B,
    expected_density = 0.15,
    sort_nodes = sort_nodes
  )
  m$rep <- rep
  m
}

model6 <- function(n, k, sort_nodes, rep = 1, diagonal = 0.8, within_between_ratio = 3, off_diagonal = 0.01) {

  B <- matrix(off_diagonal, nrow = k, ncol = k)
  B_between <- (diagonal / within_between_ratio - (k - 2) * off_diagonal)
  B[1, ] <- B_between
  B[k, 2] <- B_between
  diag(B) <- diagonal
  B <- (t(B) + B) / 2

  log_debug(
    "Within community mass is {sum(diag(B))} and between community mass is {sum(B) - sum(diag(B))}."
  )

  theta <- 1 + rexp(n, 1/8)

  m <- dcsbm(
    theta = theta,
    B = B,
    expected_density = 0.15,
    sort_nodes = sort_nodes
  )
  m$rep <- rep
  m
}
