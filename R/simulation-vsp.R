compute_alignment <- function(true_fa, fa) {

  # see https://cran.r-project.org/web/packages/RcppHungarian/vignettes/Introduction-to-RcppHungarian.html

  # the issue is that varimax estimates are only identified up to sign
  # and column order permutation. if you naively search over column
  # order permutations you get hit with O(k!) runtime. luckily column
  # ordering reduces to the assignment problem, which can be solved
  # much more efficiently with the Hungarian algorithm

  k <- fa$rank

  Z_true <- true_fa$Z
  Z_est <- cbind(fa$Z, -fa$Z)

  Y_true <- true_fa$Y
  Y_est <- cbind(fa$Y, -fa$Y)

  Cz <- matrix(0, k, 2 * k)  # to align Z factors
  Cy <- matrix(0, k, 2 * k)  # to align Y factors

  for (pop in 1:k) {
    for (est in 1:(2 * k)) {
      Cz[pop, est] <- sum((Z_true[, pop] - Z_est[, est])^2)
      Cy[pop, est] <- sum((Y_true[, pop] - Y_est[, est])^2)
    }
  }

  # find the column matching that minimizes Frobenius loss of the
  # vsp estimates

  list(
    z_solved = HungarianSolver(Cz),
    y_solved = HungarianSolver(Cy),
    Cz = Cz,
    Cy = Cy,
    Z_est = Z_est,
    Y_est = Y_est
  )
}

align_factors <- function(fa, alignment) {

  z_est_index <- alignment$z_solved$pairs[, 2]
  y_est_index <- alignment$y_solved$pairs[, 2]

  fa$Z <- alignment$Z_est[, z_est_index]
  fa$Y <- alignment$Y_est[, y_est_index]

  # also need to align B
  for (index in seq_along(z_est_index)) {
    if (z_est_index[index] > fa$rank) {
      z_est_index[index] <- z_est_index[index] - fa$rank
    }
    if (y_est_index[index] > fa$rank) {
      y_est_index[index] <- y_est_index[index] - fa$rank
    }
  }
  fa$B <- fa$B[z_est_index, y_est_index]
  fa
}
