# in practice we will pass A, rank, max_iter, and additional to all of these
# functions and they will ignore any arguments they don't need

adaptive_initialize <- function(A, rank, ...) {

  fastadi::adaptive_initialize(
    A, rank = rank,
    p_hat = 0.5,
    alpha_method = "approximate",
    ...
  )
}

zero_imputed_svd <- function(A, rank, ...) {
  s <- svds(A, rank)
  as_svd_like(s)
}

citation_impute_svd <- function(A, rank, ...) {
  citation_impute2(
    A, rank = rank,
    check_interval = NULL,
    initialization = "svd",
    ...
  )
}

citation_impute_approximate <- function(A, rank, ...) {
  citation_impute2(
    A, rank = rank,
    check_interval = NULL,
    initialization = "approximate",
    ...
  )
}

citation_impute_preclipped <- function(A, rank, ell_z, ell_y, ...) {

  n <- nrow(A)

  # here we set rows and columns of a sparse matrix to zero without
  # adding zeros explicitly, for memory efficiency

  # set the first ell_y columns of A to zero, and the last ell_z rows of
  # A to zero. if A is symmetric, the storage format is different
  # and things will go sideways here

  m <- as(A, "TsparseMatrix")

  ind <- m@i < (n - ell_z) & m@j >= ell_y
  m@x <- m@x[ind]
  m@i <- m@i[ind]
  m@j <- m@j[ind]

  A <- as(m, "CsparseMatrix")

  citation_impute2(
    A, rank = rank,
    check_interval = NULL,
    initialization = "approximate",
    ...
  )
}


citation_impute_preclipped2 <- function(A, rank, ell_z = 100000, ell_y = 50000, ...) {

  n <- nrow(A)

  # here we set rows and columns of a sparse matrix to zero without
  # adding zeros explicitly, for memory efficiency

  # set the first ell_y columns of A to zero, and the last ell_z rows of
  # A to zero. if A is symmetric, the storage format is different
  # and things will go sideways here

  m <- as(A, "TsparseMatrix")

  ind <- m@i < (n - ell_z) & m@j >= ell_y
  m@x <- m@x[ind]
  m@i <- m@i[ind]
  m@j <- m@j[ind]

  A <- as(m, "CsparseMatrix")

  citation_impute2(
    A, rank = rank,
    check_interval = NULL,
    initialization = "approximate",
    ...
  )
}

clip_vsp <- function(fa, ell_y, ell_z) {

  n <- nrow(fa$u)

  fa$u[(n - ell_z):n, ] <- 0
  fa$Z[(n - ell_z):n, ] <- 0

  fa$v[1:ell_y, ] <- 0
  fa$Y[1:ell_y, ] <- 0

  fa
}
