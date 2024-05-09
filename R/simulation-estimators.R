clip <- function(M, l) {

  # M[, 1:l] <- 0
  # M[(nrow(M) - l + 1):nrow(M), ] <- 0
  # M

  ell_y <- round(l)
  ell_z <- round(l)

  n <- nrow(M)
  m <- as(M, "TsparseMatrix")

  ind <- m@i < (n - ell_z) & m@j >= ell_y
  m@x <- m@x[ind]
  m@i <- m@i[ind]
  m@j <- m@j[ind]

  as(m, "CsparseMatrix")
}

full_svds <- function(A, k, ...) {
  s <- svds(A, k)
  as_svd_like(s)
}

zero_imputed_svds <- function(A, k, ...) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  s <- svds(A_obs, k)
  as_svd_like(s)
}

symmetric_svd <- function(A, k, ...) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  A_sym <- (A_obs + t(A_obs)) / 2
  s <- svds(A_sym, k)
  as_svd_like(s)
}

adapt_impute <- function(A, k, max_iter, ...) {
  n <- nrow(A)
  mask <- rsparsematrix(n, n, nnz = n^2 / 2, rand.x = NULL)
  A_obs <- A * mask
  fastadi::adaptive_impute(
    A_obs,
    rank = k,
    max_iter = max_iter,
    initialization = "approximate",
    additional = 20
  )
}

cite_impute <- function(A, k, max_iter, ...) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  fastadi::citation_impute2(
    clip(A_obs, nrow(A) / 10),
    rank = k,
    max_iter = max_iter,
    initialization = "approximate",
    additional = 20
  )
}

cite_impute2 <- function(A, k, max_iter, ...) {
  A_obs <- as(triu(A) * 1, "dgCMatrix")
  fastadi::citation_impute2(
    A_obs,
    rank = k,
    max_iter = max_iter,
    initialization = "approximate",
    additional = 20
  )
}
