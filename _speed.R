library(targets)
library(tarchetypes)
library(logger)

log_threshold(WARN, namespace = "fastadi")

options(
  tidyverse.quiet = TRUE,
  Ncpus = 1
)

tar_option_set(
  packages = c(
   "fastRG", "Matrix", "softImpute", "fastadi", "bench", "here", "tidyverse"
  )
)

Sys.getenv("MKL_NUM_THREADS")
Sys.getenv("OPENBLAS_NUM_THREADS")

create_data <- function(n, k, expected_degree) {

  diagonal <- 0.8
  off_diagonal <- 0.01

  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal

  theta_in <- 1 + rexp(n, 1/8)
  theta_out <- 1 + rexp(n, 1/8)

  EA <- directed_dcsbm(
    theta_in = theta_in / sum(theta_in),
    theta_out = theta_out / sum(theta_out),
    B = B,
    expected_density = expected_degree / n,
    sort_nodes = FALSE
  )

  as(triu(sample_sparse(EA)), "dgCMatrix")
}

explicitly_observe_zeroes_sparse <- function(X, eps = 0.00001) {
  stopifnot("X must be square" = NROW(X) == NCOL(X))

  n <- NROW(X)

  indices <- expand.grid(i = 1:n, j = 1:n) %>%
    dplyr::filter(i <= j) %>%
    dplyr::left_join(summary(X), by = c("i", "j")) %>%
    dplyr::mutate(x = dplyr::if_else(is.na(x), eps, x))

  sparseMatrix(i = indices$i, j = indices$j, x = indices$x, dims = c(n, n))
}

explicitly_observe_zeroes_dense <- function(X) {
  Y <- as.matrix(X)
  Y[lower.tri(Y)] <- NA
  Y
}

##### performance stuff --------------------------------------------------------

sparse_als <- function(A, rank = 5, lambda = 0.001, iter = 5) {
  Ai <- as(explicitly_observe_zeroes_sparse(A), "Incomplete")
  softImpute(Ai, rank.max = rank, type = "als", maxit = iter, trace.it = FALSE, final.svd = FALSE, lambda = lambda)
}

sparse_svd <- function(A, rank = 5, lambda = 0.001, iter = 5) {
  Ai <- as(explicitly_observe_zeroes_sparse(A), "Incomplete")
  softImpute(Ai, rank.max = rank, type = "svd", maxit = iter, trace.it = FALSE, final.svd = FALSE, lambda = lambda)
}

dense_als <- function(A, rank = 5, lambda = 0.001, iter = 5) {
  Ad <- explicitly_observe_zeroes_dense(A)
  softImpute(Ad, rank.max = rank, type = "als", maxit = iter, trace.it = FALSE, final.svd = FALSE, lambda = lambda)
}

dense_svd <- function(A, rank = 5, lambda = 0.001, iter = 5) {
  Ad <- explicitly_observe_zeroes_dense(A)
  softImpute(Ad, rank.max = rank, type = "svd", maxit = iter, trace.it = FALSE, final.svd = FALSE, lambda = lambda)
}

citation_svd <- function(A, rank = 5, iter = 5) {
  citation_impute(A, rank = rank, check_interval = NULL, max_iter = iter)
}

run_performance_comparison <- function(n, k, iter) {

  A <- create_data(n, k, k)

  bench::mark(
    sparse_als(A, rank = k, iter = iter),
    sparse_svd(A, rank = k, iter = iter),
    dense_als(A, rank = k, iter = iter),
    dense_svd(A, rank = k, iter = iter),
    citation_svd(A, rank = k, iter = iter),
    check = FALSE,
    max_iterations = 10
  ) %>%
    mutate(n = n, k = k, iter = iter)
}

plot_memory_use <- function(performance) {

  plot <- performance %>%
    bind_rows() %>%
    mutate(
      meth = as.character(expression),
      meth2 = case_when(
        stringr::str_detect(meth, "citation_svd") ~ "CitationImpute",
        stringr::str_detect(meth, "dense_svd") ~ "softImpute\n(dense)",
        stringr::str_detect(meth, "sparse_svd") ~ "softImpute\n(sparse)",
        stringr::str_detect(meth, "dense_als") ~ "softALS\n(dense)",
        stringr::str_detect(meth, "sparse_als") ~ "softALS\n(sparse)"
      ),
      method = forcats::fct_reorder2(meth2, n, mem_alloc)
    ) %>%
    ggplot(aes(n, mem_alloc, color = method)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 1) +
    geom_point() +
    scale_color_viridis_d(direction = -1) +
    scale_y_continuous(labels = scales::label_bytes("auto_binary")) +
    theme_minimal(12) +
    labs(
      title = "Memory used",
      x = "Nodes",
      color = ""
    ) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom"
    )

  if (!dir.exists(here("figures/performance")))
    dir.create(here("figures/performance"))

  path_pdf <- here("figures/performance/memory-use.pdf")

  ggsave(
    path_pdf,
    plot = plot,
    width = 6,
    height = 6
  )

  path_png <- here("figures/performance/memory-use.png")

  ggsave(
    path_png,
    plot = plot,
    width = 6,
    height = 6
  )

  c(path_pdf, path_png)
}

plot_computation_time <- function(performance) {

  plot <- performance %>%
    bind_rows() %>%
    mutate(
      meth = as.character(expression),
      meth2 = case_when(
        stringr::str_detect(meth, "citation_svd") ~ "CitationImpute",
        stringr::str_detect(meth, "dense_svd") ~ "softImpute\n(dense)",
        stringr::str_detect(meth, "sparse_svd") ~ "softImpute\n(sparse)",
        stringr::str_detect(meth, "dense_als") ~ "softALS\n(dense)",
        stringr::str_detect(meth, "sparse_als") ~ "softALS\n(sparse)"
      ),
      method = forcats::fct_reorder2(meth2, n, mem_alloc)
    )  %>%
    ggplot(aes(n, median, color = method)) +
    # stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
    geom_line(size = 1) +
    geom_point() +
    scale_color_viridis_d(direction = -1) +
    scale_y_continuous() +
    theme_minimal(12) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(
      title = "Computation time",
      y = "Runtime (seconds)",
      x = "Nodes",
      color = ""
    ) +
    theme(
      legend.position = "bottom"
    )

  if (!dir.exists(here("figures/performance")))
    dir.create(here("figures/performance"))

  path_pdf <- here("figures/performance/computation-time.pdf")

  ggsave(
    path_pdf,
    plot = plot,
    width = 6,
    height = 6
  )

  path_png <- here("figures/performance/computation-time.png")

  ggsave(
    path_png,
    plot = plot,
    width = 6,
    height = 6
  )

  c(path_pdf, path_png)

}

list(

  tar_target(n, c(50, 500, 1000, 2500, 5000, 10000)),

  tar_target(k, 5),

  tar_target(iter, 1),

  tar_target(
    grid,
    tibble::tibble(
      n = n,
      k = k,
      iter = iter
    ),
    pattern = cross(n, k, iter)
  ),

  tar_target(
    performance,
    run_performance_comparison(grid$n, grid$k, grid$iter),
    pattern = map(grid),
    iteration = "list"
  ),

  tar_target(
    memory_plot,
    plot_memory_use(performance),
    format = "file"
  ),

  tar_target(
    time_plot,
    plot_computation_time(performance),
    format = "file"
  )
)







