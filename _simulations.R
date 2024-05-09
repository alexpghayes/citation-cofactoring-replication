library(targets)
library(tarchetypes)
library(crew)

tar_option_set(
  controller = crew_controller_local(workers = 8),
  packages = c(
    "fastadi", "fastRG", "furrr", "GGally", "ggraph", "glue", "here",
    "igraph", "invertiforms", "jsonlite", "kableExtra", "knitr", "LRMF3",
    "logger", "Matrix", "Rcpp", "RcppHungarian", "RSpectra", "scales", "stringr",
    "targets",
    "tidygraph", "tidyr", "tidytext", "tidyverse", "unglue", "vsp", "wordspace"
  )
)

source(here::here("R/simulation-estimators.R"))
source(here::here("R/simulation-loss.R"))
source(here::here("R/simulation-models.R"))
source(here::here("R/simulation-vsp.R"))


##### simulation options -------------------------------------------------------

models_chr <- c(
  "model4",
  "model5",
  "model6"
)

estimators_chr <- c(
  "full_svds",
  "zero_imputed_svds",
  "symmetric_svd",
  "cite_impute",
  "cite_impute2"
)

param_max_iter <- tar_target(max_iter, 15)

# TODO: set to 200 to match results from paper, but you can check that
# everything runs with just a single replication
param_num_reps <- tar_target(num_reps, 1)

param_n_seq <- tar_target(n_seq, c(100, 182, 331, 603, 1099, 2000))

param_rank <- tar_target(k, c(3, 6, 9))
param_sort_nodes <- tar_target(sort_nodes, FALSE)

##### define computations ------------------------------------------------------

### simulation helpers ---------------------------------------------------------

num_estimators <- length(estimators_chr)

model_syms <- rlang::syms(models_chr)
estimator_syms <- rlang::syms(estimators_chr)

first_n <- function(object, n) {
  UseMethod("first_n")
}

first_n.directed_dcsbm <- function(object, n) {
  object$n <- n
  object$d <- n
  object$X <- object$X[1:n, , drop = FALSE]
  object$Y <- object$Y[1:n, , drop = FALSE]
  object$z_out <- object$z_out[1:n]
  object$z_in <- object$z_in[1:n]
  object$theta_out <- object$theta_out[1:n]
  object$theta_in <- object$theta_in[1:n]
  object
}

first_n.undirected_dcsbm <- function(object, n) {
  object$n <- n
  object$X <- object$X[1:n, , drop = FALSE]
  object$z <- object$z[1:n]
  object$theta <- object$theta[1:n]
  object
}

first_n.Matrix <- function(object, n) {
  object[1:n, 1:n, drop = FALSE]
}

# return: list of subsetted distributions
construct_distribution_sequence <- function(dist, n_seq) {
  purrr::map(n_seq, \(n) first_n(dist, n))
}

target_losses <- tar_map(

  values = list(estimator = estimator_syms),

  tar_target(
    subspace_estimate,
    purrr::map(realization_seq, \(.x) estimator(.x, k = parameters$k, max_iter = max_iter)),
    pattern = map(realization_seq, parameters)
  ),

  tar_target(
    vsp_estimate,
    map(subspace_estimate, vsp, rank = parameters$k),
    pattern = map(subspace_estimate, parameters),
    iteration = "list"
  ),

  tar_target(
    alignment,
    map2(distribution_vsp, vsp_estimate, compute_alignment),
    pattern = map(distribution_vsp, vsp_estimate),
    iter = "list"
  ),

  tar_target(
    vsp_aligned,
    map2(vsp_estimate, alignment, align_factors),
    pattern = map(vsp_estimate, alignment),
    iter = "list"
  ),

  tar_target(
    loss,
    purrr::pmap_dfr(
      list(svd, subspace_estimate, vsp_aligned, distribution_vsp),
      \(a, b, c, d) loss_helper(
        svd_true = a,
        svd_estimate = b,
        true_fa = c,
        fa = d,
        params = parameters
      )
    ),
    pattern = map(svd, subspace_estimate, vsp_aligned, distribution_vsp, parameters)
  )
)

target_big <- tar_map(
  unlist = FALSE,
  values = list(model = model_syms),

  tar_target(
    max_distribution,
    model(n = max_n, k = parameters$k, sort_nodes = parameters$sort_nodes, rep = parameters$rep),
    pattern = map(parameters),
    iteration = "list"
  ),

  tar_target(
    distribution_seq,
    purrr::map(n_seq, \(n) first_n(max_distribution, n)),
    pattern = map(max_distribution),
    iteration = "list"
  ),

  tar_target(
    max_realization,
    sample_sparse(max_distribution),
    pattern = map(max_distribution),
    iteration = "list"
  ),

  tar_target(
    realization_seq,
    purrr::map(n_seq, \(n) first_n(max_realization, n)),
    pattern = map(max_realization)
  ),

  tar_target(
    svd,
    purrr::map(distribution_seq, svds),
    pattern = map(distribution_seq),
    iter = "list"
  ),

  tar_target(
    distribution_vsp,
    purrr::map(svd, \(s) vsp(as_svd_like(s), rank = length(s$d))),
    pattern = map(svd),
    iteration = "list"
  ),

  target_losses
)


target_combined <- tar_combine(
  combined_losses,
  target_big[-c(1:(6 + (length(target_losses) - 1) * num_estimators))],
  command = dplyr::bind_rows(!!!.x, .id = "estimator_model")
)

list(

  param_max_iter,
  param_num_reps,
  param_n_seq,
  param_rank,
  param_sort_nodes,

  tar_target(
    max_n,
    max(n_seq)
  ),

  tar_target(
    rep,
    1:num_reps,
  ),

  tar_target(
    parameters,
    tibble::tibble(
      k = k,
      sort_nodes = sort_nodes,
      rep = rep
    ),
    pattern = cross(k, sort_nodes, rep)
  ),

  target_big,

  target_combined,

  tar_target(summarized_loss, summarize_loss(combined_losses)),

  tar_group_by(loss_by_model, summarized_loss, model),

  tar_target(
    pre_post_clip_plot,
    make_pre_post_clip_plot(loss_by_model),
    pattern = map(loss_by_model),
    format = "file"
  ),

  tar_target(
    consistency_plot,
    make_consistency_plot(loss_by_model),
    pattern = map(loss_by_model),
    format = "file"
  )
)

