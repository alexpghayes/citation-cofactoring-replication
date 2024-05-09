library(crew)
library(targets)
library(tarchetypes)

library(logger)

log_threshold(INFO) # control citation_impute() output

tar_option_set(
  # controller = crew_controller_local(workers = 8),
  packages = c(
    "fastadi", "fastRG", "furrr", "GGally", "ggraph", "glue", "here",
    "igraph", "invertiforms", "jsonlite", "kableExtra", "knitr",
    "LRMF3", "Matrix", "Rcpp", "RcppHungarian", "RSpectra", "scales", "stringr",
    "targets", "tidygraph", "tidytext", "tidyverse", "unglue", "vsp", "wordspace"
  ),
  imports = "vsp"
)

source(here::here("R/analysis-create-graph.R"))
source(here::here("R/analysis-estimators.R"))
source(here::here("R/analysis-interpret-vsp.R"))
source(here::here("R/analysis-forward-citation.R"))
source(here::here("R/analysis-keywords.R"))

estimators_chr <- c(
  "zero_imputed_svd",
  "citation_impute_approximate",
  "citation_impute_preclipped2"
)

estimators_sym <- rlang::syms(estimators_chr)

list(
  # NOTE: if you change the rank vector, you must update the index argument for the rank
  #       of interest in the forward citation table targets
  tar_target(rank, c(5, 10, 20, 30, 40)),

  tar_target(additional, 100),
  tar_target(max_iter, 50),

  tar_target(data_directory, here::here("data-raw/2024-extract/")),

  tar_target(
    hyperparameters,
    tibble::tibble(
      rank = rank,
      additional = additional,
      max_iter = max_iter
    ),
    pattern = cross(rank, additional, max_iter)
  ),

  tar_target(
    full_graph,
    create_graph(data_directory),
    deployment = "main"
  ),

  tar_target(graph, largest_connected_component(full_graph)),

  tar_target(node_data, get_node_data(graph)),
  tar_target(title_words, get_title_words(node_data)),

  tar_target(A, as_adjacency_matrix(graph, sparse = TRUE)),

  tar_target(
    ell,
    tibble::tibble(
      y = c(1, 10000, 25000, 50000, 30000, 60000, 70000, 50000),
      z = c(1, 10000, 25000, 50000, 80000, 60000, 70000, 100000)
    )
  ),

  tar_map(

    values = tibble::tibble(
      estimator = estimators_sym,
      estimator_name = estimators_chr,
    ),

    tar_target(
      subspace,
      estimator(
        A = A,
        rank = hyperparameters$rank,
        additional = hyperparameters$additional,
        max_iter = hyperparameters$max_iter
      ),
      pattern = map(hyperparameters),
      deployment = "main"
    ),

    tar_target(
      fa,
      vsp(
        subspace,
        rank = hyperparameters$rank,
        rownames = rownames(A),
        colnames = colnames(A)
      ),
      pattern = map(subspace, hyperparameters),
      iteration = "list"
    ),

    tar_target(
      vsp_clipped,
      clip_vsp(fa, ell_y = ell$y, ell_z = ell$z),
      pattern = cross(fa, ell),
      iteration = "list"
    ),

    tar_target(
      hub_csvs,
      extract_hubs(vsp_clipped, node_data, hyperparameters, ell$y, ell$z, estimator_name),
      pattern = map(vsp_clipped, cross(hyperparameters, ell)),
      format = "file"
    ),

    tar_target(
      keyword_csvs,
      create_keyword_csvs(vsp_clipped, title_words, hyperparameters, ell$y, ell$z, estimator_name, num_best = 6),
      pattern = map(vsp_clipped, cross(hyperparameters, ell)),
      format = "file",
      cue = tar_cue(file = FALSE)  # don't re-build keywords csvs when analysts manually name factors in them!
    ),

    tar_target(
      keyword_tex, # don't edit the tex by hand, only the csvs!
      create_keyword_tables(keyword_csvs, hyperparameters, ell$y, ell$z, estimator_name),
      pattern = map(keyword_csvs, cross(hyperparameters, ell)),
      format = "file",
      cue = tar_cue(mode = "always") # depends on keyword_csvs
    ),

    tar_target(
      hub_tex,
      create_hub_tables(keyword_csvs, hub_csvs, hyperparameters, ell$y, ell$z, estimator_name),
      pattern = map(keyword_csvs, hub_csvs, cross(hyperparameters, ell)),
      format = "file",
      cue = tar_cue(mode = "always") # depends on keyword_csvs
    ),

    tar_target(
      diagnostic_plots,
      extract_diagnostic_plots(vsp_clipped, hyperparameters, ell$y, ell$z, estimator_name),
      pattern = map(vsp_clipped, cross(hyperparameters, ell)),
      format = "file"
    ),

    # use this to find the correct index to slice() into below
    tar_target(
      dummy,
      investigate(vsp_clipped, hyperparameters, ell$y, ell$z, estimator_name),
      pattern = slice(map(vsp_clipped, cross(hyperparameters, ell)), index = 32),
      cue = tar_cue(mode = "always")
    ),

    tar_target(
      imputed_degrees,
      impute_degrees(vsp_clipped, node_data, estimator_name, num_threads = 6),
      deployment = "main",
      pattern = slice(vsp_clipped, index = 32)
    ),

    tar_target(
      forward_table,
      create_forward_citation_table(imputed_degrees, hyperparameters, ell$y, ell$z, estimator_name, num_top_imputations = 15),
      pattern = map(imputed_degrees, slice(cross(hyperparameters, ell), index = 32)),
      format = "file"
    ),

    tar_target(
      backward_table,
      create_backward_citation_table(imputed_degrees, hyperparameters, ell$y, ell$z, estimator_name, num_top_imputations = 15),
      pattern = map(imputed_degrees, slice(cross(hyperparameters, ell), index = 32)),
      format = "file"
    )

  ),

  tar_target(
    hyperparameters_preclipped,
    tibble(
      rank = 30,
      additional = 100,
      max_iter = 50
    )
  ),

  # if these end up looking really good, we'll probably want to compute
  # a full grid of ell x hyperparameters, which will take forever to run
  tar_target(
    subspace_preclipped,
    citation_impute_preclipped(
      A = A,
      rank = hyperparameters_preclipped$rank,
      additional = hyperparameters_preclipped$additional,
      max_iter = hyperparameters_preclipped$max_iter,
      ell_z = ell$z,
      ell_y = ell$y
    ),
    pattern = map(ell)
  ),

  tar_target(
    fa_preclipped,
    vsp(
      subspace_preclipped,
      rank = hyperparameters_preclipped$rank,
      rownames = rownames(A),
      colnames = colnames(A)
    ),
    pattern = map(subspace_preclipped),
    iteration = "list"
  ),

  tar_target(
    hub_csvs_preclipped,
    extract_hubs(fa_preclipped, node_data, hyperparameters_preclipped, ell$y, ell$z, "citation_impute_preclipped"),
    pattern = map(fa_preclipped, ell),
    format = "file"
  ),

  tar_target(
    keyword_csvs_preclipped,
    create_keyword_csvs(fa_preclipped, title_words, hyperparameters_preclipped, ell$y, ell$z, "citation_impute_preclipped", num_best = 6),
    pattern = map(fa_preclipped, ell),
    format = "file",
    cue = tar_cue(file = FALSE)  # don't re-build keywords csvs when analysts manually name factors in them!
  ),

  tar_target(
    diagnostic_plots_preclipped,
    extract_diagnostic_plots(fa_preclipped, hyperparameters_preclipped, ell$y, ell$z, "citation_impute_preclipped"),
    pattern = map(fa_preclipped, ell),
    format = "file"
  )

)

