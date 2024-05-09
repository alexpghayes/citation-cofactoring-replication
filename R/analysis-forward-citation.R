Rcpp::sourceCpp(here::here("src/estimate-degrees.cpp"))

lower_triangle_abs_rowsums <- function(s, num_threads) {
  drop(lower_triangle_abs_rowsums_impl(s$u, s$d, s$v, num_threads))
}

lower_triangle_abs_colsums <- function(s, num_threads) {
  drop(lower_triangle_abs_colsums_impl(s$u, s$d, s$v, num_threads))
}

impute_degrees <- function(mf, node_data, estimator, num_threads) {

  # NOTE: based on how we've defined A

  # rowSums is (absolute) out-degree
  # colSums is (absolute) in-degree

  # skip some extraneous computation
  if (estimator != "citation_impute_preclipped2") {
    out <- node_data |>
      mutate(
        in_imputed = NA,
        out_imputed = NA,
        total_in_degree = NA,
        total_out_degree = NA
      )

    return(out)
  }

  imputed_out_degree <- lower_triangle_abs_rowsums(mf, num_threads)
  imputed_in_degree <- lower_triangle_abs_colsums(mf, num_threads)

  node_data |>
    mutate(
      in_imputed = imputed_in_degree,
      out_imputed = imputed_out_degree,
      total_in_degree = in_degree + imputed_in_degree,
      total_out_degree = out_degree + imputed_out_degree
    )
}

create_forward_citation_table <- function(imputed_degrees, hyperparameters, ell_y, ell_z, estimator, num_top_imputations) {

  base_path <- here(glue("estimates/forward_citations/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  path <- glue("{base_path}/forward_citation_table.tex")

  imputed_degrees %>%
    arrange(desc(in_imputed)) %>%
    slice(1:num_top_imputations) %>%
    mutate(
      Title = glue("{stringr::str_to_sentence(title)} ({year})"),
      in_imputed = round(in_imputed)
    ) %>%
    rename(
      `Cited by` = in_degree,
      `Imputed` = in_imputed,
      Year = year
    ) %>%
    select(Title, `Imputed`, `Cited by`) %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Imputed incoming citations (identified edges only)",
      label = glue("{estimator}-forward-citations")
    ) %>%
    kable_styling(
      latex_options = "scale_down",
      font_size = 7.6
    ) %>%
    column_spec(1, width = "46em") %>%
    writeLines(path)

  path
}


create_backward_citation_table <- function(imputed_degrees, hyperparameters, ell_y, ell_z, estimator, num_top_imputations) {

  base_path <- here(glue("estimates/forward_citations/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  path <- glue("{base_path}/backward_citation_table.tex")

  imputed_degrees %>%
    arrange(desc(out_imputed)) %>%
    slice(1:num_top_imputations) %>%
    mutate(
      Title = glue("{stringr::str_to_sentence(title)} ({year})"),
      out_imputed = round(out_imputed)
    ) %>%
    rename(
      `Cites` = out_degree,
      `Imputed` = out_imputed,
      Year = year
    ) %>%
    select(Title, `Imputed`, `Cites`) %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Imputed outgoing citations (identified edges only)",
      label = glue("{estimator}-backward-citations")
    ) %>%
    kable_styling(
      latex_options = "scale_down",
      font_size = 7.6
    ) %>%
    column_spec(1, width = "46em") %>%
    writeLines(path)

  path
}

