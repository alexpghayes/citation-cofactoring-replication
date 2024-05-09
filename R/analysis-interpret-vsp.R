extract_hubs <- function(fa, node_data, hyperparameters, ell_y, ell_z, estimator) {

  base_path <- here(glue("estimates/hubs/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  y_hubs <- fa |>
    get_varimax_y() |>
    tidyr::gather(factor, loading, dplyr::contains("y"), -id) |>
    dplyr::group_by(factor) |>
    dplyr::slice_max(order_by = loading, n = 10, with_ties = FALSE) |>
    left_join(node_data, by = c("id" = "name")) |>
    select(-authors, -month) |>
    arrange(factor, desc(loading))

  z_hubs <- fa |>
    get_varimax_z() |>
    tidyr::gather(factor, loading, dplyr::contains("z"), -id) |>
    dplyr::group_by(factor) |>
    dplyr::slice_max(order_by = loading, n = 10, with_ties = FALSE) |>
    left_join(node_data, by = c("id" = "name")) |>
    select(-authors, -month) |>
    arrange(factor, desc(loading))

  y_hub_csv_path <- glue("{base_path}/y-hubs.csv")
  z_hub_csv_path <- glue("{base_path}/z-hubs.csv")

  write_csv(y_hubs, y_hub_csv_path)
  write_csv(z_hubs, z_hub_csv_path)

  c(y_hub_csv_path, z_hub_csv_path)
}

clean_title <- function(title) {
  title |>
    str_to_title() |>
    str_replace_all("<I>|<E>", " $") |>
    str_replace_all("<\\/I>|<\\/E>", "$ ") |>
    str_replace_all("<Sub>", "$_{") |>
    str_replace_all("<Sup>", "$^{") |>
    str_replace_all("</Sub>|</Sup>", "}$")
}

create_hub_tables <- function(keyword_paths, hub_paths, hyperparameters, ell_y, ell_z, estimator) {

  base_path <- here(glue("estimates/hubs/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  y_keyword_path <- keyword_paths[1]
  z_keyword_path <- keyword_paths[2]

  y_hub_path <- hub_paths[1]
  z_hub_path <- hub_paths[2]

  y_tex_path <- glue("{base_path}/y_hubs.tex")
  z_tex_path <- glue("{base_path}/z_hubs.tex")

  y_hub_path |>
    read_csv(show_col_types = FALSE) |>
    arrange(factor) |>
    slice(1:5, .by = factor) |>
    select(factor, title, in_degree, out_degree) |>
    mutate_at(vars(title), clean_title) |>
    rename(
      "ID" = "factor",
      "Title" = "title",
      "Cites" = "out_degree",
      "Cited by" = "in_degree"
    ) |>
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Y (incoming citation) factor hubs",
      label = glue("{estimator}-y-hubs"),
      escape = FALSE
    ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position", "repeat_header"),
      font_size = 8
    ) |>
    column_spec(2, width = "35em") |>
    writeLines(y_tex_path)

  z_hub_path |>
    read_csv(show_col_types = FALSE) |>
    arrange(factor) |>
    slice(1:5, .by = factor) |>
    select(factor, title, in_degree, out_degree) |>
    mutate_at(vars(title), clean_title) |>
    rename(
      "ID" = "factor",
      "Title" = "title",
      "Cites" = "out_degree",
      "Cited by" = "in_degree"
    ) |>
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Z (outgoing citation) factor hubs",
      label = glue("{estimator}-z-hubs"),
      escape = FALSE
    ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position", "repeat_header"),
      font_size = 8
    ) |>
    column_spec(2, width = "35em") |>
    writeLines(z_tex_path)

  c(y_tex_path, z_tex_path)
}

extract_keywords <- function(fa, node_data, hyperparameters) {

  key <- unglue(tar_name(), "{junk}_{estimator=adaptive_initialize|zero_imputed_svd|citation_impute_svd|citation_impute_approximate|citation_impute_preclipped}_{matrix=A|L|A_permuted}_{hash=[a-z0-9]*}")[[1]]

  base_path <- glue("estimates/keywords/{key$matrix}-rank-{hyperparameters$rank}/")

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  path <- here(glue("{base_path}/{key$estimator}-keywords.pdf"))

  rmarkdown::render(
    input = here("inst/factor-keyword-template.Rmd"),
    output_file = path,
    params = list(
      node_data = node_data,
      fa = fa,
      rank = key$rank,
      max_iter = hyperparameters$max_iter,
      matrix = key$matrix,
      estimator = key$estimator
    ),
    quiet = TRUE
  )

  path
}

plot_B_estimate <- function(X, ...) {
  as_tibble(as.matrix(X), rownames = "row") |>
    tidyr::gather(col, value, -row) |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(
      fill = "Citation\npropensity", x = "", y = ""
    ) +
    theme_minimal(20) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90)
    )
}


investigate <- function(fa, hyperparameters, ell_y, ell_z, estimator) {
  log_info(glue("estimator: {estimator}, rank: {hyperparameters$rank}, ell_y: {ell_y}, ell_z: {ell_z}"))
  1L
}

extract_diagnostic_plots <- function(fa, hyperparameters, ell_y, ell_z, estimator) {

  base_path <- here(glue("estimates/diagnostic_plots/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  k <- fa$rank

  y_plot <- plot_varimax_y_pairs(fa, factors = 1:k)
  z_plot <- plot_varimax_z_pairs(fa, factors = 1:k)

  b_plot <- plot_mixing_matrix(fa) +
    theme_minimal(18) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient2() +
    labs(
      fill = "Citation\npropensity", x = "", y = ""
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )

  y_factor_path <- glue("{base_path}/factors-y.png")
  z_factor_path <- glue("{base_path}/factors-z.png")
  ipr_path <- glue("{base_path}/ipr-pairs-plots.png")
  b_path <- glue("{base_path}/mixing-matrix.pdf")
  scrplot_path <- glue("{base_path}/screeplot.png")

  ggsave(y_factor_path, y_plot, width = min(k, 15), height = min(k, 15))
  ggsave(z_factor_path, z_plot, width = min(k, 15), height = min(k, 15))
  ggsave(b_path, b_plot, dpi = 300)

  c(
    y_factor_path,
    z_factor_path,
    b_path
  )
}
