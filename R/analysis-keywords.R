get_title_words <- function(node_data) {
  words <- node_data %>%
    unnest_tokens(word, title) %>%
    cast_sparse(name, word)

  appearances <- colSums(words)

  words[, appearances > 5]
}

create_keyword_csvs <- function(fa, title_words, hyperparameters, ell_y, ell_z, estimator, num_best) {

  y_keywords <- bff(fa$Y, title_words, num_best = num_best)
  z_keywords <- bff(fa$Z, title_words, num_best = num_best)

  y_factor_words <- y_keywords %>%
    mutate(
      id = colnames(fa$Y),
      words = pmap_chr(select(y_keywords, -factor), paste, sep = ", ")
    ) %>%
    select(factor, id, words)

  z_factor_words <- z_keywords %>%
    mutate(
      id = colnames(fa$Z),
      words = pmap_chr(select(z_keywords, -factor), paste, sep = ", ")
    ) %>%
    select(factor, id, words)

  base_path <- here(glue("estimates/keywords/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  y_path <- glue("{base_path}/y_keywords.csv")
  z_path <- glue("{base_path}/z_keywords.csv")

  write_csv(y_factor_words, y_path)
  write_csv(z_factor_words, z_path)

  c(y_path, z_path)
}

create_keyword_tables <- function(keyword_csv_paths, hyperparameters, ell_y, ell_z, estimator) {

  base_path <- here(glue("estimates/keywords/A-rank-{hyperparameters$rank}/{estimator}-elly-{ell_y}-ellz-{ell_z}/"))

  y_csv_path <- keyword_csv_paths[1]
  z_csv_path <- keyword_csv_paths[2]

  y_tex_path <- glue("{base_path}/y_keywords.tex")
  z_tex_path <- glue("{base_path}/z_keywords.tex")

  y_csv_path %>%
    read_csv(show_col_types = FALSE) %>%
    arrange(id) %>%
    select(factor, words, id) %>%
    rename("Factor Name" = "factor", "ID" = "id", "Top words" = "words") %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Keywords for Y (incoming citation) factors",
      label = glue("{estimator}-y-keywords")
    ) %>%
    kable_styling(
      latex_options = "scale_down",
      font_size = 8
    ) %>%
    column_spec(2, width = "34em") %>%
    writeLines(y_tex_path)

  z_csv_path %>%
    read_csv(show_col_types = FALSE) %>%
    arrange(id) %>%
    select(factor, words, id) %>%
    rename("Factor Name" = "factor", "ID" = "id", "Top words" = "words") %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      caption = "Keywords for Z (outgoing citation) factors",
      label = glue("{estimator}-z-keywords")
    ) %>%
    kable_styling(
      latex_options = "scale_down",
      font_size = 8
    ) %>%
    column_spec(2, width = "34em") %>%
    writeLines(z_tex_path)

  c(y_tex_path, z_tex_path)
}
