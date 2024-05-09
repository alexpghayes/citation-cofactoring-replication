read_from_path <- function(path) {

  # some files are empty and cause problems here,
  # in particular files that are 55 bytes in size it seems, which
  # seems to be some ESCI files

  tbl <- path |>
    gzfile() |>
    stream_in(verbose = FALSE) |>
    as_tibble()

  if (NROW(tbl) < 1) {
    return(tibble::tibble())
  }

  tbl |>
    select(
      id,
      title,
      source_title,  # journal
      names,         # authors
      pub_year,
      pub_month,
      references,
      abstract_text
    ) |>
    mutate(
      authors = map(names, "display_name"),
      cites = map(references, "id"),
      # some months are coded as APR-MAY pairs,
      # take the first month in these cases
      pub_month = substring(pub_month, 1, 3)
    ) |>
    select(
      -names,
      -references
    )
}

create_graph <- function(data_path) {

  plan(multisession, workers = 20)
  on.exit(plan(sequential))

  raw_data_paths <- list.files(
    data_path,
    pattern = ".json.gz",
    full.names = TRUE
  )

  data <- future_map_dfr(raw_data_paths, read_from_path)

  node_data <- data |>
    mutate(
      year = as.numeric(pub_year),
      month = match(tolower(pub_month), tolower(month.abb))
    ) |>
    select(-cites, -pub_year, -pub_month)

  edge_data <- data |>
    rename(from = id, to = cites) |>
    select(from, to) |>
    unnest(to) |>
    na.omit()

  # ~262k papers, ~2M edges, measured by tidygraph

  graph <- edge_data |>
    as_tbl_graph() |>
    activate(nodes) |>
    left_join(node_data, by = c("name" = "id")) |>
    filter(
      name %in% node_data$id
    )

  graph
}

largest_connected_component <- function(graph) {

  # subset to largest connected component
  comps <- components(graph, mode = "weak")
  largest_component_index <- which.max(comps$csize)

  weak_component <- graph |>
    filter(comps$membership == largest_component_index)

  # 197k nodes, 1.5 M edges, measured by tidygraph (in largest component)

  # time_centile of 1 is the earliest (i.e. papers from 2020)

  chronological_component <- weak_component |>
    arrange(desc(year), desc(month)) |>
    mutate(
      index = row_number(),
      time_centile = cut(index, breaks = 100, labels = FALSE)
    )

  chronological_component
}

get_node_data <- function(graph) {
  graph |>
    activate(nodes) |>
    mutate(
      in_degree = centrality_degree(mode = "in"),
      out_degree = centrality_degree(mode = "out")
    ) |>
    as_tibble()
}
