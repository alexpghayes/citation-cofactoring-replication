#' Compute sin-theta distance between subspaces
#'
#' @param u An orthogonal basis for a k-dimensional subspace of
#'   n-dimensional space.
#' @param v An orthogonal basis for a k-dimensional subspace of
#'   n-dimensional space.
#'
#' @return
#' @keywords internal
subspace_loss <- function(u, v, ell) {
  # see [1] Vu and Lei 2013 section 2.3
  # and [2] Rohe, Chatterjee, Yu 2011 Annals of Statistics page 1908

  u <- normalize.cols(u, tol = 1e-10)
  v <- normalize.cols(v, tol = 1e-10)

  s <- svd(crossprod(u, v))
  loss <- ncol(u) - sum(s$d^2)  # sign indeterminate sometimes numerically
  loss * sign(loss)
}

clip_u <- function(u, ell) {
  n <- nrow(u)
  u[1:(n - ell), ]
}

clip_v <- function(v, ell) {
  n <- nrow(v)
  v[(ell + 1):n, ]
}


loss_helper <- function(svd_true, svd_estimate, true_fa, fa, params) {

  n <- nrow(svd_estimate$u)
  k <- ncol(svd_estimate$u)
  ell <- n / 10

  tibble(
    # B_loss = sqrt(norm(true_fa$B - fa$B, type = "F")),
    u_loss = subspace_loss(svd_estimate$u, svd_true$u, ell),
    v_loss = subspace_loss(svd_estimate$v, svd_true$v),
    z_loss = sqrt(norm(as.matrix(fa$Z - true_fa$Z), type = "F")^2 / prod(dim(fa$Z))),
    y_loss = sqrt(norm(as.matrix(fa$Y - true_fa$Y), type = "F")^2 / prod(dim(fa$Y))),
    u_loss_clip = subspace_loss(clip_u(svd_estimate$u, ell), clip_u(svd_true$u, ell)),
    v_loss_clip = subspace_loss(clip_v(svd_estimate$v, ell), clip_v(svd_true$v, ell)),
    z_loss_clip = sqrt(norm(as.matrix(clip_u(fa$Z - true_fa$Z, ell)), type = "F")^2 / prod(dim(clip_u(fa$Z, ell)))),
    y_loss_clip = sqrt(norm(as.matrix(clip_v(fa$Y - true_fa$Y, ell)), type = "F")^2 / prod(dim(clip_v(fa$Y, ell)))),
    n = n
  ) |>
    bind_cols(params)
}

summarize_loss <- function(combined_losses) {
  combined_losses |>
    separate(
      col = estimator_model,
      into = c("estimator", "model"),
      sep = "_(?=model[0-9]*$)"
    ) |>
    mutate(
      estimator = str_replace(estimator, "loss_", ""),
      subspace_loss = u_loss + v_loss,
      subspace_loss_clip = u_loss_clip + v_loss_clip,
      factor_loss = z_loss + y_loss,
      factor_loss_clip = z_loss_clip + y_loss_clip
    ) |>
    pivot_longer(
      c(subspace_loss, subspace_loss_clip, factor_loss, factor_loss_clip),
      names_to = "loss_type",
      values_to = "loss"
    ) |>
    summarize(
      mean = mean(loss),
      sd = sd(loss),
      .by = c(model, estimator, n, k, loss_type)
    ) |>
    mutate(
      estimator_nice = factor(
        case_when(
          estimator == "full_svds" ~ "Fully Observed SVD",
          estimator == "symmetric_svd" ~ "Symmetrized SVD",
          estimator == "zero_imputed_svds" ~ "Zero Imputed SVD",
          estimator == "cite_impute" ~ "CitationImpute (Pre-clip)",
          estimator == "cite_impute2" ~ "CitationImpute (Post-clip)",
          estimator == "adapt_impute" ~ "AdaptiveImpute",
          TRUE ~ "Unknown Estimator"
        ),
        levels = c(
          "Zero Imputed SVD",
          "Symmetrized SVD",
          "CitationImpute (Pre-clip)",
          "CitationImpute (Post-clip)",
          "AdaptiveImpute",
          "Fully Observed SVD",
          "Unknown Estimator"
        )
      ),
      estimator_nice_short = factor(
        case_when(
          estimator == "full_svds" ~ "Fully Observed SVD",
          estimator == "symmetric_svd" ~ "Symmetrized SVD",
          estimator == "zero_imputed_svds" ~ "Zero Imputed SVD",
          estimator == "cite_impute" ~ "CitationImpute",  # "official" version is post-clip
          estimator == "adapt_impute" ~ "AdaptiveImpute",
          TRUE ~ "Unknown Estimator"
        ),
        levels = c(
          "Zero Imputed SVD",
          "Symmetrized SVD",
          "CitationImpute",
          "AdaptiveImpute",
          "Fully Observed SVD",
          "Unknown Estimator"
        )
      ),
      model_nice = case_when(
        model == "model4" ~ "Model 4",
        model == "model5" ~ "Model 5",
        model == "model6" ~ "Model 6",
        TRUE ~ "Unknown Model"
      )
    )
}


make_pre_post_clip_plot <- function(summarized_loss) {

  model <- summarized_loss$model[[1]]

  plot <- summarized_loss |>
    mutate(
      loss_type_nice = factor(
        case_when(
          loss_type == "factor_loss" ~ "All factors",
          loss_type == "factor_loss_clip" ~ "Clipped factors",
          loss_type == "subspace_loss" ~ "Full subspace",
          loss_type == "subspace_loss_clip" ~ "Clipped subspace",
          TRUE ~ "Unknown Loss"
        ),
        levels = c(
          "Clipped subspace",
          "Clipped factors",
          "Full subspace",
          "All factors",
          "Unknown Loss"
        )
      )
    ) |>
    filter(loss_type == "factor_loss_clip") |>
    ggplot() +
    aes(
      x = n,
      y = mean,
      color = estimator_nice,
      label = estimator_nice
    ) +
    geom_line() +
    geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd)) +
    scale_color_viridis_d(direction = -1) +
    labs(
      y = "Loss (log scale)",
      x = "Number of nodes (log scale)",
      color = "Estimator"
    ) +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    facet_grid(
      cols = vars(k),
      labeller = labeller(.cols = label_both),
      scales = "free_y"
    ) +
    theme_minimal(16) +
    theme(
      panel.spacing = unit(0.7, "lines")
    )


  width <- 10
  height <- 8 * 9/16

  base_path <- here("figures/simulations")

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  plot_path_pdf <- glue("{base_path}/pre_post_{model}.pdf")

  ggsave(
    plot_path_pdf,
    plot = plot,
    width = width,
    height = height,
    dpi = 600
  )

  c(plot_path_pdf)
}

make_consistency_plot <- function(summarized_loss) {
  model <- summarized_loss$model[[1]]

  plot <- summarized_loss |>
    filter(stringr::str_detect(loss_type, "clip"), estimator != "cite_impute2") |>
    mutate(
      loss_type = factor(
        case_when(
          loss_type == "factor_loss_clip" ~ "Clipped factors",
          loss_type == "subspace_loss_clip" ~ "Clipped subspace",
          TRUE ~ "Unknown Loss"
        ),
        levels = c(
          "Clipped factors",
          "Clipped subspace",
          "Unknown Loss"
        )
      )
    ) |>
    # filter(model == "model4") |>
    ggplot() +
    aes(
      x = n,
      y = mean,
      color = estimator_nice_short,
      label = estimator_nice_short
    ) +
    geom_line() +
    geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd)) +
    scale_color_viridis_d(direction = -1) +
    labs(
      y = "Loss (log scale)",
      x = "Number of nodes (log scale)",
      color = "Estimator"
    ) +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    facet_grid(
      cols = vars(k),
      rows = vars(loss_type),
      scales = "free_y",
      labeller = labeller(k = label_both)
    ) +
    theme_minimal(16) +
    theme(
      panel.spacing = unit(0.7, "lines")
    )


  width <- 10
  height <- 8 * 9/16


  base_path <- here("figures/simulations")

  if (!dir.exists(base_path))
    dir.create(base_path, recursive = TRUE)

  plot_path_pdf <- glue("{base_path}/consistency_{model}.pdf")

  ggsave(
    plot_path_pdf,
    plot = plot,
    width = width,
    height = height,
    dpi = 600
  )

  c(plot_path_pdf)
}

