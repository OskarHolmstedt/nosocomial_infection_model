# MCMC convergence diagnostics: five-panel figure.
#
# Panel 1 (top-left)    — trace of n_imports with running mean; stationarity.
# Panel 2 (top-right)   — ACF of n_imports with ESS estimate; mixing speed.
# Panel 3 (mid-left)    — trace + running mean of mean infection time; second
#                          scalar to confirm stationarity is not artefactual.
# Panel 4 (mid-right)   — ancestry raster for the n_raster most uncertain cases;
#                          rapid colour switching = good mixing, flatlines = stuck.
# Panel 5 (bottom, full width) — per-case posterior entropy bar chart, sorted
#                          high to low. If outbreak_data is provided, bars are
#                          coloured green (mode matches truth) or red (wrong).
#
# outbreak_data: optional; adds truth reference lines and mode-accuracy colours.
# thin:     plot every `thin`-th iteration in the raster.
# n_raster: number of cases shown in the raster panel.
# acf_lags: maximum lag for the ACF panel.
plot_mcmc_diagnostics <- function(Y, outbreak_data = NULL,
                                   thin = 10, n_raster = 8, acf_lags = 100) {
  library(patchwork)

  N       <- nrow(Y$anc)
  n_cases <- ncol(Y$anc)
  iters   <- seq_len(N)

  # ── Per-case posterior entropy (shared by raster + bar chart) ──────────────
  entropy <- apply(Y$anc, 2, function(x) {
    p <- prop.table(table(factor(x), useNA = "always"))
    -sum(p * log(p + 1e-10))
  })

  # ── Scalar 1: n_imports ────────────────────────────────────────────────────
  n_imp  <- rowSums(is.na(Y$anc))
  df_imp <- data.frame(iter = iters, n_imp = n_imp,
                        running = cumsum(n_imp) / iters)

  true_imp <- if (!is.null(outbreak_data))
    sum(is.na(outbreak_data$ObsRec$Anc2)) else NA_real_

  p_trace <- ggplot(df_imp, aes(x = iter)) +
    geom_line(aes(y = n_imp), linewidth = 0.25, color = "steelblue", alpha = 0.7) +
    geom_line(aes(y = running), linewidth = 0.9, color = "steelblue") +
    { if (!is.na(true_imp))
        geom_hline(yintercept = true_imp, linetype = "dashed",
                   color = "#e07b39", linewidth = 0.7) } +
    labs(x = "Iteration", y = "N community imports",
         title = "Trace + running mean: community imports",
         subtitle = if (!is.na(true_imp))
           sprintf("orange dashed = truth (%g)", true_imp) else NULL)

  # ── Scalar 2: mean infection time ──────────────────────────────────────────
  mean_inf <- rowMeans(Y$inf)
  df_inf   <- data.frame(iter = iters, mean_inf = mean_inf,
                          running = cumsum(mean_inf) / iters)

  true_mean_inf <- if (!is.null(outbreak_data))
    mean(outbreak_data$ObsRec$Infc) else NA_real_

  p_inf <- ggplot(df_inf, aes(x = iter)) +
    geom_line(aes(y = mean_inf), linewidth = 0.25, color = "steelblue", alpha = 0.7) +
    geom_line(aes(y = running), linewidth = 0.9, color = "steelblue") +
    { if (!is.na(true_mean_inf))
        geom_hline(yintercept = true_mean_inf, linetype = "dashed",
                   color = "#e07b39", linewidth = 0.7) } +
    labs(x = "Iteration", y = "Mean infection time",
         title = "Trace + running mean: infection times")

  # ── ACF of n_imports ───────────────────────────────────────────────────────
  acf_obj <- acf(n_imp, lag.max = acf_lags, plot = FALSE)
  ci      <- qnorm(0.975) / sqrt(N)
  df_acf  <- data.frame(lag = as.numeric(acf_obj$lag[-1]),
                         acf = as.numeric(acf_obj$acf[-1]))

  pos_acf <- acf_obj$acf[acf_obj$acf > 0]
  ess     <- round(N / (1 + 2 * sum(pos_acf)))

  p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_col(fill = "steelblue", width = 0.8) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
               color = "#e07b39", linewidth = 0.6) +
    labs(x = "Lag", y = "ACF",
         title = "Autocorrelation: community imports",
         subtitle = sprintf("Approximate ESS = %d / %d", ess, N))

  # ── Ancestry raster ────────────────────────────────────────────────────────
  top_idx  <- order(entropy, decreasing = TRUE)[seq_len(min(n_raster, n_cases))]
  thin_idx <- seq(1, N, by = thin)
  anc_sub  <- Y$anc[thin_idx, top_idx, drop = FALSE]
  colnames(anc_sub) <- paste0("case ", top_idx)

  df_raster <- as.data.frame(anc_sub) |>
    dplyr::mutate(iter = thin_idx) |>
    tidyr::pivot_longer(-iter, names_to = "case", values_to = "ancestor") |>
    dplyr::mutate(
      case     = factor(case, levels = paste0("case ", top_idx)),
      ancestor = factor(ancestor)
    )

  p_raster <- ggplot(df_raster, aes(x = iter, y = case, fill = ancestor)) +
    geom_tile() +
    scale_fill_viridis_d(na.value = "grey88", guide = "none") +
    labs(x = "Iteration", y = "Case (high \u2192 low entropy)",
         title = sprintf("Ancestry traces: %d most uncertain cases", min(n_raster, n_cases)),
         subtitle = "grey = community import; colour = hospital ancestor")

  # ── Per-case entropy bar chart ─────────────────────────────────────────────
  # Compute posterior mode per case for correctness colouring
  mode_anc <- apply(Y$anc, 2, function(x) {
    tab <- sort(table(x, useNA = "always"), decreasing = TRUE)
    val <- names(tab)[1L]
    if (is.na(val) || val == "NA") NA_integer_ else as.integer(val)
  })

  df_ent <- data.frame(
    case    = seq_len(n_cases),
    entropy = entropy
  )

  if (!is.null(outbreak_data)) {
    true_anc     <- outbreak_data$ObsRec$Anc2
    df_ent$correct <- dplyr::case_when(
      is.na(mode_anc) & is.na(true_anc)   ~ "correct",
      is.na(mode_anc) | is.na(true_anc)   ~ "wrong",
      mode_anc == true_anc                 ~ "correct",
      TRUE                                 ~ "wrong"
    )
    fill_scale <- scale_fill_manual(
      values = c(correct = "#27ae60", wrong = "#e74c3c"),
      name = "Mode", guide = "none"
    )
  } else {
    df_ent$correct <- "unknown"
    fill_scale <- scale_fill_manual(values = c(unknown = "steelblue"), guide = "none")
  }

  # Sort cases by entropy descending
  df_ent <- df_ent[order(df_ent$entropy, decreasing = TRUE), ]
  df_ent$case <- factor(df_ent$case, levels = df_ent$case)

  p_entropy <- ggplot(df_ent, aes(x = case, y = entropy, fill = correct)) +
    geom_col(width = 0.8) +
    fill_scale +
    labs(x = "Case (sorted by entropy)", y = "Posterior entropy",
         title = "Per-case posterior entropy",
         subtitle = if (!is.null(outbreak_data))
           "green = mode matches truth; red = wrong" else NULL) +
    theme(axis.text.x = element_text(size = 7))

  # ── Combine: two rows of pairs + entropy spanning full width ───────────────
  (p_trace | p_acf) / (p_inf | p_raster) / p_entropy +
    plot_layout(heights = c(1, 1, 0.6))
}

# Per-case mode vs. truth scatter: visualises the quality of wrong predictions.
#
# Each case is a point. Correct cases sit on the diagonal (mode = truth).
# Wrong cases sit above it — horizontal distance to the diagonal shows
# how close the truth came to winning; vertical distance shows how confident
# the model was in the wrong answer.
# Point size encodes difficulty: larger = fewer eligible ancestors = harder.
plot_mode_vs_truth <- function(Y, outbreak_data) {
  sc <- ancestry_score(outbreak_data$ObsRec$Anc2, Y$anc,
                       outbreak_data$ObsRec$Adm, outbreak_data$ObsRec$PTest)

  p_mode_vals <- apply(Y$anc, 2, function(x) {
    tab <- sort(table(x, useNA = "always"), decreasing = TRUE)
    as.numeric(tab[1]) / sum(tab)
  })

  true_anc   <- outbreak_data$ObsRec$Anc2
  mode_anc   <- posterior_mode_anc(Y$anc)
  correct    <- ifelse(is.na(true_anc),
                       is.na(mode_anc),
                       !is.na(mode_anc) & mode_anc == true_anc)

  df <- data.frame(
    case     = seq_len(length(sc$p_true)),
    p_true   = sc$p_true,
    p_mode   = p_mode_vals,
    n_elig   = sc$n_elig,
    correct  = correct
  )

  ggplot(df, aes(x = p_true, y = p_mode, color = correct, size = 1 / n_elig)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(xend = p_true, yend = p_true),
                 color = "grey75", linewidth = 0.4) +
    geom_point(alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = case), size = 3,
                              show.legend = FALSE, max.overlaps = 20) +
    scale_color_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                       labels = c("TRUE" = "correct", "FALSE" = "wrong"),
                       name = "Mode") +
    scale_size_continuous(range = c(2, 6),
                          name = "1 / n_eligible\n(harder \u2192 larger)") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "p_true  (posterior mass on true ancestor)",
         y = "p_mode  (posterior mass on modal ancestor)",
         title = "Mode probability vs. truth probability, per case",
         subtitle = paste("Diagonal = correct; above = wrong mode outcompeted truth;",
                          "segment shows gap to diagonal"))
}

# Plot all posterior ancestry edges with support >= min_prob.
# Width and opacity both scale with posterior probability.
# Green edges = ancestor matches truth; red = wrong ancestor.
# Community (NA) posterior samples are not drawn (no visible source node).
plot_timeline_posterior <- function(outbreak_data, Y,
                                    min_prob  = 0.05,
                                    gen_scale = 2.5,
                                    ...) {
  obs      <- outbreak_data$ObsRec
  n        <- nrow(obs)
  N        <- nrow(Y$anc)
  true_anc <- obs$Anc2

  # Base layout from the true observed tree; strip transmission edges to replace them
  out        <- plot_timeline(outbreak_data, observed = TRUE, show_tree = TRUE,
                              gen_scale = gen_scale, ...)
  is_tx      <- !grepl("^\\.", out$x$edges$from) & !grepl("^\\.", out$x$edges$to)
  out$x$edges <- out$x$edges[!is_tx, , drop = FALSE]

  # Build posterior edge table
  edge_rows <- lapply(seq_len(n), function(i) {
    real <- Y$anc[, i][!is.na(Y$anc[, i])]
    if (length(real) == 0L) return(NULL)
    tab   <- table(real)
    probs <- as.numeric(tab) / N
    ancs  <- as.integer(names(tab))
    keep  <- probs >= min_prob
    if (!any(keep)) return(NULL)
    data.frame(
      from    = obs$Id[ancs[keep]],
      to      = obs$Id[i],
      p       = probs[keep],
      correct = !is.na(true_anc[i]) & ancs[keep] == true_anc[i],
      stringsAsFactors = FALSE
    )
  })
  post_edges <- do.call(rbind, edge_rows)

  if (!is.null(post_edges) && nrow(post_edges) > 0) {
    alpha <- 0.2 + 0.8 * post_edges$p   # keep even low-prob edges somewhat visible
    post_edges$color     <- ifelse(post_edges$correct,
                                   mapply(hex_to_rgba, "#27ae60", alpha),  # green = correct
                                   mapply(hex_to_rgba, "#e74c3c", alpha))  # red   = wrong
    post_edges$width     <- gen_scale * (0.3 + 2.7 * post_edges$p)
    post_edges$arrows    <- NA_character_
    post_edges$smooth    <- FALSE
    post_edges$dashes    <- FALSE
    post_edges$arrows.to <- TRUE
    post_edges$p         <- NULL
    post_edges$correct   <- NULL
    out$x$edges <- dplyr::bind_rows(out$x$edges, post_edges)
  }
  out
}

# Overlay the posterior mode transmission tree (dashed, curved) on top of the
# true observed tree (solid, colored). Curved mode edges stay visually separate
# even when they agree with the true edge, and visibly diverge when they don't.
# mode_color / mode_alpha control the dashed overlay appearance.
plot_timeline_mcmc <- function(outbreak_data, Y,
                               mode_color  = "#222222",
                               mode_alpha  = 0.45,
                               edge_colors = c(room = "#e74c3c", ward = "#f39c12", hospital = "#3498db"),
                               gen_scale   = 2.5,
                               ...) {
  mode_anc <- posterior_mode_anc(Y$anc)
  obs      <- outbreak_data$ObsRec

  # Base: true observed tree (solid colored edges)
  out <- plot_timeline(outbreak_data, observed = TRUE, show_tree = TRUE,
                       edge_colors = edge_colors, gen_scale = gen_scale, ...)

  # Build mode edges for all non-community inferred ancestors
  has_mode <- !is.na(mode_anc)
  if (any(has_mode)) {
    idx      <- which(has_mode)
    mode_edges <- data.frame(
      from      = obs$Id[mode_anc[idx]],
      to        = obs$Id[idx],
      arrows    = NA_character_,
      color     = hex_to_rgba(mode_color, alpha = mode_alpha),
      width     = gen_scale,
      smooth    = TRUE,   # curved so they separate visually from straight true edges
      dashes    = TRUE,
      arrows.to = TRUE,
      stringsAsFactors = FALSE
    )
    out$x$edges <- dplyr::bind_rows(out$x$edges, mode_edges)
  }
  out
}

# Overlay the posterior mode transmission tree on the true observed tree.
# True tree edges are shown faded (spatial colors, low alpha) as reference.
# Mode tree edges are colored by agreement: green = correct, red = wrong.
# Cases correctly inferred as community-imported produce no edge (neither tree has one).
plot_timeline_mode <- function(outbreak_data, Y,
                               truth_alpha  = 0.20,
                               gen_scale    = 2.5,
                               edge_colors  = c(room = "#e74c3c", ward = "#f39c12", hospital = "#3498db"),
                               match_color  = "#27ae60",   # green — mode matches truth
                               miss_color   = "#e74c3c",   # red   — mode disagrees
                               mode_color    = "#2c3e50",  # neutral edge color when color_edges = FALSE
                               color_edges   = TRUE,       # if FALSE, draw mode edges in neutral color instead of green/red
                               color_nodes   = TRUE,       # if FALSE, skip green/red node border colouring
                               inf_times = "truth",        # infection times for node positions:
                                                           #   "truth" — ground truth Infc
                                                           #   "mode"  — posterior mode of Y$inf
                                                           #   "mean"  — posterior mean of Y$inf
                                                           #   named numeric vector — custom, keyed by patient Id
                               ...) {
  obs      <- outbreak_data$ObsRec
  true_anc <- obs$Anc2
  mode_anc <- posterior_mode_anc(Y$anc)

  # Resolve infection times for node positions
  inf_override <- if (is.character(inf_times)) {
    switch(inf_times,
      truth = NULL,
      mode  = setNames(apply(Y$inf, 2, function(x) as.numeric(names(which.max(table(x))))), obs$Id),
      mean  = setNames(colMeans(Y$inf), obs$Id),
      stop('inf_times must be "truth", "mode", "mean", or a named numeric vector')
    )
  } else if (is.numeric(inf_times)) {
    inf_times  # passed through directly to inf_override
  } else {
    stop('inf_times must be "truth", "mode", "mean", or a named numeric vector')
  }

  # Base: full layout + true tree in spatial colors
  out   <- plot_timeline(outbreak_data, observed = TRUE, show_tree = TRUE,
                         edge_colors = edge_colors, gen_scale = gen_scale,
                         inf_override = inf_override, ...)

  # ── Per-case posterior probabilities ────────────────────────────────────────
  # p_true[i]: fraction of samples where the true ancestor was sampled
  # p_mode[i]: fraction of samples on the modal (most frequent) ancestor
  N_samp <- nrow(Y$anc)
  p_true_vec <- vapply(seq_len(nrow(obs)), function(i) {
    ta <- true_anc[i]
    if (is.na(ta)) return(NA_real_)
    mean(!is.na(Y$anc[, i]) & Y$anc[, i] == ta)
  }, numeric(1))
  p_mode_vec <- vapply(seq_len(nrow(obs)), function(i)
    max(table(Y$anc[, i], useNA = "always")) / N_samp, numeric(1))

  # ── Identify true transmission edges ────────────────────────────────────────
  is_tx <- !grepl("^\\.", out$x$edges$from) & !grepl("^\\.", out$x$edges$to)

  # Fade true tree; width = gen_scale scaled by posterior support for truth
  true_tx       <- out$x$edges[is_tx, , drop = FALSE]
  true_tx$color <- hex_to_rgba("#999999", alpha = truth_alpha)
  to_idx        <- match(true_tx$to, obs$Id)
  true_tx$width <- pmax(gen_scale * 4 * p_true_vec[to_idx], 0.5)
  out$x$edges[is_tx, ] <- true_tx

  # ── Mode edges: width scaled by posterior support for modal ancestor ─────────
  has_mode <- !is.na(mode_anc)
  if (any(has_mode)) {
    idx      <- which(has_mode)
    correct  <- !is.na(true_anc[idx]) & mode_anc[idx] == true_anc[idx]
    edge_col <- if (color_edges) ifelse(correct, match_color, miss_color) else mode_color
    mode_edges <- data.frame(
      from      = obs$Id[mode_anc[idx]],
      to        = obs$Id[idx],
      color     = edge_col,
      width     = pmax(gen_scale * 4 * p_mode_vec[idx], 0.5),
      arrows    = NA_character_,
      smooth    = FALSE,
      dashes    = FALSE,
      arrows.to = TRUE,
      stringsAsFactors = FALSE
    )
    out$x$edges <- dplyr::bind_rows(out$x$edges, mode_edges)
  }

  # Node borders: green = mode correctly classified community vs hospital, red = misclassified.
  # This makes community-inference errors visible even when no mode edge is drawn.
  if (color_nodes) {
    true_community <- is.na(true_anc)
    mode_community <- is.na(mode_anc)
    border_col     <- ifelse(true_community == mode_community, match_color, miss_color)
    node_idx       <- match(obs$Id, out$x$nodes$id)
    valid          <- !is.na(node_idx)
    out$x$nodes$color.border[node_idx[valid]]           <- border_col[valid]
    out$x$nodes$color.highlight.border[node_idx[valid]] <- border_col[valid]
    out$x$nodes$borderWidth[node_idx[valid]]            <- 4
  }

  out
}

hex_to_rgba <- function(hex, alpha = 0.3) {
  r <- strtoi(substr(hex, 2L, 3L), 16L)
  g <- strtoi(substr(hex, 4L, 5L), 16L)
  b <- strtoi(substr(hex, 6L, 7L), 16L)
  sprintf("rgba(%d,%d,%d,%.2f)", r, g, b, alpha)
}

# Blend hex color toward white — avoids transparency stacking on overlapping elements
lighten_hex <- function(hex, amount = 0.65) {
  r <- strtoi(substr(hex, 2L, 3L), 16L)
  g <- strtoi(substr(hex, 4L, 5L), 16L)
  b <- strtoi(substr(hex, 6L, 7L), 16L)
  sprintf("#%02X%02X%02X",
          round(r + (255L - r) * amount),
          round(g + (255L - g) * amount),
          round(b + (255L - b) * amount))
}

# Gantt-style timeline: each case is a rectangle spanning admission to discharge,
# with a vertical mark at infection time and arrows for transmission edges.
plot_gantt <- function(epi, positions, node_color = "Room", num_days = NULL,
                       tick_by = 5, node_height = 0.6) {
  ll <- epi$linelist
  el <- epi$contacts

  # Bed rows only (rowname != Room column)
  bed_pos <- positions[rownames(positions) != positions$Room, ]
  bed_idx <- setNames(seq_len(nrow(bed_pos)), rownames(bed_pos))

  ll$y    <- bed_idx[ll$Bed]
  ll$ymin <- ll$y - node_height / 2
  ll$ymax <- ll$y + node_height / 2

  # Cap Dis at num_days for still-admitted patients
  if (is.null(num_days)) num_days <- max(ll$Dis, na.rm = TRUE)
  ll$Dis <- ifelse(is.na(ll$Dis), num_days, ll$Dis)

  # Edge coordinates: source Infc → dest Infc, at respective bed y-positions
  coords <- ll[, c("id", "Infc", "y")]
  el_data <- el |>
    dplyr::filter(!is.na(from)) |>
    dplyr::left_join(coords, by = c("from" = "id")) |>
    dplyr::rename(x_from = Infc, y_from = y) |>
    dplyr::left_join(coords, by = c("to" = "id")) |>
    dplyr::rename(x_to = Infc, y_to = y)

  x_breaks <- seq(0, ceiling(num_days / tick_by) * tick_by, by = tick_by)

  ggplot(ll) +
    # Stay rectangles
    geom_rect(aes(xmin = Adm, xmax = Dis,
                  ymin = ymin, ymax = ymax,
                  fill = .data[[node_color]]),
              color = "white", linewidth = 0.3, alpha = 0.85) +
    # Infection time marker
    geom_segment(aes(x = Infc, xend = Infc, y = ymin, yend = ymax),
                 color = "black", linewidth = 0.5) +
    # Transmission edges
    geom_segment(data = el_data,
                 aes(x = x_from, xend = x_to, y = y_from, yend = y_to),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 color = "gray35", linewidth = 0.4) +
    scale_x_continuous(breaks = x_breaks, expand = expansion(add = 1)) +
    scale_y_continuous(breaks = seq_len(nrow(bed_pos)),
                       labels = rownames(bed_pos),
                       expand = expansion(add = 0.5)) +
    labs(x = "Day", y = "Bed", fill = node_color) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#eeeeee"))
}

# Timeline plot: pins x by infection time, y by bed index (clusters by room/ward).
# Adds x- and y-axis lines with labels, ghost stay bars, and guide lines per bed.
# Takes OutbreakData directly (theta and positions extracted internally).
# observed = TRUE uses ObsTree/ObsRec, FALSE uses FullTree/CaseRec.
# flip = FALSE: time on x, beds on y.  flip = TRUE: beds on x, time on y (top = day 0).
plot_timeline <- function(outbreak_data, observed = TRUE, flip = FALSE, show_tree = TRUE,
                          show_all_stays = FALSE,
                          node_color = "Room", edge_label = "V3",
                          x_scale = 40, y_scale = 40, tick_by = 5, height = "600px",
                          node_font_size = 1, edge_font_size = 1,
                          ghost_width = 20, bar_gap = 8,
                          edge_colors = c(room = "#e74c3c", ward = "#f39c12", hospital = "#3498db"),
                          uniform_edge_color = NULL,  # if set, all transmission edges use this color
                          edge_width_by = "auto",     # "auto" (weight-based), "gen" (generation), "fixed" (uniform gen_scale)
                          gen_scale = 2.5,
                          inf_override = NULL) {

  # --- Unpack OutbreakData ---
  positions <- outbreak_data$theta$Positions
  tree      <- if (observed) outbreak_data$ObsTree else outbreak_data$FullTree
  linelist  <- if (observed) outbreak_data$ObsRec  else outbreak_data$CaseRec
  EL        <- cbind(as_edgelist(tree), E(tree)$weight)
  EL        <- EL[!is.na(EL[, 1]), , drop = FALSE]
  epi       <- make_epicontacts(linelist, EL, directed = TRUE)

  out   <- epicontacts:::vis_epicontacts(epi, node_color = node_color, edge_label = edge_label)
  ll    <- epi$linelist
  nodes <- out$x$nodes

  # vis_epicontacts drops linelist cases that have no edges (e.g. community-imported cases).
  # Add them back, inheriting colors from the node_color grouping of existing nodes.
  missing_ids <- setdiff(ll$id, nodes$id)
  if (length(missing_ids) > 0) {
    color_map    <- setNames(nodes$color.background, nodes[[node_color]])
    missing_ll   <- ll[ll$id %in% missing_ids, , drop = FALSE]
    fill_color   <- color_map[missing_ll[[node_color]]]
    fill_color[is.na(fill_color)] <- "#cccccc"
    missing_ll$label                      <- missing_ll$id
    missing_ll$title                      <- NA_character_
    missing_ll$color.background           <- fill_color
    missing_ll$color.highlight.background <- fill_color
    missing_ll$color.border               <- "black"
    missing_ll$color.highlight.border     <- "black"
    missing_ll$size                       <- 20
    missing_ll$borderWidth                <- 2
    nodes <- dplyr::bind_rows(nodes, missing_ll)
  }

  # Bed rows only (rowname != Room column)
  bed_pos <- positions[rownames(positions) != positions$Room, ]
  n_beds  <- nrow(bed_pos)
  bed_idx <- setNames(seq_len(n_beds), rownames(bed_pos))
  bed_of_node <- setNames(ll$Bed, ll$id)
  infc_of_node <- if (!is.null(inf_override)) {
    # Use provided override, falling back to ground truth for any missing ids
    base <- setNames(ll$Infc, ll$id)
    base[names(inf_override)] <- inf_override
    base
  } else {
    setNames(ll$Infc, ll$id)
  }

  num_days  <- max(ll$Dis, na.rm = TRUE)
  tick_days <- seq(0, ceiling(num_days / tick_by) * tick_by, by = tick_by)
  node_hex  <- setNames(nodes$color.background, nodes$id)

  # Helper: invisible anchor dot
  dot <- function(id, x, y) {
    if (length(id) == 0L || length(x) == 0L || length(y) == 0L) {
      id <- character()
      x <- y <- numeric()
    }
    n <- max(length(id), length(x), length(y))
    data.frame(
      id = rep_len(id, n),
      label = rep("", n),
      x = rep_len(x, n),
      y = rep_len(y, n),
      shape = rep("dot", n),
      size = rep(1, n),
      color.background = rep("#ffffff", n),
      color.border = rep("#ffffff", n),
      stringsAsFactors = FALSE
    )
  }

  if (!flip) {
    # --- Normal: time → x, bed index → y ---
    bed_coords  <- bed_idx[bed_of_node[nodes$id]] * y_scale
    nodes$x     <- infc_of_node[nodes$id] * x_scale
    nodes$y     <- bed_coords

    node_y   <- setNames(nodes$y, nodes$id)
    dis_x    <- (ifelse(is.na(ll$Dis), num_days, ll$Dis) + 1) * x_scale - bar_gap

    # All-patient background bars (horizontal, gray)
    if (show_all_stays) {
      fr      <- outbreak_data$FullRec
      fr      <- fr[as.character(fr$Bed) %in% names(bed_idx), ]
      fr_y    <- bed_idx[as.character(fr$Bed)] * y_scale
      fr_disx <- (pmin(fr$Dis, num_days) + 1) * x_scale - bar_gap
      bg_left  <- dot(paste0(".bg.L.", fr$Id), fr$Adm * x_scale, fr_y)
      bg_right <- dot(paste0(".bg.R.", fr$Id), fr_disx,           fr_y)
      bg_edges <- data.frame(from = paste0(".bg.L.", fr$Id), to = paste0(".bg.R.", fr$Id),
                             arrows = "", color = "#eeeeee", width = ghost_width, smooth = FALSE)
    }

    # Ghost stay bars (horizontal)
    ghost_left  <- dot(paste0(".ghost.L.", ll$id), ll$Adm * x_scale, node_y[ll$id])
    ghost_right <- dot(paste0(".ghost.R.", ll$id), dis_x,            node_y[ll$id])

    # PTest markers: short vertical edges crossing each ghost bar
    ll_pt       <- ll[!is.na(ll$PTest), ]
    ptest_x     <- ll_pt$PTest * x_scale
    ptest_y     <- node_y[ll_pt$id]
    marker_top  <- dot(paste0(".pt.T.", ll_pt$id), ptest_x, ptest_y - ghost_width / 2)
    marker_bot  <- dot(paste0(".pt.B.", ll_pt$id), ptest_x, ptest_y + ghost_width / 2)
    marker_edges <- data.frame(
      from = paste0(".pt.T.", ll_pt$id), to = paste0(".pt.B.", ll_pt$id),
      arrows = "", color = "#222222", width = 1.5, smooth = FALSE)

    # Day tick labels below last bed
    x_right    <- max(tick_days) * x_scale
    y_tick     <- n_beds * y_scale + 30L
    tick_nodes <- data.frame(
      id = paste0(".tick.", tick_days), label = as.character(tick_days),
      x = tick_days * x_scale, y = y_tick,
      shape = "text", font.size = 18, font.color = "#444444",
      stringsAsFactors = FALSE)

    # Bed labels on left
    label_nodes <- data.frame(
      id = paste0(".ylabel.", rownames(bed_pos)), label = rownames(bed_pos),
      x = -50, y = seq_len(n_beds) * y_scale,
      shape = "text", font.size = 14, font.color = "#444444",
      stringsAsFactors = FALSE)

    # Guide line anchors (horizontal, one per bed row)
    guide_L <- dot(paste0(".row.L.", seq_len(n_beds)), 0,       seq_len(n_beds) * y_scale)
    guide_R <- dot(paste0(".row.R.", seq_len(n_beds)), x_right, seq_len(n_beds) * y_scale)

    # Axis spine anchors
    spine <- dot(c(".spine.a", ".spine.b"), 0, c(y_scale, y_tick))

    bg_nodes <- if (show_all_stays) dplyr::bind_rows(bg_left, bg_right) else NULL
    out$x$nodes <- dplyr::bind_rows(nodes, tick_nodes, label_nodes,
                                    guide_L, guide_R, bg_nodes,
                                    ghost_left, ghost_right, spine,
                                    marker_top, marker_bot)

    axis_edge  <- data.frame(from = paste0(".tick.", tick_days[1]),
                             to   = paste0(".tick.", tick_days[length(tick_days)]),
                             arrows = "", color = "#aaaaaa", width = 1, smooth = FALSE)
    spine_edge <- data.frame(from = ".spine.a", to = ".spine.b",
                             arrows = "", color = "#aaaaaa", width = 1, smooth = FALSE)
    guide_edges <- data.frame(from = paste0(".row.L.", seq_len(n_beds)),
                              to   = paste0(".row.R.", seq_len(n_beds)),
                              arrows = "", color = "#dddddd", width = 0.5, smooth = FALSE)

  } else {
    # --- Flipped: bed index → x, time → y (y=0 at top = day 0) ---
    bed_coords  <- bed_idx[bed_of_node[nodes$id]] * x_scale
    nodes$x     <- bed_coords
    nodes$y     <- infc_of_node[nodes$id] * y_scale

    node_x   <- setNames(nodes$x, nodes$id)
    dis_y    <- (ifelse(is.na(ll$Dis), num_days, ll$Dis) + 1) * y_scale - bar_gap

    # All-patient background bars (vertical, gray)
    if (show_all_stays) {
      fr      <- outbreak_data$FullRec
      fr      <- fr[as.character(fr$Bed) %in% names(bed_idx), ]
      fr_x    <- bed_idx[as.character(fr$Bed)] * x_scale
      fr_disy <- (pmin(fr$Dis, num_days) + 1) * y_scale - bar_gap
      bg_left  <- dot(paste0(".bg.L.", fr$Id), fr_x, fr$Adm * y_scale)
      bg_right <- dot(paste0(".bg.R.", fr$Id), fr_x, fr_disy)
      bg_edges <- data.frame(from = paste0(".bg.L.", fr$Id), to = paste0(".bg.R.", fr$Id),
                             arrows = "", color = "#eeeeee", width = ghost_width, smooth = FALSE)
    }

    # Ghost stay bars (vertical)
    ghost_left  <- dot(paste0(".ghost.L.", ll$id), node_x[ll$id], ll$Adm * y_scale)
    ghost_right <- dot(paste0(".ghost.R.", ll$id), node_x[ll$id], dis_y)

    # PTest markers: short horizontal edges crossing each ghost bar
    ll_pt       <- ll[!is.na(ll$PTest), ]
    ptest_y     <- ll_pt$PTest * y_scale
    ptest_x     <- node_x[ll_pt$id]
    marker_top  <- dot(paste0(".pt.T.", ll_pt$id), ptest_x - ghost_width / 2, ptest_y)
    marker_bot  <- dot(paste0(".pt.B.", ll_pt$id), ptest_x + ghost_width / 2, ptest_y)
    marker_edges <- data.frame(
      from = paste0(".pt.T.", ll_pt$id), to = paste0(".pt.B.", ll_pt$id),
      arrows = "", color = "#222222", width = 1.5, smooth = FALSE)

    # Day tick labels to the left
    y_bottom   <- max(tick_days) * y_scale
    x_tick     <- -50L
    tick_nodes <- data.frame(
      id = paste0(".tick.", tick_days), label = as.character(tick_days),
      x = x_tick, y = tick_days * y_scale,
      shape = "text", font.size = 18, font.color = "#444444",
      stringsAsFactors = FALSE)

    # Bed labels below the plot
    label_nodes <- data.frame(
      id = paste0(".xlabel.", rownames(bed_pos)), label = rownames(bed_pos),
      x = seq_len(n_beds) * x_scale, y = y_bottom + 30L,
      shape = "text", font.size = 14, font.color = "#444444",
      stringsAsFactors = FALSE)

    # Guide line anchors (vertical, one per bed column)
    guide_L <- dot(paste0(".row.L.", seq_len(n_beds)), seq_len(n_beds) * x_scale, 0)
    guide_R <- dot(paste0(".row.R.", seq_len(n_beds)), seq_len(n_beds) * x_scale, y_bottom)

    # Axis spine anchors
    spine <- dot(c(".spine.a", ".spine.b"), c(x_scale, (n_beds + 1) * x_scale), y_bottom + 30L)

    bg_nodes <- if (show_all_stays) dplyr::bind_rows(bg_left, bg_right) else NULL
    out$x$nodes <- dplyr::bind_rows(nodes, tick_nodes, label_nodes,
                                    guide_L, guide_R, bg_nodes,
                                    ghost_left, ghost_right, spine,
                                    marker_top, marker_bot)

    axis_edge  <- data.frame(from = ".spine.a", to = ".spine.b",
                             arrows = "", color = "#aaaaaa", width = 1, smooth = FALSE)
    spine_edge <- data.frame(from = paste0(".tick.", tick_days[1]),
                             to   = paste0(".tick.", tick_days[length(tick_days)]),
                             arrows = "", color = "#aaaaaa", width = 1, smooth = FALSE)
    guide_edges <- data.frame(from = paste0(".row.L.", seq_len(n_beds)),
                              to   = paste0(".row.R.", seq_len(n_beds)),
                              arrows = "", color = "#dddddd", width = 0.5, smooth = FALSE)
  }

  # Ghost edges (same regardless of orientation)
  ghost_edges <- data.frame(
    from = paste0(".ghost.L.", ll$id), to = paste0(".ghost.R.", ll$id),
    arrows = "", color = hex_to_rgba(node_hex[ll$id], alpha = 0.35),
    width = ghost_width, smooth = FALSE)

  bg_e <- if (show_all_stays) bg_edges else NULL

  if (show_tree) {
    # Transmission edge color and width
    node_room <- setNames(ll$Room, ll$id)
    node_ward <- setNames(ll$Ward, ll$id)
    tx        <- out$x$edges
    tx$color  <- if (!is.null(uniform_edge_color)) {
      uniform_edge_color
    } else {
      ifelse(node_room[tx$from] == node_room[tx$to], edge_colors["room"],
      ifelse(node_ward[tx$from] == node_ward[tx$to], edge_colors["ward"],
             edge_colors["hospital"]))
    }
    tx$width  <- if (edge_width_by == "gen") {
      gen_scale * setNames(ll$Gen, ll$id)[tx$to]
    } else if (edge_width_by == "fixed") {
      gen_scale
    } else {
      gen_scale * (1 + as.numeric(tx[[edge_label]]))
    }
    out$x$edges <- dplyr::bind_rows(bg_e, ghost_edges, marker_edges, tx, axis_edge, spine_edge, guide_edges)
  } else {
    # Background only: drop case nodes and transmission edges
    out$x$nodes <- out$x$nodes[!out$x$nodes$id %in% ll$id, ]
    out$x$edges <- dplyr::bind_rows(bg_e, ghost_edges, marker_edges, axis_edge, spine_edge, guide_edges)
  }

  if (!show_tree) {
    out$x$legend <- NULL
    out$x$groups <- NULL
  }

  out |>
    visPhysics(enabled = FALSE) |>
    visNodes(fixed = list(x = TRUE, y = TRUE), font = list(size = node_font_size)) |>
    visEdges(font = list(size = edge_font_size)) |>
    visLegend(width = 0.1) |>
    visOptions(height = height)
}

# Save a visNetwork timeline plot to crisp PNG or vector PDF.
#
# Uses chromote directly (headless Chrome) for precise control:
#   - sets device pixel ratio for HiDPI output
#   - calls vis.js network.fit() to zoom the camera to the node bounding box,
#     replicating the "Viewer zoom" appearance
#   - waits for the fit to complete before screenshotting
#
# file:        filename; if no directory given, saved to figures/
# width/height: viewport size in CSS pixels
# scale:       device pixel ratio — 3 gives poster-quality resolution
# delay:       seconds to wait for the widget to finish its initial render
save_timeline <- function(plot, file, width = 1400, scale = 3, delay = 2.0) {
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("png", "pdf"))
    stop('file extension must be "png" or "pdf" (pdf gives vector output)')

  if (!grepl(.Platform$file.sep, file, fixed = TRUE)) {
    dir.create("figures", showWarnings = FALSE)
    file <- file.path("figures", file)
  }

  # Paper height must match the widget container height exactly.
  # Any HTML patching of the container breaks canvas rendering in printToPDF.
  height_px <- as.integer(gsub("[^0-9]", "", plot$x$height))

  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)

  htmlwidgets::saveWidget(plot, tmp, selfcontained = TRUE, title = "")

  html <- readLines(tmp, warn = FALSE)
  html <- sub("<head>",
    '<head><style>html,body{margin:0;padding:0;overflow:hidden;}</style>',
    html, fixed = TRUE)
  writeLines(html, tmp)

  b <- chromote::ChromoteSession$new()
  on.exit(b$close(), add = TRUE)

  b$Emulation$setDeviceMetricsOverride(
    width = width, height = height_px,
    deviceScaleFactor = scale, mobile = FALSE
  )

  b$Page$navigate(paste0("file://", normalizePath(tmp, mustWork = TRUE)))
  Sys.sleep(delay)
  b$Runtime$evaluate('document.readyState')
  Sys.sleep(0.4)

  if (ext == "png") {
    b$screenshot(filename = file)
    # Auto-trim uniform whitespace border
    if (requireNamespace("magick", quietly = TRUE))
      magick::image_read(file) |> magick::image_trim() |> magick::image_write(file)
  } else {
    # PDF: whitespace is unavoidable from vis.js's internal fit() padding;
    # crop manually in Illustrator/Inkscape or with pdfcrop (TeX Live)
    pdf_raw <- b$Page$printToPDF(
      paperWidth      = width     / 96,
      paperHeight     = height_px / 96,
      marginTop       = 0, marginBottom = 0, marginLeft = 0, marginRight = 0,
      printBackground = TRUE,
      pageRanges      = "1"
    )
    writeBin(jsonlite::base64_dec(pdf_raw$data), file)
  }

  invisible(file)
}
