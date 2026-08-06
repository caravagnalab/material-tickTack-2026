rm(list = ls())
library(tidyverse)
library(BradleyTerry2)
library(data.table)
library(parallel)

# ── Constants ──────────────────────────────────────────────────────────────────
MIN_FRAC           <- 0.10
TIE_FACTOR         <- 0.5
EPS                <- 1e-3
MIN_PAIR_OBS       <- 3
MIN_CONNECTIONS    <- 3
CI_WIDTH_THRESHOLD <- 4

# ── Output directories ─────────────────────────────────────────────────────────
dir.res         <- "../Analysis_results/BTM"
sub.dir.res     <- file.path(dir.res, "drivers")
sub.dir.res.arm <- file.path(dir.res, "arm_level")
dir.create(sub.dir.res,     recursive = TRUE)
dir.create(sub.dir.res.arm, recursive = TRUE)

# ── Helper: build pairwise win/loss table ──────────────────────────────────────
build_pairwise <- function(dt, items, tie_factor, eps) {
  all_comb   <- combn(seq_along(items), 2)
  pairs_list <- lapply(seq_len(ncol(all_comb)), function(i) sort(items[all_comb[, i]]))
  
  n_cores <- max(1L, round(parallel::detectCores() / 2))
  
  results <- parallel::mclapply(pairs_list, function(pair) {
    p1 <- pair[1L]; p2 <- pair[2L]
    d1 <- dt[.(p1), .(sample_id, r1 = clock_rank)]
    d2 <- dt[.(p2), .(sample_id, r2 = clock_rank)]
    m  <- d1[d2, on = "sample_id", nomatch = NULL]
    if (nrow(m) == 0L) return(NULL)
    m[, `:=`(win  = r1 > r2,
             tie  = r1 == r2,
             loss = r1 < r2)]
    m[, .(p1 = p1, p2 = p2,
          # NOTE: round() removed — BTm accepts non-integer weights, and
          # rounding distorts counts when ties are common
          w1 = sum(win)  + tie_factor * sum(tie) + eps,
          w2 = sum(loss) + tie_factor * sum(tie) + eps)]
  }, mc.cores = n_cores)
  
  data.table::rbindlist(results, use.names = TRUE)
}

# ── Helper: filter pairwise table for sparsity ────────────────────────────────
filter_pairwise <- function(df_summarized, min_pair_obs, min_connections, eps) {
  # Drop pairs with too few co-occurrences
  df_summarized <- df_summarized %>%
    dplyr::filter(w1 + w2 >= min_pair_obs + eps)
  
  if (nrow(df_summarized) == 0) return(NULL)
  
  # Drop poorly connected items (appear in too few pairs)
  item_connections <- table(c(df_summarized$p1, df_summarized$p2))
  well_connected   <- names(item_connections[item_connections >= min_connections])
  
  df_summarized <- df_summarized %>%
    dplyr::filter(p1 %in% well_connected, p2 %in% well_connected)
  
  if (nrow(df_summarized) == 0) return(NULL)
  
  df_summarized
}

# ── Helper: fit BT model, return abilities + raw fit ──────────────────────────
fit_bt_model <- function(df_summarized, id_col = "gene") {
  all_items     <- sort(unique(c(df_summarized$p1, df_summarized$p2)))
  df_summarized <- df_summarized %>%
    mutate(p1 = factor(p1, levels = all_items),
           p2 = factor(p2, levels = all_items))
  
  fit <- tryCatch(
    BTm(cbind(w1, w2), p1, p2, data = df_summarized),
    error = function(e) {
      message("  BTm failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  
  abilities <- BTabilities(fit) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(id_col) %>%
    mutate(conf_low  = ability - 1.96 * s.e.,
           conf_high = ability + 1.96 * s.e.) %>%
    arrange(ability)
  
  list(abilities = abilities, fit = fit)
}

# ── Helper: flag wide CIs ─────────────────────────────────────────────────────
flag_unreliable <- function(abilities_df, ci_width_threshold, id_col) {
  abilities_df <- abilities_df %>%
    dplyr::mutate(ci_width   = conf_high - conf_low,
                  unreliable = ci_width > ci_width_threshold)
  
  if (any(abilities_df$unreliable)) {
    message("  Warning: ", sum(abilities_df$unreliable),
            " item(s) with CI width > ", ci_width_threshold, ": ",
            paste(abilities_df[[id_col]][abilities_df$unreliable], collapse = ", "))
  }
  
  abilities_df
}

# ── Helper: forest plot ───────────────────────────────────────────────────────
make_forest_plot <- function(abilities_df, id_col, tumour_type) {
  ggplot(abilities_df,
         aes(x = reorder(.data[[id_col]], ability), y = ability)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                  width = 0.3, color = "steelblue") +
    geom_point(aes(color = unreliable), size = 2.5) +
    scale_color_manual(
      values = c("FALSE" = "firebrick", "TRUE" = "gray60"),
      labels = c("FALSE" = "reliable",  "TRUE" = "sparse"),
      name   = NULL
    ) +
    coord_flip() +
    labs(
      title    = paste0("Mutational timeline — ", tumour_type),
      subtitle = "Higher scores = later events  |  Lower scores = earlier events",
      x        = id_col,
      y        = "Relative temporal ability (BT log-odds)"
    ) +
    theme_minimal()
}

# ── Helper: build and save metadata ──────────────────────────────────────────
make_meta <- function(tumour_type, n_samples, n_items_initial,
                      abilities_df, df_summarized) {
  list(
    tumour_type     = tumour_type,
    n_samples       = n_samples,
    n_items_initial = n_items_initial,
    n_items_fit     = nrow(abilities_df),
    n_pairs_fit     = nrow(df_summarized),
    n_unreliable    = sum(abilities_df$unreliable),
    thresholds      = list(
      MIN_FRAC           = MIN_FRAC,
      MIN_PAIR_OBS       = MIN_PAIR_OBS,
      MIN_CONNECTIONS    = MIN_CONNECTIONS,
      CI_WIDTH_THRESHOLD = CI_WIDTH_THRESHOLD
    ),
    timestamp = Sys.time()
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# Gene-level events
# ══════════════════════════════════════════════════════════════════════════════
df <- readRDS("../Data/Drivers.rds") %>%
  dplyr::filter(class != "WGD", timed)

# Exploratory histogram — uncomment to review before full run
# df %>%
#   dplyr::filter(ttype == "BRCA") %>%
#   dplyr::group_by(sample_id) %>%
#   dplyr::summarise(k = max(best_K)) %>%
#   ggplot(aes(x = k)) + geom_histogram()

results_gene <- lapply(unique(df$ttype), function(tumour_type) {
  message("Gene-level BT: ", tumour_type)
  
  df_tt <- df %>%
    dplyr::filter(ttype == tumour_type, class != "WGD") %>%
    dplyr::group_by(sample_id) %>%
    # dense_rank within sample: ties get equal rank, preserving the co-occurrence
    # structure. Cross-sample comparison relies on relative rank, not absolute clock.
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>%
    dplyr::ungroup()
  
  n_samples <- dplyr::n_distinct(df_tt$sample_id)
  
  frequent_genes <- df_tt %>%
    dplyr::distinct(sample_id, gene) %>%
    dplyr::count(gene) %>%
    dplyr::filter(n / n_samples >= MIN_FRAC) %>%
    dplyr::pull(gene)
  
  if (length(frequent_genes) <= 5) {
    message("  Skipped — fewer than 6 frequent genes")
    return(NULL)
  }
  
  dt <- df_tt %>%
    dplyr::filter(gene %in% frequent_genes) %>%
    dplyr::select(sample_id, gene, clock_rank) %>%
    dplyr::distinct() %>%
    data.table::as.data.table()
  data.table::setkey(dt, gene, sample_id)
  
  df_summarized <- build_pairwise(dt, frequent_genes, TIE_FACTOR, EPS)
  
  if (nrow(df_summarized) == 0) {
    message("  Skipped — no co-occurring pairs")
    return(NULL)
  }
  
  df_summarized <- filter_pairwise(df_summarized, MIN_PAIR_OBS, MIN_CONNECTIONS, EPS)
  
  if (is.null(df_summarized)) {
    message("  Skipped — no pairs surviving sparsity filters")
    return(NULL)
  }
  
  bt_out <- fit_bt_model(df_summarized, id_col = "gene")
  if (is.null(bt_out)) return(NULL)
  
  abilities_df <- flag_unreliable(bt_out$abilities, CI_WIDTH_THRESHOLD, id_col = "gene")
  meta         <- make_meta(tumour_type, n_samples, length(frequent_genes),
                            abilities_df, df_summarized)
  
  p <- make_forest_plot(abilities_df, id_col = "gene", tumour_type = tumour_type)
  
  ggsave(file.path(sub.dir.res, paste0(tumour_type, ".pdf")),
         plot = p, width = 8, height = min(0.25 * nrow(abilities_df), 48))
  
  saveRDS(abilities_df,      file.path(sub.dir.res, paste0(tumour_type, "_abilities.rds")))
  saveRDS(df_summarized,     file.path(sub.dir.res, paste0(tumour_type, "_pairwise.rds")))
  saveRDS(bt_out$fit,        file.path(sub.dir.res, paste0(tumour_type, "_fit.rds")))
  saveRDS(meta,              file.path(sub.dir.res, paste0(tumour_type, "_meta.rds")))
  
  list(plot = p, abilities = abilities_df, fit = bt_out$fit,
       pairwise = df_summarized, meta = meta)
})
names(results_gene) <- unique(df$ttype)

# ══════════════════════════════════════════════════════════════════════════════
# Arm-level events
# ══════════════════════════════════════════════════════════════════════════════
df_arm <- readRDS("../Data/Segments.rds") %>%
  dplyr::filter(class != "WGD")

results_arm <- lapply(unique(df_arm$ttype), function(tumour_type) {
  message("Arm-level BT: ", tumour_type)
  
  df_tt <- df_arm %>%
    dplyr::filter(ttype == tumour_type) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>%
    dplyr::ungroup()
  
  n_samples <- dplyr::n_distinct(df_tt$sample)
  
  frequent_arm_levels <- df_tt %>%
    dplyr::distinct(sample, chr, arm_level_class) %>%
    dplyr::count(chr, arm_level_class) %>%
    dplyr::filter(n / n_samples >= MIN_FRAC,
                  arm_level_class != "focal") %>%
    dplyr::mutate(event = paste0(chr, "_", arm_level_class)) %>%
    dplyr::pull(event)
  
  if (length(frequent_arm_levels) <= 3) {
    message("  Skipped — fewer than 4 frequent arm-level events")
    return(NULL)
  }
  
  dt <- df_tt %>%
    dplyr::mutate(event = paste0(chr, "_", arm_level_class)) %>%
    dplyr::filter(event %in% frequent_arm_levels) %>%
    dplyr::select(sample_id = sample, event, clock_rank) %>%
    dplyr::distinct() %>%
    data.table::as.data.table()
  data.table::setkey(dt, event, sample_id)
  
  df_summarized <- build_pairwise(dt, frequent_arm_levels, TIE_FACTOR, EPS)
  
  if (nrow(df_summarized) == 0) {
    message("  Skipped — no co-occurring pairs")
    return(NULL)
  }
  
  df_summarized <- filter_pairwise(df_summarized, MIN_PAIR_OBS, MIN_CONNECTIONS, EPS)
  
  if (is.null(df_summarized)) {
    message("  Skipped — no pairs surviving sparsity filters")
    return(NULL)
  }
  
  bt_out <- fit_bt_model(df_summarized, id_col = "event")
  if (is.null(bt_out)) return(NULL)
  
  abilities_df <- flag_unreliable(bt_out$abilities, CI_WIDTH_THRESHOLD, id_col = "event")
  meta         <- make_meta(tumour_type, n_samples, length(frequent_arm_levels),
                            abilities_df, df_summarized)
  
  p <- make_forest_plot(abilities_df, id_col = "event", tumour_type = tumour_type)
  
  ggsave(file.path(sub.dir.res.arm, paste0(tumour_type, "_arm_level.pdf")),
         plot = p, width = 6, height = max(2, nrow(abilities_df) * 0.25))
  
  saveRDS(abilities_df,      file.path(sub.dir.res.arm, paste0(tumour_type, "_abilities.rds")))
  saveRDS(df_summarized,     file.path(sub.dir.res.arm, paste0(tumour_type, "_pairwise.rds")))
  saveRDS(bt_out$fit,        file.path(sub.dir.res.arm, paste0(tumour_type, "_fit.rds")))
  saveRDS(meta,              file.path(sub.dir.res.arm, paste0(tumour_type, "_meta.rds")))
  
  list(plot = p, abilities = abilities_df, fit = bt_out$fit,
       pairwise = df_summarized, meta = meta)
})
names(results_arm) <- unique(df_arm$ttype)