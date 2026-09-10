#######################################################################
# ANALISI DI COMPLETEZZA MAGNITUDO (Mc) DEL CATALOGO ITALIANO
#
# Risposta ai commenti R1-5 e R3-11 (soglia M0, distinta da Mc).
#
# METODI (standard, letteratura consolidata)
# --------------------------------------------
# 1) MAXC (Maximum Curvature): Mc = moda della distribuzione di
#    frequenza-magnitudo non cumulata, con correzione standard +0.2
#    (Wiemer & Wyss 2000; Woessner & Wiemer 2005).
# 2) Stabilita' del b-value (b-value stability / MBASS): per ogni
#    soglia candidata Mc*, si stima b(Mc*) (stimatore ML di
#    Aki/Utsu, con correzione di binning dM/2) e lo si confronta con
#    la media di b calcolata sulle 5 soglie successive (Mc*..Mc*+0.5);
#    Mc e' la soglia piu' bassa per cui questo scarto e' entro una
#    tolleranza (default 0.03) -- convenzione standard di ZMAP
#    (Wiemer & Wyss 2000; Cao & Gao 2002).
# Errore standard del b-value: formula di Shi & Bolt (1982).
#
# Entrambi i metodi sono applicati (a) all'intero catalogo, (b) a una
# griglia spaziale 4x4, per rispondere esplicitamente alla richiesta
# del revisore di verificare l'adeguatezza di M0 sia globalmente sia
# spazialmente.
#
# NOTA: M0 (soglia usata per la stima ETAS) e Mc (magnitudo di
# completezza del catalogo) sono concettualmente distinti (Sornette &
# Werner 2005a,b) -- qui stimiamo Mc in modo indipendente e lo
# confrontiamo con M0=2.5, non li confondiamo.
#
# DIPENDENZE: solo base R + ggplot2 (+ opzionalmente viridis e
# rnaturalearth/sf per la mappa, se disponibili). Serve in sessione
# l'oggetto catalog.withcov (colonne magn1, long, lat).
#
# USO:
#   Sys.setenv(COMPLETENESS_MODE = "run")   # default
#   source("catalog_completeness_analysis.R")
#######################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

cmp_log <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
cmp_dir_create <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
`%||%` <- function(a, b) if (!is.null(a)) a else b


## ================================================================
## 1. FUNZIONI DI STIMA: RISOLUZIONE MAGNITUDO, b-VALUE, Mc
## ================================================================

## Rileva la risoluzione (bin width) delle magnitudo, provando le
## precisioni piu' comuni nei cataloghi sismici. Stampa il valore
## trovato: controllare che sia sensato per i vostri dati.
detect_dm <- function(mags, candidates = c(0.01, 0.02, 0.1, 0.2)) {
  candidates <- sort(candidates, decreasing = TRUE)  # dal piu' grossolano al piu' fine
  for (dm in candidates) {
    resid <- abs(mags / dm - round(mags / dm))
    if (all(resid < 1e-6)) return(dm)
  }
  ## fallback: minima distanza osservata tra valori distinti ordinati
  u <- sort(unique(mags))
  d <- diff(u); d <- d[d > 1e-6]
  if (!length(d)) return(0.1)
  round(min(d), 3)
}

## Stimatore ML del b-value (Aki 1965 / Utsu 1965), con correzione di
## binning (Mc - dM/2), ed errore standard di Shi & Bolt (1982).
estimate_b_value <- function(mags, mc, dM) {
  m <- mags[mags >= mc - 1e-9]
  n <- length(m)
  if (n < 2) return(list(b = NA_real_, se = NA_real_, n = n, a_at_mc = NA_real_))
  mean_m <- mean(m)
  b <- log10(exp(1)) / (mean_m - (mc - dM / 2))
  se <- 2.30 * b^2 * sqrt(sum((m - mean_m)^2) / (n * (n - 1)))
  a_at_mc <- log10(n) + b * mc
  list(b = b, se = se, n = n, a_at_mc = a_at_mc)
}

## MAXC: moda della distribuzione di frequenza-magnitudo non cumulata,
## con la correzione standard +0.2 (Woessner & Wiemer 2005).
estimate_mc_maxc <- function(mags, dM, correction = 0.2) {
  if (length(mags) < 10) return(list(mc_raw = NA_real_, mc_corrected = NA_real_))
  breaks <- seq(floor(min(mags) / dM) * dM - dM / 2,
               ceiling(max(mags) / dM) * dM + dM / 2, by = dM)
  h <- hist(mags, breaks = breaks, plot = FALSE)
  mc_raw <- h$mids[which.max(h$counts)]
  list(mc_raw = mc_raw, mc_corrected = mc_raw + correction)
}

## Stabilita' del b-value (b-value stability / MBASS): per ciascuna
## soglia candidata, calcola lo scarto tra b(Mc) e la media di b sulle
## soglie successive entro forward_range unita' DI MAGNITUDO (non
## numero di passi!) -- default 0.5, la convenzione classica ZMAP con
## dM=0.1 (5 passi x 0.1 = 0.5). Con dM piu' fine (es. 0.01) il numero
## di passi corrispondenti viene scalato automaticamente: usare un
## numero fisso di PASSI invece che un range fisso in magnitudo
## restringerebbe la finestra "in avanti" fino a renderla sensibile al
## rumore locale invece che a un vero plateau.
estimate_mc_bvalue_stability <- function(mags, dM, mc_candidates = NULL,
                                         forward_range = 0.5, db_threshold = 0.03,
                                         min_n = 50) {
  if (is.null(mc_candidates)) {
    mc_candidates <- seq(floor(min(mags) / dM) * dM, quantile(mags, 0.95, na.rm = TRUE), by = dM)
  }
  rows <- lapply(mc_candidates, function(mc) {
    fit <- estimate_b_value(mags, mc, dM)
    data.frame(mc = mc, b = fit$b, se = fit$se, n = fit$n)
  })
  curve <- do.call(rbind, rows)

  curve$b_avg_forward <- NA_real_
  curve$db <- NA_real_
  for (i in seq_len(nrow(curve))) {
    fwd_idx <- which(curve$mc >= curve$mc[i] - 1e-9 & curve$mc <= curve$mc[i] + forward_range + 1e-9)
    if (length(fwd_idx) >= 2) {
      curve$b_avg_forward[i] <- mean(curve$b[fwd_idx], na.rm = TRUE)
      curve$db[i] <- abs(curve$b_avg_forward[i] - curve$b[i])
    }
  }

  ok <- curve[!is.na(curve$db) & curve$db <= db_threshold & curve$n >= min_n, ]
  mc_stability <- if (nrow(ok)) min(ok$mc) else NA_real_

  list(curve = curve, mc_stability = mc_stability)
}


## ================================================================
## 2. GRIGLIA SPAZIALE 4x4
## ================================================================

build_spatial_grid <- function(long, lat, nx = 4, ny = 4) {
  long_breaks <- seq(min(long, na.rm = TRUE), max(long, na.rm = TRUE), length.out = nx + 1)
  lat_breaks  <- seq(min(lat, na.rm = TRUE), max(lat, na.rm = TRUE), length.out = ny + 1)
  cell_x <- cut(long, breaks = long_breaks, include.lowest = TRUE, labels = FALSE)
  cell_y <- cut(lat, breaks = lat_breaks, include.lowest = TRUE, labels = FALSE)
  list(cell_x = cell_x, cell_y = cell_y, long_breaks = long_breaks, lat_breaks = lat_breaks)
}


## ================================================================
## 3. ANALISI COMPLETA: GLOBALE + PER CELLA
## ================================================================

run_completeness_analysis <- function(catalog, m0 = 2.5, nx = 4, ny = 4,
                                      min_n_cell = 50, dM = NULL,
                                      db_threshold = 0.03, forward_range = 0.5,
                                      maxc_correction = 0.2) {
  mags <- catalog$magn1
  if (is.null(dM)) dM <- detect_dm(mags)
  cmp_log("Risoluzione magnitudo rilevata (dM): %.3f", dM)

  ## ---- globale ----
  global_maxc <- estimate_mc_maxc(mags, dM, correction = maxc_correction)
  global_stability <- estimate_mc_bvalue_stability(mags, dM, db_threshold = db_threshold,
                                                    forward_range = forward_range, min_n = min_n_cell)
  global_b_at_m0 <- estimate_b_value(mags, m0, dM)

  mc_ref_global <- max(global_maxc$mc_corrected, global_stability$mc_stability, na.rm = TRUE)
  global_summary <- data.frame(
    region = "GLOBAL", n_total = length(mags), n_above_m0 = sum(mags >= m0),
    mc_maxc_raw = global_maxc$mc_raw, mc_maxc_corrected = global_maxc$mc_corrected,
    mc_bvalue_stability = global_stability$mc_stability,
    b_at_m0 = global_b_at_m0$b, se_b_at_m0 = global_b_at_m0$se,
    m0_adequate = m0 >= mc_ref_global,
    stringsAsFactors = FALSE
  )

  ## ---- per cella (griglia nx x ny) ----
  grid <- build_spatial_grid(catalog$long, catalog$lat, nx, ny)
  cell_rows <- list()
  for (iy in seq_len(ny)) {
    for (ix in seq_len(nx)) {
      idx <- which(grid$cell_x == ix & grid$cell_y == iy)
      cell_label <- sprintf("cell_%d_%d", ix, iy)
      base_row <- data.frame(
        region = cell_label, cell_x = ix, cell_y = iy,
        long_min = grid$long_breaks[ix], long_max = grid$long_breaks[ix + 1],
        lat_min = grid$lat_breaks[iy], lat_max = grid$lat_breaks[iy + 1],
        n_total = length(idx), n_above_m0 = sum(mags[idx] >= m0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      if (length(idx) < min_n_cell) {
        cell_rows[[cell_label]] <- cbind(base_row, data.frame(
          mc_maxc_raw = NA_real_, mc_maxc_corrected = NA_real_, mc_bvalue_stability = NA_real_,
          b_at_m0 = NA_real_, se_b_at_m0 = NA_real_, m0_adequate = NA,
          note = "dati_insufficienti", stringsAsFactors = FALSE
        ))
        next
      }

      m_cell <- mags[idx]
      maxc_cell <- estimate_mc_maxc(m_cell, dM, correction = maxc_correction)
      stab_cell <- estimate_mc_bvalue_stability(m_cell, dM, db_threshold = db_threshold,
                                                forward_range = forward_range, min_n = min_n_cell)
      b_cell <- estimate_b_value(m_cell, m0, dM)
      mc_ref_cell <- max(maxc_cell$mc_corrected, stab_cell$mc_stability, na.rm = TRUE)

      cell_rows[[cell_label]] <- cbind(base_row, data.frame(
        mc_maxc_raw = maxc_cell$mc_raw, mc_maxc_corrected = maxc_cell$mc_corrected,
        mc_bvalue_stability = stab_cell$mc_stability,
        b_at_m0 = b_cell$b, se_b_at_m0 = b_cell$se,
        m0_adequate = if (is.finite(mc_ref_cell)) m0 >= mc_ref_cell else NA,
        note = NA_character_, stringsAsFactors = FALSE
      ))
    }
  }
  cell_summary <- do.call(rbind, cell_rows)
  rownames(cell_summary) <- NULL

  list(
    dM = dM, m0 = m0, nx = nx, ny = ny, min_n_cell = min_n_cell,
    global = global_summary, by_cell = cell_summary,
    global_stability_curve = global_stability$curve,
    mags = mags, catalog = catalog, grid = grid
  )
}


## ================================================================
## 4. GRAFICI
## ================================================================

## Distribuzione di frequenza-magnitudo (cumulata e non), con la retta
## di Gutenberg-Richter fittata sopra M0 e le stime di Mc annotate.
plot_fmd <- function(result) {
  mags <- result$mags; dM <- result$dM; m0 <- result$m0
  breaks <- seq(floor(min(mags) / dM) * dM - dM / 2,
               ceiling(max(mags) / dM) * dM + dM / 2, by = dM)
  h <- hist(mags, breaks = breaks, plot = FALSE)
  df <- data.frame(mag = h$mids, noncum = h$counts, cum = rev(cumsum(rev(h$counts))))
  df_long <- rbind(
    data.frame(mag = df$mag, count = df$noncum, type = "Non-cumulative"),
    data.frame(mag = df$mag, count = df$cum, type = "Cumulative")
  )
  df_long <- df_long[df_long$count > 0, ]

  b_fit <- estimate_b_value(mags, m0, dM)
  gr_line <- data.frame(mag = seq(m0, max(mags), length.out = 50))
  gr_line$count <- 10^(b_fit$a_at_mc - b_fit$b * gr_line$mag)

  ggplot(df_long, aes(x = mag, y = count, color = type, shape = type)) +
    geom_point(size = 2) +
    geom_line(data = gr_line, aes(x = mag, y = count), color = "black",
             linetype = "dashed", inherit.aes = FALSE) +
    geom_vline(xintercept = m0, linetype = "dotted", color = "red") +
    geom_vline(xintercept = result$global$mc_maxc_corrected, linetype = "dotted", color = "darkgreen") +
    geom_vline(xintercept = result$global$mc_bvalue_stability, linetype = "dotted", color = "blue") +
    scale_y_log10() +
    labs(title = "Frequency-magnitude distribution (intero catalogo)",
        subtitle = sprintf("b(M0=%.2f) = %.3f +/- %.3f | Mc MAXC (corretto, verde) = %.2f | Mc b-stability (blu) = %.2f | M0 (rosso) = %.2f",
                           m0, b_fit$b, b_fit$se, result$global$mc_maxc_corrected,
                           result$global$mc_bvalue_stability, m0),
        x = "Magnitudo", y = "Numero di eventi (scala log)", color = NULL, shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9))
}

## Curva di stabilita' del b-value: b(Mc) +/- SE al variare della
## soglia candidata, con M0 e la soglia di stabilita' stimata annotate.
plot_bvalue_stability <- function(result) {
  curve <- result$global_stability_curve
  ggplot(curve, aes(x = mc, y = b)) +
    geom_ribbon(aes(ymin = b - se, ymax = b + se), alpha = 0.15) +
    geom_line() + geom_point(size = 1.5) +
    geom_vline(xintercept = result$m0, linetype = "dotted", color = "red") +
    geom_vline(xintercept = result$global$mc_bvalue_stability, linetype = "dashed", color = "blue") +
    labs(title = "Stabilita' del b-value al variare della soglia candidata",
        subtitle = sprintf("Rosso tratteggiato = M0 (%.2f) | Blu tratteggiato = Mc stimata da stabilita' (%.2f)",
                           result$m0, result$global$mc_bvalue_stability),
        x = "Soglia candidata Mc", y = "b-value (+/- SE)") +
    theme_minimal(base_size = 12)
}

## Mappa a griglia 4x4 dei valori di Mc stimati per cella, con annotazione
## del numero di eventi; celle con dati insufficienti mostrate in grigio.
plot_spatial_mc_grid <- function(result, fill_var = c("mc_maxc_corrected", "mc_bvalue_stability")) {
  fill_var <- match.arg(fill_var)
  df <- result$by_cell
  df$fill_val <- df[[fill_var]]
  df$label <- ifelse(is.na(df$fill_val),
                     sprintf("n.d.\n(n=%d)", df$n_total),
                     sprintf("%.2f\n(n=%d)", df$fill_val, df$n_total))

  p <- ggplot(df) +
    geom_rect(aes(xmin = long_min, xmax = long_max, ymin = lat_min, ymax = lat_max, fill = fill_val),
             color = "grey40") +
    geom_text(aes(x = (long_min + long_max) / 2, y = (lat_min + lat_max) / 2, label = label),
             size = 3, color = "black")

  ## overlay costa italiana se rnaturalearth/sf sono disponibili (opzionale)
  if (requireNamespace("rnaturalearth", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    italy_sf <- tryCatch(rnaturalearth::ne_countries(country = "Italy", scale = "medium", returnclass = "sf"),
                         error = function(e) NULL)
    if (!is.null(italy_sf)) {
      p <- p + ggplot2::geom_sf(data = italy_sf, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE)
    }
  }

  p +
    scale_fill_viridis_c(name = fill_var, option = "plasma", na.value = "grey85") +
    coord_sf(xlim = range(df$long_min, df$long_max), ylim = range(df$lat_min, df$lat_max), expand = FALSE) +
    labs(title = sprintf("Variazione spaziale di Mc (griglia %dx%d) -- %s", result$nx, result$ny, fill_var),
        subtitle = sprintf("M0 = %.2f. Celle con Mc > M0 sono potenzialmente inadeguate (vedi tabella per il dettaglio numerico)", result$m0),
        x = "Longitudine", y = "Latitudine") +
    theme_minimal(base_size = 12)
}


## ================================================================
## 5. OUTPUT: TABELLE, SALVATAGGIO, RIASSUNTO TESTUALE
## ================================================================

print_completeness_summary <- function(result) {
  cat("\n============================================================\n")
  cat("RISULTATI GLOBALI\n")
  cat("============================================================\n")
  print(result$global, row.names = FALSE, digits = 3)

  cat("\n============================================================\n")
  cat(sprintf("RISULTATI PER CELLA (griglia %dx%d, min_n=%d)\n", result$nx, result$ny, result$min_n_cell))
  cat("============================================================\n")
  print(result$by_cell[, c("region", "n_total", "n_above_m0", "mc_maxc_corrected",
                           "mc_bvalue_stability", "m0_adequate", "note")],
       row.names = FALSE, digits = 3)

  n_inadequate <- sum(!isTRUE(result$by_cell$m0_adequate) & is.finite(result$by_cell$mc_maxc_corrected) |
                       (!is.na(result$by_cell$m0_adequate) & !result$by_cell$m0_adequate), na.rm = TRUE)
  n_insufficient <- sum(result$by_cell$note == "dati_insufficienti", na.rm = TRUE)
  n_ok <- sum(result$by_cell$m0_adequate, na.rm = TRUE)

  cat(sprintf("\nSintesi: %d/%d celle con dati sufficienti; tra queste, %d con M0 adeguata, %d con M0 potenzialmente inadeguata.\n",
             result$nx * result$ny - n_insufficient, result$nx * result$ny, n_ok,
             (result$nx * result$ny - n_insufficient) - n_ok))
  if (n_insufficient > 0) {
    cat(sprintf("%d celle avevano meno di %d eventi e sono state escluse dalla stima locale.\n",
               n_insufficient, result$min_n_cell))
  }
  invisible(result)
}

save_completeness_results <- function(result, out_dir = "completeness_analysis_outputs") {
  cmp_dir_create(out_dir)
  cmp_dir_create(file.path(out_dir, "figures"))

  utils::write.csv(result$global, file.path(out_dir, "table_global_completeness.csv"), row.names = FALSE)
  utils::write.csv(result$by_cell, file.path(out_dir, "table_by_cell_completeness.csv"), row.names = FALSE)
  utils::write.csv(result$global_stability_curve, file.path(out_dir, "table_bvalue_stability_curve.csv"), row.names = FALSE)
  saveRDS(result, file.path(out_dir, "completeness_result_full.rds"))

  p1 <- plot_fmd(result)
  p2 <- plot_bvalue_stability(result)
  p3 <- plot_spatial_mc_grid(result, "mc_maxc_corrected")
  p4 <- plot_spatial_mc_grid(result, "mc_bvalue_stability")

  ggsave(file.path(out_dir, "figures", "fig_fmd.png"), p1, width = 8, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "figures", "fig_bvalue_stability.png"), p2, width = 8, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "figures", "fig_spatial_mc_maxc.png"), p3, width = 8, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "figures", "fig_spatial_mc_bvalue.png"), p4, width = 8, height = 7, dpi = 300)

  cmp_log("Tabelle e grafici salvati in: %s", normalizePath(out_dir))
  invisible(list(fmd = p1, bvalue_stability = p2, spatial_maxc = p3, spatial_bvalue = p4))
}


## ================================================================
## 6. ESECUZIONE
## ================================================================

COMPLETENESS_MODE <- tolower(Sys.getenv("COMPLETENESS_MODE", unset = "run"))
COMPLETENESS_OUT_DIR <- Sys.getenv("COMPLETENESS_OUT_DIR", unset = "completeness_analysis_outputs")

if (identical(COMPLETENESS_MODE, "run")) {
  if (!exists("catalog.withcov")) {
    stop("Non trovo 'catalog.withcov' in sessione. Caricalo (library(etasFLP); data(catalog.withcov)) prima di sourceare questo script.")
  }
  completeness_result <- run_completeness_analysis(
    catalog = catalog.withcov, m0 = 2.5, nx = 4, ny = 4, min_n_cell = 50
  )
  print_completeness_summary(completeness_result)
  plots_completeness <- save_completeness_results(completeness_result, out_dir = COMPLETENESS_OUT_DIR)

} else {
  cmp_log("COMPLETENESS_MODE='%s': funzioni caricate, nessuna analisi lanciata.", COMPLETENESS_MODE)
}
