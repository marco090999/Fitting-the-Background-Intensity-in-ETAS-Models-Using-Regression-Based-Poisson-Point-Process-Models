#######################################################################
# SENSITIVITY ANALYSIS SU M0: 3 MODELLI x 4 SOGLIE, CATALOGO REALE
#
# Risposta a R1-5 / R3-11 ("test sensitivity to several thresholds").
#
# Rifitta i tre modelli esatti usati per Tabella 1/2 del paper
# (etas.class.p1, etas.class.trig.p2 con flp=TRUE, etas.par.p9),
# facendo variare SOLO magn.threshold (M0) tra {2.5, 2.8, 3.0, 3.2},
# le soglie informate dall'analisi di completezza (Fase 1: MAXC
# globale=2.76, b-value stability globale=3.05, range per cella
# 2.76-2.93). magn.threshold.back resta fisso a 3.9 in ogni fit:
# verificato nel codice sorgente di etasclass()/etasclass.par() che
# questo parametro e' inerte data la convenzione di chiamata usata qui
# (w e hdef sempre forniti esplicitamente), quindi non introduce alcuna
# variazione confondente.
#
# DIPENDENZE: richiede in sessione le definizioni di etasclass() e
# etasclass.par() (dallo stesso script sorgente del piano di
# simulazione) e l'oggetto catalog.withcov.
#
# USO:
#   Sys.setenv(M0SENS_MODE = "run")
#   source("m0_sensitivity_application.R")
#######################################################################

suppressPackageStartupMessages(library(ggplot2))

m0s_log <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("etasclass", mode = "function") || !exists("etasclass.par", mode = "function")) {
  stop("Non trovo etasclass()/etasclass.par() in sessione. Carica prima le definizioni ",
       "dal file sorgente del piano di simulazione (es. etas_parametric_final_plan_no_smooth.R).")
}
if (!exists("catalog.withcov")) {
  stop("Non trovo 'catalog.withcov' in sessione.")
}

M0S_THRESHOLDS <- c(2.5, 2.8, 3.0, 3.2)
M0S_BACK <- 3.9   # fisso: verificato inerte, vedi commento in testa al file
M0S_OUT_DIR <- Sys.getenv("M0SENS_OUT_DIR", unset = "m0_sensitivity_outputs")


## ================================================================
## 1. FIT DEI TRE MODELLI A UNA DATA SOGLIA M0
## ================================================================

fit_three_models_at_m0 <- function(m0, cat_data = catalog.withcov, m0_back = M0S_BACK) {
  ## IMPORTANTE -- due lunghezze diverse per argomenti diversi, dovuto
  ## all'ordine interno di elaborazione di etasclass()/etasclass.par():
  ## - hvarx/hvary vengono cbind-ate al catalogo PRIMA del filtro per
  ##   magn.threshold (dentro cat.select(), che restituisce SEMPRE il
  ##   catalogo non filtrato con una colonna logica 'ind'); devono
  ##   quindi avere la lunghezza del catalogo ORIGINALE non filtrato.
  ## - w (e wp per etasclass.par) vengono invece validati DOPO che il
  ##   filtro e' stato applicato: devono avere la lunghezza del
  ##   catalogo GIA' FILTRATO per la soglia m0 corrente.
  ## A M0=2.5 le due lunghezze coincidono per caso (nessun evento viene
  ## filtrato), il che ha mascherato la distinzione finora.
  n_full <- nrow(cat_data)
  n_filtered <- sum(cat_data$magn1 >= m0, na.rm = TRUE)
  m0s_log("  Eventi attesi sopra M0=%.2f: n=%d (catalogo non filtrato: %d)", m0, n_filtered, n_full)

  common_hvarx <- replicate(n_full, 1)
  common_hvary <- replicate(n_full, 1)
  common_w <- replicate(n_filtered, 1)

  m0s_log("  Fitting ETAS classico (M0=%.2f)...", m0)
  fit_classic <- tryCatch(
    etasclass(cat.orig = cat_data, magn.threshold = m0, magn.threshold.back = m0_back,
             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
             w = common_w, hvarx = common_hvarx, hvary = common_hvary,
             formula1 = "time ~ magnitude - 1", declustering = TRUE,
             thinning = FALSE, flp = FALSE, ndeclust = 15, onlytime = FALSE,
             is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
             compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36),
    error = function(e) { m0s_log("    FALLITO: %s", conditionMessage(e)); NULL }
  )

  m0s_log("  Fitting ETAS-FLP (M0=%.2f)...", m0)
  fit_flp <- tryCatch(
    etasclass(cat.orig = cat_data, magn.threshold = m0, magn.threshold.back = m0_back,
             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5,
             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
             w = common_w, hvarx = common_hvarx, hvary = common_hvary,
             formula1 = "time ~ magnitude + distmin - 1", declustering = TRUE,
             thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
             is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
             compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36),
    error = function(e) { m0s_log("    FALLITO: %s", conditionMessage(e)); NULL }
  )

  m0s_log("  Fitting ETAS-P (M0=%.2f)...", m0)
  fit_param <- tryCatch(
    etasclass.par(cat.orig = cat_data, time.update = FALSE, magn.threshold = m0, magn.threshold.back = m0_back,
                 tmax = max(cat_data$time), long.range = range(cat_data$long), lat.range = range(cat_data$lat),
                 mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                 params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
                 formula1 = "time ~ magnitude + distmin - 1",
                 offset = 0, hdef = c(1, 1), w = common_w, hvarx = common_hvarx, hvary = common_hvary,
                 declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                 ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE,
                 description = "", cat.back = NULL, back.smooth = 1,
                 sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS",
                 compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                 formula.bg = ~ s(x, y, k = 25) + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                 type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"),
                 mult.bg = 4, ncube.bg = NULL, verbose.bg = FALSE, offset_k.bg = FALSE,
                 grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2),
    error = function(e) { m0s_log("    FALLITO: %s", conditionMessage(e)); NULL }
  )

  list(classic = fit_classic, flp = fit_flp, param = fit_param)
}


## ================================================================
## 2. ESTRAZIONE DEI PARAMETRI IN FORMATO LUNGO
## ================================================================

extract_trigger_params_m0 <- function(fit_obj, model_name, m0) {
  if (is.null(fit_obj) || is.null(fit_obj$params.MLtot)) return(data.frame())
  p <- fit_obj$params.MLtot
  se <- fit_obj$sqm %||% setNames(rep(NA_real_, length(p)), names(p))
  data.frame(
    model = model_name, m0 = m0, parameter = names(p),
    estimate = as.numeric(p), se = as.numeric(se[names(p)]),
    stringsAsFactors = FALSE
  )
}

extract_bg_coefs_m0 <- function(fit_obj, model_name, m0) {
  if (is.null(fit_obj) || is.null(fit_obj$model.bg$mod_global)) return(data.frame())
  mod <- fit_obj$model.bg$mod_global
  ptab <- tryCatch(summary(mod)$p.table, error = function(e) NULL)
  if (is.null(ptab)) return(data.frame())
  data.frame(
    model = model_name, m0 = m0, coefficient = rownames(ptab),
    estimate = ptab[, "Estimate"], se = ptab[, "Std. Error"],
    stringsAsFactors = FALSE
  )
}

extract_aic_m0 <- function(fits_list, m0) {
  get_final_aic <- function(fit) {
    if (is.null(fit) || is.null(fit$AIC.iter) || !length(fit$AIC.iter)) return(NA_real_)
    as.numeric(utils::tail(fit$AIC.iter, 1))
  }
  data.frame(
    model = c("classic", "flp", "param"), m0 = m0,
    aic = c(get_final_aic(fits_list$classic), get_final_aic(fits_list$flp), get_final_aic(fits_list$param)),
    stringsAsFactors = FALSE
  )
}


## ================================================================
## 3. ORCHESTRAZIONE: LOOP SULLE 4 SOGLIE
## ================================================================

run_m0_sensitivity <- function(thresholds = M0S_THRESHOLDS, cat_data = catalog.withcov,
                               out_dir = M0S_OUT_DIR, force = FALSE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  raw_fits_file <- file.path(out_dir, "m0_sensitivity_raw_fits.rds")

  ## Resume: se esistono gia' risultati salvati, riusali e completa
  ## solo le soglie mancanti o non ancora riuscite (a meno di force=TRUE).
  ## Ogni fit costa qui parecchi minuti (vedi log: ~20 min per un singolo
  ## modello a M0=2.5), quindi ripartire sempre da zero e' molto costoso.
  raw_fits <- list()
  if (!isTRUE(force) && file.exists(raw_fits_file)) {
    raw_fits <- readRDS(raw_fits_file)
    m0s_log("Trovati risultati precedenti in %s: %s", raw_fits_file,
           paste(names(raw_fits), collapse = ", "))
  }

  is_threshold_done <- function(m0) {
    key <- as.character(m0)
    if (is.null(raw_fits[[key]])) return(FALSE)
    fits <- raw_fits[[key]]
    ## considerata "fatta" solo se TUTTI e tre i modelli sono riusciti;
    ## altrimenti la si rifa' (potrebbe essere stata un fallimento parziale)
    !is.null(fits$classic) && !is.null(fits$flp) && !is.null(fits$param)
  }

  trig_rows <- list(); bg_rows <- list(); aic_rows <- list()

  for (m0 in thresholds) {
    if (is_threshold_done(m0) && !isTRUE(force)) {
      m0s_log("=== M0 = %.2f: gia' presente e completa, salto ===", m0)
    } else {
      m0s_log("=== M0 = %.2f ===", m0)
      fits <- fit_three_models_at_m0(m0, cat_data = cat_data)
      raw_fits[[as.character(m0)]] <- fits
      ## checkpoint incrementale: salva subito dopo ogni soglia, cosi'
      ## un fallimento successivo non fa perdere il lavoro gia' fatto
      saveRDS(raw_fits, raw_fits_file)
    }
  }

  ## ricostruisce le tabelle da TUTTE le soglie presenti in raw_fits
  ## (sia quelle appena calcolate sia quelle riprese da un run precedente)
  for (m0_key in names(raw_fits)) {
    m0 <- as.numeric(m0_key)
    fits <- raw_fits[[m0_key]]
    trig_rows[[length(trig_rows) + 1L]] <- extract_trigger_params_m0(fits$classic, "classic", m0)
    trig_rows[[length(trig_rows) + 1L]] <- extract_trigger_params_m0(fits$flp, "flp", m0)
    trig_rows[[length(trig_rows) + 1L]] <- extract_trigger_params_m0(fits$param, "param", m0)
    bg_rows[[length(bg_rows) + 1L]] <- extract_bg_coefs_m0(fits$param, "param", m0)
    aic_rows[[length(aic_rows) + 1L]] <- extract_aic_m0(fits, m0)
  }

  rb <- function(x) { x <- x[vapply(x, function(d) is.data.frame(d) && nrow(d) > 0, logical(1))]; if (length(x)) do.call(rbind, x) else data.frame() }
  results <- list(
    trigger_params = rb(trig_rows),
    bg_coefs = rb(bg_rows),
    aic = rb(aic_rows),
    raw_fits = raw_fits
  )

  saveRDS(results, file.path(out_dir, "m0_sensitivity_results.rds"))
  utils::write.csv(results$trigger_params, file.path(out_dir, "table_trigger_params_by_m0.csv"), row.names = FALSE)
  utils::write.csv(results$bg_coefs, file.path(out_dir, "table_bg_coefs_by_m0.csv"), row.names = FALSE)
  utils::write.csv(results$aic, file.path(out_dir, "table_aic_by_m0.csv"), row.names = FALSE)

  m0s_log("Fatto. Risultati salvati in: %s", normalizePath(out_dir))
  results
}


## ================================================================
## 4. GRAFICI
## ================================================================

MODEL_LABELS_M0 <- c(classic = "ETAS", flp = "ETAS-FLP", param = "ETAS-P")

## Traiettoria dei parametri di triggering chiave al variare di M0
plot_trigger_trajectory <- function(results, params_to_show = c("k0", "c", "p", "d", "q", "distmin")) {
  d <- results$trigger_params
  d <- d[d$parameter %in% params_to_show, ]
  if (!nrow(d)) return(NULL)
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))

  ggplot(d, aes(x = m0, y = estimate, color = model)) +
    geom_ribbon(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se, fill = model),
               alpha = 0.15, color = NA) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ parameter, scales = "free_y") +
    scale_color_discrete(labels = MODEL_LABELS_M0, name = NULL) +
    scale_fill_discrete(labels = MODEL_LABELS_M0, name = NULL) +
    labs(title = "Stabilita' dei parametri di triggering al variare di M0",
        subtitle = "Catalogo reale; bande = stima +/- 1.96 SE (condizionale)",
        x = "Soglia M0", y = "Stima") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}

## Traiettoria del coefficiente distmin nel background (solo ETAS-P)
plot_bg_coef_trajectory <- function(results) {
  d <- results$bg_coefs
  d <- d[d$coefficient == "distmin", ]
  if (!nrow(d)) return(NULL)

  ggplot(d, aes(x = m0, y = estimate)) +
    geom_ribbon(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se), alpha = 0.2, fill = "steelblue") +
    geom_line(color = "steelblue") + geom_point(size = 2, color = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = "Coefficiente distmin nel background (ETAS-P) al variare di M0",
        subtitle = "Banda = stima +/- 1.96 SE (condizionale, dal GAM)",
        x = "Soglia M0", y = expression(theta[distmin])) +
    theme_minimal(base_size = 12)
}

## AIC dei tre modelli al variare di M0 (solo per riferimento: il
## conteggio dei parametri effettivo differisce tra kernel/GAM, quindi
## il confronto assoluto tra modelli va interpretato con la stessa
## cautela discussata per la Tabella 1 del main text)
plot_aic_by_m0 <- function(results) {
  d <- results$aic
  if (!nrow(d) || all(is.na(d$aic))) return(NULL)
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))

  ggplot(d, aes(x = m0, y = aic, color = model)) +
    geom_line() + geom_point(size = 2) +
    scale_color_discrete(labels = MODEL_LABELS_M0, name = NULL) +
    labs(title = "AIC dei tre modelli al variare di M0",
        x = "Soglia M0", y = "AIC") +
    theme_minimal(base_size = 12)
}


## ================================================================
## 5. ESPLORAZIONE / RIASSUNTO
## ================================================================

print_m0_sensitivity_summary <- function(results) {
  cat("\n============================================================\n")
  cat("PARAMETRI DI TRIGGERING PER MODELLO E SOGLIA M0\n")
  cat("============================================================\n")
  d <- results$trigger_params
  print(d[order(d$model, d$parameter, d$m0), ], row.names = FALSE, digits = 4)

  cat("\n============================================================\n")
  cat("COEFFICIENTI DEL BACKGROUND (ETAS-P)\n")
  cat("============================================================\n")
  print(results$bg_coefs, row.names = FALSE, digits = 4)

  cat("\n============================================================\n")
  cat("AIC PER MODELLO E SOGLIA M0\n")
  cat("============================================================\n")
  print(results$aic, row.names = FALSE, digits = 6)

  ## variazione relativa del coefficiente distmin (triggering) tra
  ## M0=2.5 e le altre soglie, per modello -- la sintesi numerica piu'
  ## diretta per la lettera di risposta
  cat("\n============================================================\n")
  cat("VARIAZIONE RELATIVA DI 'distmin' (TRIGGERING) RISPETTO A M0=2.5\n")
  cat("============================================================\n")
  dt <- results$trigger_params[results$trigger_params$parameter == "distmin", ]
  for (mdl in unique(dt$model)) {
    sub <- dt[dt$model == mdl, ]
    base <- sub$estimate[sub$m0 == min(sub$m0)]
    sub$rel_change_vs_base <- (sub$estimate - base) / base
    cat(sprintf("\n-- %s --\n", MODEL_LABELS_M0[mdl] %||% mdl))
    print(sub[, c("m0", "estimate", "se", "rel_change_vs_base")], row.names = FALSE, digits = 3)
  }
  invisible(results)
}


## ================================================================
## 6. ESECUZIONE
## ================================================================

M0SENS_MODE <- tolower(Sys.getenv("M0SENS_MODE", unset = "none"))

if (identical(M0SENS_MODE, "run")) {
  m0_sensitivity_results <- run_m0_sensitivity()
  print_m0_sensitivity_summary(m0_sensitivity_results)

  p1 <- plot_trigger_trajectory(m0_sensitivity_results)
  p2 <- plot_bg_coef_trajectory(m0_sensitivity_results)
  p3 <- plot_aic_by_m0(m0_sensitivity_results)

  if (!is.null(p1)) { print(p1); ggsave(file.path(M0S_OUT_DIR, "fig_trigger_trajectory.png"), p1, width = 10, height = 6, dpi = 300) }
  if (!is.null(p2)) { print(p2); ggsave(file.path(M0S_OUT_DIR, "fig_bg_coef_trajectory.png"), p2, width = 7, height = 5, dpi = 300) }
  if (!is.null(p3)) { print(p3); ggsave(file.path(M0S_OUT_DIR, "fig_aic_by_m0.png"), p3, width = 7, height = 5, dpi = 300) }

} else {
  m0s_log("M0SENS_MODE='%s': funzioni caricate, nessun fit lanciato.", M0SENS_MODE)
}
