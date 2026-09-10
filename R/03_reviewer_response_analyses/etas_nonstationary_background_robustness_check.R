#######################################################################
# ROBUSTNESS CHECK: TEMPORAL NONSTATIONARITY IN THE ETAS BACKGROUND
#
# Risposta ai punti 1 e 2 di Reviewer 1.
#
# COSA FA QUESTO SCRIPT
# ----------------------
# Genera cataloghi sintetici da un processo ETAS in cui il tasso di
# background NON e' costante nel tempo ma segue un trend lineare:
#
#     lambda(s,t) = mu(t) * f0(s)  +  triggering
#     mu(t) = mu0 * g(t),   g(t) = 1 + kappa * (t - tmid) / tau
#
# con g(t) centrata in modo che la sua media su [tmin,tmax] sia
# ESATTAMENTE 1 per qualunque kappa. Questo garantisce che:
#   - il numero atteso di eventi di background resti mu0 * tau, cioe'
#     la stessa quantita' che gia' calibrate per gli scenari stazionari
#     (make_count_regimes_full / count_lookup) resta valida SENZA
#     modifiche;
#   - la forma spaziale marginale vera resta ESATTAMENTE f0(s), quindi
#     le funzioni di verita'/metrica spaziale gia' esistenti
#     (make_grid_eval, true_background_grid_df, estimate_*_grid,
#     fit_background_metrics, compare_background_shape_param) restano
#     riutilizzabili senza modifiche.
#
# I tre stimatori (ETAS classico, ETAS-FLP, ETAS-P) vengono fittati
# SENZA ALCUNA MODIFICA rispetto agli script originali: tutti e tre
# assumono background omogeneo nel tempo. Lo scopo e' misurare, non
# correggere, l'effetto di questa assunzione violata.
#
# Il piano copre due scenari spaziali (per isolare l'effetto puramente
# temporale da quello combinato con una covariata spaziale strutturata,
# analoga al ruolo di "distanza dalla faglia" nell'applicazione reale):
#   - "constant" : f0(s) piatta
#   - "cov1"     : f0(s) = stessa covariata Z1 gia' usata nel piano
#                  principale (sin/cos), stesso beta vero
#
# DIPENDENZE
# ----------
# Questo file NON duplica le funzioni di stima/fit gia' scritte. Carica
# (senza autorun) le definizioni da:
#   1) etas_parametric_final_plan_no_smooth.R
#        -> etasclass(), etasclass.par(), fit_etas_classic_full(),
#           fit_etas_parametric_full(), background_fit_spec(),
#           make_grid_eval(), estimate_classic_background_grid(),
#           estimate_parametric_background_grid(),
#           true_background_grid_df(), fit_parameter_comparison(),
#           fit_trigger_components(), fit_classification_metrics(),
#           get_cat_sim(), get_sim_core(), make_eqcat_from_sim(),
#           saveRDS_atomic(), safe_run(), TRUE_BASE_FULL, ecc.
#   2) add_etas_flp_to_existing_results_NO_OVERWRITE_v3.R
#        -> fit_etas_flp_full() (= etasclass(..., flp = TRUE, ...))
#
# USO CONSIGLIATO
# ----------------
#   setwd("<cartella con i due script originali>")
#   Sys.setenv(NONSTAT_MODE = "test")   # 1 replica, per verificare
#   source("etas_nonstationary_background_robustness_check.R")
#
#   Sys.setenv(NONSTAT_MODE = "run")    # piano completo
#   source("etas_nonstationary_background_robustness_check.R")
#
#   Sys.setenv(NONSTAT_MODE = "collect_and_plot")  # solo analisi/grafici
#   source("etas_nonstationary_background_robustness_check.R")
#
# Variabili opzionali:
#   NONSTAT_BASE_DIR        default: getwd()
#   NONSTAT_PLAN_SCRIPT      default: etas_parametric_final_plan_no_smooth.R in BASE_DIR
#   NONSTAT_FLP_SCRIPT       default: add_etas_flp_to_existing_results_NO_OVERWRITE_v3.R in BASE_DIR
#   NONSTAT_OUT_ROOT         default: BASE_DIR/etas_nonstationary_outputs_v1
#   NONSTAT_KAPPA            default: 1.6  (rapporto fine/inizio finestra = 9x)
#   NONSTAT_NREP             default: 30
#   NONSTAT_NBINS            default: 10   (bin temporali per la diagnostica)
#######################################################################


## ================================================================
## 0. UTILITY DI BASE (indipendenti, come nello script FLP)
## ================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

ns_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ns_log <- function(..., .level = "INFO") {
  cat(sprintf("[%s] [%s] ", ns_now(), .level), sprintf(...), "\n", sep = "")
  flush.console()
}

ns_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ns_env_chr <- function(name, default = NULL) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) default else x
}

ns_env_num <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) default else out
}

ns_env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) default else out
}

ns_append_csv_row <- function(row, file) {
  ns_dir_create(dirname(file))
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  utils::write.table(
    row, file = file, sep = ",", row.names = FALSE,
    col.names = !file.exists(file), append = file.exists(file)
  )
  invisible(file)
}

ns_saveRDS_atomic <- function(object, file, compress = "gzip") {
  ns_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_ns_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) { file.copy(tmp, file, overwrite = TRUE); unlink(tmp) }
  invisible(file)
}


## ================================================================
## 1. CARICAMENTO SICURO DELLE DEFINIZIONI (nessun autorun)
## ================================================================

NONSTAT_BASE_DIR <- ns_env_chr("NONSTAT_BASE_DIR", getwd())

ns_guess_script <- function(explicit_env, candidates) {
  explicit <- ns_env_chr(explicit_env, NULL)
  all_cand <- unique(c(explicit, candidates))
  all_cand <- all_cand[!is.na(all_cand) & nzchar(all_cand)]
  hit <- all_cand[file.exists(all_cand)]
  if (!length(hit)) {
    stop("File non trovato. Provati: ", paste(all_cand, collapse = "\n  "),
         "\nImposta Sys.setenv(", explicit_env, " = '/percorso/completo/file.R').")
  }
  normalizePath(hit[1], mustWork = TRUE)
}

NONSTAT_PLAN_SCRIPT <- ns_guess_script(
  "NONSTAT_PLAN_SCRIPT",
  file.path(NONSTAT_BASE_DIR, "etas_parametric_final_plan_no_smooth.R")
)

NONSTAT_FLP_SCRIPT <- ns_guess_script(
  "NONSTAT_FLP_SCRIPT",
  file.path(NONSTAT_BASE_DIR, "add_etas_flp_to_existing_results_NO_OVERWRITE_v3.R")
)

## Carica SOLO le definizioni del piano principale (fino all'AUTORUN SWITCH),
## esattamente con lo stesso meccanismo gia' usato nello script FLP.
ns_load_plan_definitions <- function(plan_script = NONSTAT_PLAN_SCRIPT, envir = .GlobalEnv) {
  if (exists("etasclass.par", envir = envir, inherits = TRUE) &&
      exists("fit_etas_parametric_full", envir = envir, inherits = TRUE) &&
      exists("COUNT_MAIN_1000", envir = envir, inherits = TRUE)) {
    ns_log("Definizioni del piano principale gia' presenti in ambiente: non ricarico.")
    return(invisible(TRUE))
  }

  ns_log("Carico le definizioni da: %s", plan_script)
  txt <- readLines(plan_script, warn = FALSE)

  idx <- grep("^##[[:space:]]*11\\.[[:space:]]*AUTORUN SWITCH", txt)
  if (!length(idx)) idx <- grep("^print_plan_sizes\\(\\)", txt)
  if (!length(idx)) {
    stop("Non trovo il marcatore di AUTORUN SWITCH nello script originale. ",
         "Mi rifiuto di sourceare l'intero file perche' potrebbe lanciare simulazioni.")
  }

  keep <- seq_len(idx[1] - 1L)
  code <- paste(txt[keep], collapse = "\n")
  eval(parse(text = code), envir = envir)
  invisible(TRUE)
}

ns_load_plan_definitions()

## Carica il wrapper ETAS-FLP in modalita' "none" (solo definizioni, nessuna
## scansione/fit sui vecchi output). Passa dalla stessa dispatch logic del
## file originale, che a MODE == "none" si ferma dopo aver caricato tutto.
ns_load_flp_definitions <- function(flp_script = NONSTAT_FLP_SCRIPT) {
  if (exists("fit_etas_flp_full", mode = "function", inherits = TRUE)) {
    ns_log("fit_etas_flp_full() gia' presente in ambiente: non ricarico lo script FLP.")
    return(invisible(TRUE))
  }
  old_mode <- Sys.getenv("ETAS_FLP_MODE", unset = NA_character_)
  old_out_root <- Sys.getenv("ETAS_OUT_ROOT", unset = NA_character_)
  ## Puntiamo ETAS_OUT_ROOT a una cartella qualunque esistente (getwd()) solo
  ## per soddisfare i controlli interni dello script FLP: in modalita' "none"
  ## non viene comunque letto o scritto nulla li' dentro.
  Sys.setenv(ETAS_OUT_ROOT = NONSTAT_BASE_DIR)
  Sys.setenv(ETAS_FLP_MODE = "none")
  Sys.setenv(ETAS_PLAN_SCRIPT = NONSTAT_PLAN_SCRIPT)
  ns_log("Carico le definizioni ETAS-FLP da: %s (modalita' 'none')", flp_script)
  source(flp_script, local = .GlobalEnv, echo = FALSE)
  if (!is.na(old_mode)) Sys.setenv(ETAS_FLP_MODE = old_mode) else Sys.unsetenv("ETAS_FLP_MODE")
  if (!is.na(old_out_root)) Sys.setenv(ETAS_OUT_ROOT = old_out_root) else Sys.unsetenv("ETAS_OUT_ROOT")
  invisible(TRUE)
}

ns_load_flp_definitions()

required_objs <- c(
  "etasclass", "etasclass.par", "fit_etas_classic_full", "fit_etas_parametric_full",
  "fit_etas_flp_full", "background_fit_spec", "make_grid_eval",
  "fit_parameter_comparison", "fit_trigger_components",
  "fit_classification_metrics", "fit_background_metrics",
  "get_cat_sim", "get_sim_core", "make_eqcat_from_sim",
  "background_features_from_xy", "bg_cov_Z1_truth_full", "exact_covariate_functions_full",
  "safe_run", "saveRDS_atomic", "COUNT_MAIN_1000", "count_lookup",
  "TRUE_BASE_FULL", "TRUE_BETACOV_FULL", "TRUE_BETA_Z1_FULL", "M0_FULL", "B_FULL",
  "LONG_RANGE", "LAT_RANGE", "TMIN_USE", "T_LAG_USE", "MAX_EVENTS",
  "NDECLUST_FIT_NEW", "ITERLIM_FIT_NEW", "MULT_BG_FIT_NEW", "GRID_N",
  "NREP_CALIB", "N_CALIB_ITERS", "N_CALIB_BRACKET", "TOL_BG_REL", "TOL_TRIG_REL"
)
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
if (length(missing_objs)) {
  stop("Mancano questi oggetti dopo il caricamento delle definizioni: ",
       paste(missing_objs, collapse = ", "),
       ". Verifica i percorsi NONSTAT_PLAN_SCRIPT / NONSTAT_FLP_SCRIPT.")
}

## ------------------------------------------------------------------
## Le funzioni seguenti (estimate_classic_background_grid,
## estimate_parametric_background_grid, true_background_grid_df, e i
## loro helper) sono definite nello script originale DOPO il marcatore
## di AUTORUN SWITCH, insieme a codice che si esegue direttamente al
## caricamento (selezione dei file rappresentativi, costruzione di un
## grafico). Il loader sicuro qui sopra si ferma PRIMA di quel blocco
## apposta, per non rischiare di eseguire codice che dipende da oggetti
## (es. class_metrics_tbl) non ancora presenti in questa sessione.
## Le ridefiniamo qui in forma identica, cosi' da poterle riusare senza
## toccare o sourceare oltre il marcatore nello script originale.
## ------------------------------------------------------------------

if (!exists("weighted_kde_grid", mode = "function", inherits = TRUE)) {
  weighted_kde_grid <- function(x, y, weights, grid, h = NULL, chunk_size = 250L) {
    ok <- is.finite(x) & is.finite(y) & is.finite(weights) & weights >= 0
    x <- x[ok]; y <- y[ok]; weights <- weights[ok]
    if (!length(x)) return(rep(NA_real_, nrow(grid)))
    if (is.null(h) || length(h) < 2 || any(!is.finite(h)) || any(h <= 0)) {
      hx <- stats::bw.nrd0(x); hy <- stats::bw.nrd0(y); h <- c(hx, hy)
    }
    hx <- as.numeric(h[1]); hy <- as.numeric(h[2])
    if (!is.finite(hx) || hx <= 0) hx <- diff(range(x, na.rm = TRUE)) / 20
    if (!is.finite(hy) || hy <= 0) hy <- diff(range(y, na.rm = TRUE)) / 20
    gx <- as.numeric(grid$x); gy <- as.numeric(grid$y)
    dens <- numeric(length(gx))
    idx <- split(seq_along(x), ceiling(seq_along(x) / chunk_size))
    for (ii in idx) {
      dx <- outer(gx, x[ii], function(a, b) stats::dnorm((a - b) / hx) / hx)
      dy <- outer(gy, y[ii], function(a, b) stats::dnorm((a - b) / hy) / hy)
      dens <- dens + rowSums(dx * dy * matrix(weights[ii], nrow = length(gx), ncol = length(ii), byrow = TRUE))
    }
    dens
  }
  assign("weighted_kde_grid", weighted_kde_grid, envir = .GlobalEnv)
}

if (!exists("normalize_density_on_grid", mode = "function", inherits = TRUE)) {
  normalize_density_on_grid <- function(dens, grid) {
    dens <- as.numeric(dens); w <- grid$w
    dens[!is.finite(dens)] <- NA_real_
    integ <- sum(w * dens, na.rm = TRUE)
    if (!is.finite(integ) || integ <= 0) return(rep(NA_real_, length(dens)))
    dens / integ
  }
  assign("normalize_density_on_grid", normalize_density_on_grid, envir = .GlobalEnv)
}

if (!exists("center_log_shape", mode = "function", inherits = TRUE)) {
  center_log_shape <- function(dens, grid, eps = 1e-300) {
    logd <- log(pmax(dens, eps))
    logd - weighted.mean(logd, w = grid$w, na.rm = TRUE)
  }
  assign("center_log_shape", center_log_shape, envir = .GlobalEnv)
}

if (!exists("estimate_classic_background_grid", mode = "function", inherits = TRUE)) {
  estimate_classic_background_grid <- function(fit, grid) {
    if (is.null(fit) || is.null(fit$cat)) return(NULL)
    df <- as.data.frame(fit$cat)
    x <- if ("xcat.work" %in% names(df)) df$xcat.work else if ("x_km" %in% names(df)) df$x_km else df$long
    y <- if ("ycat.work" %in% names(df)) df$ycat.work else if ("y_km" %in% names(df)) df$y_km else df$lat
    rho <- fit$rho.weights
    if (is.null(rho) || length(rho) != nrow(df)) rho <- rep(1, nrow(df))
    h <- fit$hdef
    dens_raw <- weighted_kde_grid(x = x, y = y, weights = rho, grid = grid, h = h)
    dens <- normalize_density_on_grid(dens_raw, grid)
    data.frame(x = grid$x, y = grid$y, dens = dens,
              log_shape = center_log_shape(dens, grid), stringsAsFactors = FALSE)
  }
  assign("estimate_classic_background_grid", estimate_classic_background_grid, envir = .GlobalEnv)
}

if (!exists("estimate_parametric_background_grid", mode = "function", inherits = TRUE)) {
  estimate_parametric_background_grid <- function(fit, grid, bg_case) {
    if (is.null(fit) || is.null(fit$model.bg) || is.null(fit$model.bg$mod_global)) return(NULL)
    mod <- fit$model.bg$mod_global
    if (bg_case == "mark_cat") {
      grid_A <- grid; grid_B <- grid
      grid_A$mark_cat <- factor("A", levels = c("A", "B"))
      grid_B$mark_cat <- factor("B", levels = c("A", "B"))
      eta_A <- tryCatch(as.numeric(predict(mod, newdata = grid_A, type = "link")), error = function(e) rep(NA_real_, nrow(grid)))
      eta_B <- tryCatch(as.numeric(predict(mod, newdata = grid_B, type = "link")), error = function(e) rep(NA_real_, nrow(grid)))
      lambda_shape <- exp(eta_A) + exp(eta_B)
    } else {
      eta <- tryCatch(as.numeric(predict(mod, newdata = grid, type = "link")), error = function(e) rep(NA_real_, nrow(grid)))
      lambda_shape <- exp(eta)
    }
    dens <- normalize_density_on_grid(lambda_shape, grid)
    data.frame(x = grid$x, y = grid$y, dens = dens,
              log_shape = center_log_shape(dens, grid), stringsAsFactors = FALSE)
  }
  assign("estimate_parametric_background_grid", estimate_parametric_background_grid, envir = .GlobalEnv)
}

if (!exists("true_background_grid_df", mode = "function", inherits = TRUE)) {
  true_background_grid_df <- function(grid) {
    dens <- grid$dens_bg_true_spatial
    data.frame(x = grid$x, y = grid$y, dens = dens,
              log_shape = grid$eta_bg_true_spatial - weighted.mean(grid$eta_bg_true_spatial, w = grid$w, na.rm = TRUE),
              stringsAsFactors = FALSE)
  }
  assign("true_background_grid_df", true_background_grid_df, envir = .GlobalEnv)
}

grid_helper_objs <- c("estimate_classic_background_grid", "estimate_parametric_background_grid", "true_background_grid_df")
missing_grid_helpers <- grid_helper_objs[!vapply(grid_helper_objs, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
if (length(missing_grid_helpers)) {
  stop("Mancano ancora questi helper di griglia dopo la ridefinizione locale: ", paste(missing_grid_helpers, collapse = ", "))
}

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Il pacchetto 'ggplot2' e' necessario per i grafici.")
  library(ggplot2)
})

ns_log("Tutte le definizioni necessarie sono state caricate correttamente.")


## ================================================================
## 2. CONFIGURAZIONE DEL PIANO
## ================================================================

NONSTAT_OUT_ROOT <- ns_env_chr("NONSTAT_OUT_ROOT", file.path(NONSTAT_BASE_DIR, "etas_nonstationary_outputs_v1"))
NONSTAT_KAPPA    <- ns_env_num("NONSTAT_KAPPA", 1.6)   # rapporto fine/inizio = (1+k/2)/(1-k/2) = 9x per k=1.6
NONSTAT_NREP     <- ns_env_int("NONSTAT_NREP", 30L)
NONSTAT_NBINS    <- ns_env_int("NONSTAT_NBINS", 10L)
NONSTAT_MODE     <- tolower(ns_env_chr("NONSTAT_MODE", "none"))
NONSTAT_SEED_BASE <- ns_env_int("NONSTAT_SEED_BASE", 87000000L)

NONSTAT_BG_CASES <- c("constant", "cov1")   # scenari spaziali abbinati al trend temporale
NONSTAT_COUNT_CASE <- "N1000_balanced_50_50"

ns_dir_create(NONSTAT_OUT_ROOT)
ns_dir_create(file.path(NONSTAT_OUT_ROOT, "fits"))
ns_dir_create(file.path(NONSTAT_OUT_ROOT, "logs"))
ns_dir_create(file.path(NONSTAT_OUT_ROOT, "figures"))

ns_log("Configurazione: kappa=%.2f | nrep=%d | nbins=%d | bg_cases=%s | count_case=%s",
       NONSTAT_KAPPA, NONSTAT_NREP, NONSTAT_NBINS,
       paste(NONSTAT_BG_CASES, collapse = ","), NONSTAT_COUNT_CASE)


## ================================================================
## 3. FORMA TEMPORALE mu(t) = mu0 * g(t)  (trend lineare centrato)
## ================================================================

## g(t) ha media ESATTAMENTE 1 su [tmin,tmax] per costruzione, qualunque
## sia kappa (l'integrale del termine lineare centrato su un intervallo
## simmetrico e' nullo). Percio' mu0 mantiene lo stesso ruolo di sempre
## nella calibrazione (N_bg atteso = mu0 * tau).
mu_t_shape_linear <- function(t, tmin, tmax, kappa) {
  tau  <- tmax - tmin
  tmid <- tmin + tau / 2
  1 + kappa * (t - tmid) / tau
}

## Valore massimo di g(t) su [tmin,tmax], usato per il rifiuto-accettazione.
mu_t_shape_linear_max <- function(kappa) 1 + abs(kappa) / 2

if (NONSTAT_KAPPA >= 2 || NONSTAT_KAPPA <= -2) {
  stop("NONSTAT_KAPPA deve essere in (-2, 2) per garantire g(t) > 0 su tutta la finestra.")
}


## ================================================================
## 4. SIMULATORE CON BACKGROUND NON-STAZIONARIO NEL TEMPO
## ================================================================
## Reimplementa la generazione degli eventi di background (tempo NON
## uniforme, campionato da g(t) via rifiuto-accettazione; spazio secondo
## bg_case) e riusa la stessa logica di branching/triggering gia'
## presente in etas.par.sim_v4 (kernel Omori-Utsu temporale, kernel
## isotropo spaziale, stessa parametrizzazione). Non richiede alcuna
## modifica al file originale.

etas_sim_timevar <- function(
    params,                 # c(mu=mu0, k0, c, p, gamma, d, q)
    m0, b,
    tmin, t.lag,
    long.range, lat.range, longlat.to.km = TRUE,
    bg_case = c("constant", "cov1"),
    kappa,                  # forza del trend lineare in mu(t)
    trig_lp_coefs = c(m_rel = 0.7),
    n_support_obs = 1500, mult_support = 4,
    max_events = 50000,
    seed = NULL,
    return_bg_info = TRUE, return_support = TRUE
) {
  bg_case <- match.arg(bg_case)
  if (!is.null(seed)) set.seed(seed)

  mu    <- unname(params["mu"])
  k0    <- unname(params["k0"])
  c_    <- unname(params["c"])
  p     <- unname(params["p"])
  gamma <- unname(params["gamma"])
  d     <- unname(params["d"])
  q     <- unname(params["q"])

  tmax <- tmin + t.lag
  unit <- if (isTRUE(longlat.to.km)) 6371.3 * pi / 180 else 1
  xmin_km <- long.range[1] * unit; xmax_km <- long.range[2] * unit
  ymin_km <- lat.range[1]  * unit; ymax_km <- lat.range[2]  * unit
  area_km2 <- (xmax_km - xmin_km) * (ymax_km - ymin_km)
  window_km <- c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km)

  ## ---- supporto spaziale (per bg_info$support, usato dalle metriche
  ##      di forma gia' esistenti, es. compare_background_shape_param) ----
  n_support_dummy <- as.integer(mult_support * n_support_obs)
  n_support_tot   <- n_support_obs + n_support_dummy
  x_sup <- runif(n_support_tot, xmin_km, xmax_km)
  y_sup <- runif(n_support_tot, ymin_km, ymax_km)
  w_sup <- rep(area_km2 / n_support_tot, n_support_tot)

  eta_bg_fun <- function(x, y) {
    if (bg_case == "constant") return(rep(0, length(x)))
    feat <- background_features_from_xy(x, y, window_km)
    TRUE_BETA_Z1_FULL * bg_cov_Z1_truth_full(feat)
  }

  eta_bg_sup <- eta_bg_fun(x_sup, y_sup)
  lambda_shape_sup <- exp(eta_bg_sup)
  dens_bg_sup <- lambda_shape_sup / sum(w_sup * lambda_shape_sup)

  ## ---- numero di eventi di background: Poisson(mu0 * tau), IDENTICO
  ##      al caso stazionario perche' mean(g) = 1 per costruzione ----
  n0 <- stats::rpois(1, mu * t.lag)

  empty_out <- function() {
    list(
      cat.pois = data.frame(), cat.sim = data.frame(),
      n0 = 0L, nson = 0L, exploded = FALSE,
      truth = list(
        params = params, m0 = m0, b = b,
        bg = list(type = bg_case, coefs = if (bg_case == "cov1") c(Z1 = TRUE_BETA_Z1_FULL) else NULL),
        mu_t = list(type = "linear", kappa = kappa, tmin = tmin, tmax = tmax),
        trig = list(type = "linear", coefs = trig_lp_coefs),
        window_km = window_km, t_window = c(tmin = tmin, tmax = tmax)
      ),
      bg_info = if (isTRUE(return_bg_info)) list(
        window_km = window_km,
        support = if (isTRUE(return_support)) data.frame(x = x_sup, y = y_sup, w = w_sup,
                                                           eta_bg_true = eta_bg_sup,
                                                           lambda_shape = lambda_shape_sup,
                                                           dens_bg = dens_bg_sup) else NULL
      ) else NULL
    )
  }
  if (n0 == 0L) return(empty_out())

  ## ---- tempi di background: rifiuto-accettazione da g(t) ----
  gmax <- mu_t_shape_linear_max(kappa)
  sample_bg_times <- function(n_need) {
    out <- numeric(0)
    batch <- max(1000L, ceiling(3 * n_need))
    guard <- 0L
    while (length(out) < n_need) {
      guard <- guard + 1L
      if (guard > 10000L) stop("Rifiuto-accettazione per i tempi di background non converge.")
      tcand <- runif(batch, tmin, tmax)
      gcand <- mu_t_shape_linear(tcand, tmin, tmax, kappa)
      acc <- runif(batch) <= (gcand / gmax)
      if (any(acc)) out <- c(out, tcand[acc])
      if (length(out) < n_need && guard %% 20L == 0L) batch <- min(batch * 2L, 200000L)
    }
    out[seq_len(n_need)]
  }
  t0 <- sample_bg_times(n0)

  ## ---- posizioni di background: rifiuto-accettazione da f0(s) ----
  sample_bg_xy <- function(n_need) {
    if (bg_case == "constant") {
      return(data.frame(x = runif(n_need, xmin_km, xmax_km), y = runif(n_need, ymin_km, ymax_km)))
    }
    eta_max <- max(eta_bg_sup, na.rm = TRUE)
    out_x <- numeric(0); out_y <- numeric(0)
    batch <- max(1000L, ceiling(3 * n_need))
    guard <- 0L
    while (length(out_x) < n_need) {
      guard <- guard + 1L
      if (guard > 10000L) stop("Rifiuto-accettazione per le posizioni di background non converge.")
      xcand <- runif(batch, xmin_km, xmax_km)
      ycand <- runif(batch, ymin_km, ymax_km)
      etac <- eta_bg_fun(xcand, ycand)
      acc <- runif(batch) <= exp(etac - eta_max)
      if (any(acc)) { out_x <- c(out_x, xcand[acc]); out_y <- c(out_y, ycand[acc]) }
      if (length(out_x) < n_need && guard %% 20L == 0L) batch <- min(batch * 2L, 200000L)
    }
    data.frame(x = out_x[seq_len(n_need)], y = out_y[seq_len(n_need)])
  }
  xy0 <- sample_bg_xy(n0)
  x0 <- xy0$x; y0 <- xy0$y

  beta_GR <- log(10) * b
  m0_bg <- m0 + rexp(n0, rate = beta_GR)
  mrel_bg <- m0_bg - m0

  feat0 <- background_features_from_xy(x0, y0, window_km)
  eta_bg0 <- eta_bg_fun(x0, y0)
  eta_trig0 <- as.numeric(trig_lp_coefs["m_rel"]) * mrel_bg

  long0 <- if (isTRUE(longlat.to.km)) x0 / unit else x0
  lat0  <- if (isTRUE(longlat.to.km)) y0 / unit else y0

  cat.pois <- data.frame(
    event_id = seq_len(n0), father_id = 0L, lgen = 0L,
    time = t0, lat = lat0, long = long0, z = 0, magn1 = m0_bg,
    x_km = x0, y_km = y0, m_rel = mrel_bg,
    x_true = x0, y_true = y0,
    x_std_true = feat0$x_std, y_std_true = feat0$y_std, r2_center_true = feat0$r2_center,
    eta_bg_true = eta_bg0, eta_trig_true = eta_trig0,
    stringsAsFactors = FALSE
  )
  if (bg_case == "cov1") cat.pois$Z1 <- bg_cov_Z1_truth_full(feat0)

  cat.pois <- cat.pois[order(cat.pois$time), , drop = FALSE]
  rownames(cat.pois) <- NULL
  cat.pois$event_id <- seq_len(nrow(cat.pois))
  cat.new <- cat.pois

  ## ---- branching / triggering: identico, per costruzione, allo schema
  ##      di etas.par.sim_v4 (il meccanismo di triggering non dipende
  ##      da come sono distribuiti nel tempo i genitori di background) ----
  ak <- k0 * c_^(1 - p) / (p - 1)
  sk <- (pi * d^(1 - q)) / (q - 1)

  sample_xy_etas <- function(n, x0c, y0c, d, q) {
    theta <- runif(n, 0, 2 * pi)
    U <- runif(n)
    R <- sqrt(d * (U^(1 / (1 - q)) - 1))
    cbind(x0c + R * cos(theta), y0c + R * sin(theta))
  }

  i <- 0L
  exploded <- FALSE
  while (i < nrow(cat.new)) {
    i <- i + 1L
    m_rel_i <- cat.new$m_rel[i]
    eta_trig_i <- as.numeric(trig_lp_coefs["m_rel"]) * m_rel_i
    cat.new$eta_trig_true[i] <- eta_trig_i

    n_exp_i <- ak * sk * exp(gamma * m_rel_i + eta_trig_i)
    if (!is.finite(n_exp_i) || n_exp_i < 0) n_exp_i <- 0
    ni <- stats::rpois(1, n_exp_i)

    if (ni > 0L) {
      t_child <- c_ * runif(ni)^(-1 / (p - 1)) - c_ + cat.new$time[i]
      xy_child <- sample_xy_etas(ni, cat.new$x_km[i], cat.new$y_km[i], d, q)
      inside <- (t_child > tmin) & (t_child < tmax) &
        (xy_child[, 1] > xmin_km) & (xy_child[, 1] < xmax_km) &
        (xy_child[, 2] > ymin_km) & (xy_child[, 2] < ymax_km)

      if (any(inside)) {
        nt <- sum(inside)
        x1 <- xy_child[inside, 1]; y1 <- xy_child[inside, 2]; t1 <- t_child[inside]
        m1 <- m0 + rexp(nt, rate = beta_GR)
        mrel1 <- m1 - m0

        feat1 <- background_features_from_xy(x1, y1, window_km)
        eta_bg1 <- eta_bg_fun(x1, y1)
        eta_trig1 <- as.numeric(trig_lp_coefs["m_rel"]) * mrel1

        long1 <- if (isTRUE(longlat.to.km)) x1 / unit else x1
        lat1  <- if (isTRUE(longlat.to.km)) y1 / unit else y1

        child_df <- data.frame(
          event_id = NA_integer_, father_id = as.integer(cat.new$event_id[i]),
          lgen = as.integer(cat.new$lgen[i] + 1L),
          time = t1, lat = lat1, long = long1, z = 0, magn1 = m1,
          x_km = x1, y_km = y1, m_rel = mrel1,
          x_true = x1, y_true = y1,
          x_std_true = feat1$x_std, y_std_true = feat1$y_std, r2_center_true = feat1$r2_center,
          eta_bg_true = eta_bg1, eta_trig_true = eta_trig1,
          stringsAsFactors = FALSE
        )
        if (bg_case == "cov1") child_df$Z1 <- bg_cov_Z1_truth_full(feat1)

        cat.new <- rbind(cat.new, child_df)
      }
    }

    if (nrow(cat.new) > max_events) { exploded <- TRUE; break }
  }

  cat.new <- cat.new[order(cat.new$time), , drop = FALSE]
  rownames(cat.new) <- NULL
  cat.new$event_id <- seq_len(nrow(cat.new))

  list(
    cat.pois = cat.pois, cat.sim = cat.new,
    n0 = as.integer(n0), nson = as.integer(nrow(cat.new) - n0), exploded = exploded,
    truth = list(
      params = params, m0 = m0, b = b,
      bg = list(type = bg_case, coefs = if (bg_case == "cov1") c(Z1 = TRUE_BETA_Z1_FULL) else NULL),
      mu_t = list(type = "linear", kappa = kappa, tmin = tmin, tmax = tmax),
      trig = list(type = "linear", coefs = trig_lp_coefs),
      window_km = window_km, t_window = c(tmin = tmin, tmax = tmax)
    ),
    bg_info = if (isTRUE(return_bg_info)) list(
      window_km = window_km,
      support = if (isTRUE(return_support)) data.frame(x = x_sup, y = y_sup, w = w_sup,
                                                         eta_bg_true = eta_bg_sup,
                                                         lambda_shape = lambda_shape_sup,
                                                         dens_bg = dens_bg_sup) else NULL
    ) else NULL
  )
}


## ================================================================
## 5. WRAPPER DI SIMULAZIONE E CALIBRAZIONE (riusa la logica esistente,
##    sostituendo solo la chiamata al simulatore)
## ================================================================

simulate_one_nonstat <- function(bg_case, kappa, count_regime, seed, mu_value, k0_value,
                                 return_bg_info = TRUE, return_support = TRUE) {
  params_now <- TRUE_BASE_FULL
  params_now["mu"] <- mu_value
  params_now["k0"] <- k0_value

  etas_sim_timevar(
    params = params_now, m0 = M0_FULL, b = B_FULL,
    tmin = TMIN_USE, t.lag = count_regime$T_lag,
    long.range = LONG_RANGE, lat.range = LAT_RANGE, longlat.to.km = TRUE,
    bg_case = bg_case, kappa = kappa,
    trig_lp_coefs = c(m_rel = TRUE_BETACOV_FULL),
    n_support_obs = 1500, mult_support = 4,
    max_events = MAX_EVENTS, seed = seed,
    return_bg_info = return_bg_info, return_support = return_support
  )
}

## Calibrazione: mu si fissa analiticamente a target_bg / T_lag (identico
## al caso stazionario, valido perche' mean(g)=1); solo k0 va calibrato
## con bracket + bisezione sul numero medio di eventi triggered, con la
## STESSA logica di calibrate_cell_poster ma chiamando il nuovo simulatore.
calibrate_cell_nonstat <- function(bg_case, kappa, count_regime,
                                   nrep_pilot = NREP_CALIB,
                                   n_iters = N_CALIB_ITERS,
                                   n_bracket = N_CALIB_BRACKET,
                                   seed_base = 77000000L,
                                   verbose = TRUE,
                                   max_exploded_rate = 0.20) {
  target_bg <- count_regime$N_bg_target
  target_trig <- count_regime$N_trig_target
  mu_cur <- max(target_bg / count_regime$T_lag, 1e-8)

  eval_id <- 0L
  logs <- list()

  eval_k0 <- function(k0_value, phase, iter) {
    eval_id <<- eval_id + 1L
    n0_vec <- numeric(nrep_pilot); nson_vec <- numeric(nrep_pilot); exploded_vec <- logical(nrep_pilot)
    for (r in seq_len(nrep_pilot)) {
      seed_r <- seed_base + 100000L * eval_id + r
      sim_r <- try(
        simulate_one_nonstat(bg_case, kappa, count_regime, seed_r, mu_cur, k0_value,
                             return_bg_info = FALSE, return_support = FALSE),
        silent = TRUE
      )
      if (inherits(sim_r, "try-error")) {
        n0_vec[r] <- NA_real_; nson_vec[r] <- NA_real_; exploded_vec[r] <- TRUE
      } else {
        n0_vec[r] <- sim_r$n0; nson_vec[r] <- sim_r$nson; exploded_vec[r] <- isTRUE(sim_r$exploded)
      }
    }
    mean_n0 <- mean(n0_vec, na.rm = TRUE); mean_nson <- mean(nson_vec, na.rm = TRUE)
    exploded_rate <- mean(exploded_vec, na.rm = TRUE)
    row <- data.frame(
      eval_id = eval_id, phase = phase, iter = iter, bg_case = bg_case, kappa = kappa,
      mu = mu_cur, k0 = k0_value, mean_n0 = mean_n0, mean_nson = mean_nson,
      target_bg = target_bg, target_trig = target_trig,
      rel_err_bg = (mean_n0 - target_bg) / target_bg,
      rel_err_trig = (mean_nson - target_trig) / target_trig,
      exploded_rate = exploded_rate
    )
    logs[[length(logs) + 1L]] <<- row
    if (isTRUE(verbose)) {
      cat(sprintf("[CALIB-NONSTAT:%s] %-10s | iter %02d | mu=%.5f k0=%.6f | n0 %.1f/%.1f | nson %.1f/%.1f | expl %.2f\n",
                  phase, bg_case, iter, mu_cur, k0_value, mean_n0, target_bg, mean_nson, target_trig, exploded_rate))
    }
    row
  }

  k0_guess <- max(count_regime$k0_init, 1e-8)
  low_k0 <- 0; high_k0 <- k0_guess
  high_eval <- eval_k0(high_k0, "bracket", 1L)
  b_iter <- 1L
  while (is.finite(high_eval$mean_nson) && high_eval$mean_nson < target_trig &&
         high_eval$exploded_rate <= max_exploded_rate && b_iter < n_bracket) {
    low_k0 <- high_k0; high_k0 <- high_k0 * 1.45
    b_iter <- b_iter + 1L
    high_eval <- eval_k0(high_k0, "bracket", b_iter)
  }

  for (it in seq_len(n_iters)) {
    mid_k0 <- (low_k0 + high_k0) / 2
    mid_eval <- eval_k0(mid_k0, "bisect", it)
    if (!is.finite(mid_eval$mean_nson) || mid_eval$exploded_rate > max_exploded_rate) {
      high_k0 <- mid_k0
    } else if (mid_eval$mean_nson < target_trig) {
      low_k0 <- mid_k0
    } else {
      high_k0 <- mid_k0
    }
    if (is.finite(mid_eval$rel_err_trig) && abs(mid_eval$rel_err_trig) < TOL_TRIG_REL &&
        is.finite(mid_eval$rel_err_bg) && abs(mid_eval$rel_err_bg) < TOL_BG_REL) break
  }

  pilot_log <- do.call(rbind, logs)
  cand <- pilot_log[is.finite(pilot_log$mean_nson) & pilot_log$exploded_rate <= max_exploded_rate, , drop = FALSE]
  if (!nrow(cand)) cand <- pilot_log[is.finite(pilot_log$mean_nson), , drop = FALSE]
  if (!nrow(cand)) stop("Calibrazione fallita per bg_case=", bg_case)
  cand$abs_rel_err_trig <- abs((cand$mean_nson - target_trig) / target_trig)
  cand$abs_rel_err_bg <- abs((cand$mean_n0 - target_bg) / target_bg)
  cand$score <- cand$abs_rel_err_trig + 0.25 * cand$abs_rel_err_bg + 2 * cand$exploded_rate
  best <- cand[which.min(cand$score), , drop = FALSE]

  if (isTRUE(verbose)) {
    cat(sprintf("[CALIB-NONSTAT:SELECT] %-10s | mu=%.5f k0=%.6f | n0 %.1f/%.1f | nson %.1f/%.1f\n",
                bg_case, best$mu, best$k0, best$mean_n0, target_bg, best$mean_nson, target_trig))
  }
  list(mu_cal = as.numeric(best$mu), k0_cal = as.numeric(best$k0), pilot_log = pilot_log)
}


## ================================================================
## 6. METRICHE NUOVE: rMISE spaziale su griglia comune (3 modelli) e
##    diagnostica di recupero temporale
## ================================================================

## rMISE, stessa formula usata nel paper (Sez. 4.3): errore quadratico
## integrato relativo, pesato dai pesi di griglia w.
rmise_bg_grid <- function(dens_hat, dens_true, w) {
  num <- sum(w * (dens_hat - dens_true)^2, na.rm = TRUE)
  den <- sum(w * dens_true^2, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

## Calcola rMISEbg per i tre modelli sulla STESSA griglia di valutazione
## spaziale (riusa make_grid_eval/estimate_*_grid gia' esistenti). Questo
## e' il numero chiave per rispondere a: "il drift temporale finisce nel
## background spaziale?" -- ci aspettiamo NO, perche' la forma marginale
## vera resta f0(s) per costruzione (separabilita' di mu(t) e f0(s)).
compute_rmise_bg_three_models <- function(fit_classic, fit_flp, fit_param, bg_case, grid_n = GRID_N) {
  grid <- make_grid_eval(bg_case, grid_n = grid_n)
  true_df <- true_background_grid_df(grid)

  out <- data.frame(model = c("classic", "flp", "param"), rmise_bg = NA_real_)

  if (!is.null(fit_classic)) {
    df_c <- tryCatch(estimate_classic_background_grid(fit_classic, grid), error = function(e) NULL)
    if (!is.null(df_c)) out$rmise_bg[out$model == "classic"] <- rmise_bg_grid(df_c$dens, true_df$dens, grid$w)
  }
  if (!is.null(fit_flp)) {
    df_f <- tryCatch(estimate_classic_background_grid(fit_flp, grid), error = function(e) NULL)
    if (!is.null(df_f)) out$rmise_bg[out$model == "flp"] <- rmise_bg_grid(df_f$dens, true_df$dens, grid$w)
  }
  if (!is.null(fit_param)) {
    df_p <- tryCatch(estimate_parametric_background_grid(fit_param, grid, bg_case = bg_case), error = function(e) NULL)
    if (!is.null(df_p)) out$rmise_bg[out$model == "param"] <- rmise_bg_grid(df_p$dens, true_df$dens, grid$w)
  }
  out
}

## Diagnostica temporale: raggruppa gli eventi in NONSTAT_NBINS bin
## temporali equispaziati su [tmin,tmax] e confronta, per ciascun
## modello:
##   - la massa di probabilita' posteriore di background stimata nel bin
##     (somma di rho_i sugli eventi del bin, normalizzata sul totale)
##   - la massa attesa VERA nel bin sotto mu(t) = mu0*g(t), cioe'
##     l'integrale di g(t) sul bin diviso l'integrale su tutta la finestra
## IMPORTANTE: fit_obj$rho.weights e' allineato a fit_obj$cat, che puo'
## avere un ordine di riga diverso da cat_sim_df (stesso motivo per cui
## fit_classification_metrics() nello script principale fa sempre il
## match tramite event_id anziche' assumere lo stesso ordine). Qui
## seguiamo la stessa convenzione: uniamo per event_id invece di
## assumere corrispondenza posizionale.
time_binned_background_recovery <- function(cat_sim_df, fit_obj, tmin, tmax, kappa, n_bins = NONSTAT_NBINS) {
  if (is.null(fit_obj) || is.null(fit_obj$rho.weights) || is.null(fit_obj$cat)) return(NULL)
  rho <- as.numeric(fit_obj$rho.weights)
  df_fit <- as.data.frame(fit_obj$cat)
  if (length(rho) != nrow(df_fit)) return(NULL)

  if (!("event_id" %in% names(df_fit)) || !("event_id" %in% names(cat_sim_df))) return(NULL)
  time_lookup <- cat_sim_df[, c("event_id", "time")]
  m <- match(df_fit$event_id, time_lookup$event_id)
  if (any(is.na(m))) return(NULL)
  times <- time_lookup$time[m]

  breaks <- seq(tmin, tmax, length.out = n_bins + 1L)
  bin_id <- cut(times, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  ## massa vera attesa per bin: integrale di g(t) sul bin / integrale totale
  bin_mid <- (breaks[-1] + breaks[-(n_bins + 1L)]) / 2
  bin_width <- diff(breaks)
  g_true_mid <- mu_t_shape_linear(bin_mid, tmin, tmax, kappa)
  true_mass_bin <- g_true_mid * bin_width
  true_mass_frac <- true_mass_bin / sum(true_mass_bin)

  fitted_mass_bin <- as.numeric(tapply(rho, factor(bin_id, levels = seq_len(n_bins)), sum, na.rm = TRUE))
  fitted_mass_bin[is.na(fitted_mass_bin)] <- 0
  fitted_mass_frac <- fitted_mass_bin / sum(fitted_mass_bin)

  n_events_bin <- as.numeric(table(factor(bin_id, levels = seq_len(n_bins))))

  data.frame(
    bin = seq_len(n_bins), t_mid = bin_mid,
    n_events = n_events_bin,
    true_mass_frac = true_mass_frac,
    fitted_mass_frac = fitted_mass_frac
  )
}


## ================================================================
## 6b. FIT ETAS-FLP CON RETRY (starting values perturbati)
## ================================================================
## Diagnosticato empiricamente: con trend temporale marcato (kappa alto),
## la selezione di bandwidth FLP puo' incontrare configurazioni locali
## degeneri nella coda della finestra dove mu(t) e' basso (pochi eventi di
## background attesi), producendo "NA/NaN/Inf in foreign function call"
## dentro nlm(). Il problema e' specifico dell'interazione trend+FLP (la
## selezione forward-predittiva della bandwidth scandisce gli eventi in
## ordine temporale) e non si presenta con kappa=0 sullo stesso simulatore,
## quindi non e' un errore di generazione del catalogo. Rimediamo con un
## retry a starting values leggermente perturbati: un punto di partenza
## diverso spesso evita la configurazione locale degenere senza alterare
## la natura del confronto (il catalogo resta lo stesso).

perturb_starts <- function(starts, sd_log = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- starts
  for (nm in c("mu", "k0", "c", "d")) {
    out[[nm]] <- starts[[nm]] * exp(rnorm(1, mean = 0, sd = sd_log))
  }
  out
}

fit_etas_flp_robust <- function(cat_sim_df, starts, ndeclust, iterlim,
                                max_attempts = 3L, sd_log = 0.15, seed_base = NULL) {
  attempt <- 0L
  last_result <- NULL
  while (attempt < max_attempts) {
    attempt <- attempt + 1L
    starts_try <- if (attempt == 1L) starts else perturb_starts(starts, sd_log = sd_log,
                                                                seed = (seed_base %||% 0L) + attempt)
    res <- safe_run(function() {
      fit_etas_flp_full(cat_sim_df, starts = starts_try, ndeclust = ndeclust, iterlim = iterlim)
    })
    res$n_attempts <- attempt
    res$starts_used <- starts_try
    last_result <- res
    if (isTRUE(res$ok)) {
      if (attempt > 1L) ns_log("  fit flp riuscito al tentativo %d/%d dopo perturbazione degli starting values", attempt, max_attempts)
      return(res)
    }
    ns_log("  fit flp tentativo %d/%d fallito: %s", attempt, max_attempts, res$error_message, .level = "WARN")
  }
  last_result
}


## ================================================================
## 7. UNA REPLICA COMPLETA: simula + fitta i tre modelli + metriche
## ================================================================

run_one_catalog_nonstat <- function(bg_case, kappa, count_regime, rep_id, mu_cal, k0_cal,
                                    out_root = NONSTAT_OUT_ROOT, seed_base = NONSTAT_SEED_BASE,
                                    force = FALSE) {
  scen_id <- paste0(bg_case, "_trend")
  file_out <- file.path(out_root, "fits", sprintf("%s__%s__rep%03d.rds", scen_id, count_regime$count_case, rep_id))

  if (!isTRUE(force) && file.exists(file_out)) {
    ns_log("Gia' presente, salto: %s", basename(file_out))
    return(invisible(file_out))
  }

  seed_rep <- seed_base + 1000L * rep_id + (if (bg_case == "cov1") 500L else 0L)
  ns_log("Simulo %s | rep=%03d | seed=%d", scen_id, rep_id, seed_rep)

  sim_result <- safe_run(function() {
    simulate_one_nonstat(bg_case, kappa, count_regime, seed_rep, mu_cal, k0_cal,
                         return_bg_info = TRUE, return_support = TRUE)
  })

  if (!isTRUE(sim_result$ok)) {
    ns_log("Simulazione fallita: %s", sim_result$error_message, .level = "ERROR")
    obj <- list(metadata = list(bg_case = scen_id, kappa = kappa, count_case = count_regime$count_case, rep = rep_id, seed = seed_rep),
               simulation = sim_result, fits = list(), metrics = list())
    ns_saveRDS_atomic(obj, file_out)
    return(invisible(file_out))
  }

  sim_obj <- sim_result$value
  cat_sim_df <- make_eqcat_from_sim(sim_obj)

  ns_log("  n_eventi=%d (n0=%d, nson=%d, exploded=%s)", nrow(cat_sim_df), sim_obj$n0, sim_obj$nson, sim_obj$exploded)

  starts <- list(mu = mu_cal, k0 = k0_cal, c = as.numeric(TRUE_BASE_FULL["c"]),
                 p = as.numeric(TRUE_BASE_FULL["p"]), gamma = 0,
                 d = as.numeric(TRUE_BASE_FULL["d"]), q = as.numeric(TRUE_BASE_FULL["q"]),
                 betacov = TRUE_BETACOV_FULL)

  fit_classic <- safe_run(function() fit_etas_classic_full(cat_sim_df, starts = starts))
  fit_flp     <- fit_etas_flp_robust(cat_sim_df, starts = starts,
                                     ndeclust = NDECLUST_FIT_NEW, iterlim = ITERLIM_FIT_NEW,
                                     max_attempts = 3L, sd_log = 0.15, seed_base = seed_rep)
  fit_param   <- safe_run(function() fit_etas_parametric_full(cat_sim_df, bg_case = bg_case,
                                                               formula_key = "correct", starts = starts))

  for (nm in c("classic", "flp", "param")) {
    fr <- get(paste0("fit_", nm))
    n_att <- fr$n_attempts %||% 1L
    if (isTRUE(fr$ok)) ns_log("  fit %-8s OK   (%.1f sec, tentativi=%d)", nm, fr$elapsed_sec, n_att)
    else ns_log("  fit %-8s FAIL dopo %d tentativi: %s", nm, n_att, fr$error_message, .level = "WARN")
  }

  fv_classic <- if (isTRUE(fit_classic$ok)) fit_classic$value else NULL
  fv_flp     <- if (isTRUE(fit_flp$ok))     fit_flp$value     else NULL
  fv_param   <- if (isTRUE(fit_param$ok))   fit_param$value   else NULL

  true_params_vec <- c(mu = as.numeric(sim_obj$truth$params["mu"]), k0 = as.numeric(sim_obj$truth$params["k0"]))

  event_truth <- data.frame(
    event_id = cat_sim_df$event_id %||% seq_len(nrow(cat_sim_df)),
    time = cat_sim_df$time,
    m_rel = cat_sim_df$m_rel %||% (cat_sim_df$magn1 - M0_FULL),
    is_background_true = cat_sim_df$father_id == 0
  )

  ## rMISE spaziale sulla stessa griglia per i tre modelli
  bg_shape_metrics <- tryCatch(
    compute_rmise_bg_three_models(fv_classic, fv_flp, fv_param, bg_case = bg_case),
    error = function(e) { ns_log("compute_rmise_bg_three_models fallita: %s", conditionMessage(e), .level = "WARN"); NULL }
  )

  ## bias sui parametri di triggering (mu, k0, c, p, d, q, betacov), stesse
  ## funzioni gia' usate nel piano principale
  trig_bias <- list(
    classic = tryCatch(fit_parameter_comparison(fv_classic, true_params_vec, cat_sim_df = event_truth), error = function(e) NULL),
    flp     = tryCatch(fit_parameter_comparison(fv_flp,     true_params_vec, cat_sim_df = event_truth), error = function(e) NULL),
    param   = tryCatch(fit_parameter_comparison(fv_param,   true_params_vec, cat_sim_df = event_truth), error = function(e) NULL)
  )

  ## bias sul coefficiente spaziale Z1 (solo per bg_case == "cov1", solo per ETAS-P)
  bg_coef_bias <- NULL
  if (bg_case == "cov1" && !is.null(fv_param)) {
    bg_coef_bias <- tryCatch(fit_background_metrics(fv_param, bg_case = "cov1",
                                                     grid_truth = make_grid_eval("cov1", grid_n = GRID_N))$coef_comparison,
                             error = function(e) NULL)
  }

  ## accuratezza di classificazione background/triggered, stessa funzione
  ## gia' usata nel piano principale
  classif <- list(
    classic = tryCatch(fit_classification_metrics(fv_classic, event_truth, threshold = 0.5), error = function(e) NULL),
    flp     = tryCatch(fit_classification_metrics(fv_flp,     event_truth, threshold = 0.5), error = function(e) NULL),
    param   = tryCatch(fit_classification_metrics(fv_param,   event_truth, threshold = 0.5), error = function(e) NULL)
  )

  ## diagnostica temporale: la parte NUOVA e piu' importante per rispondere
  ## al revisore
  time_binned <- list(
    classic = tryCatch(time_binned_background_recovery(cat_sim_df, fv_classic, TMIN_USE, TMIN_USE + count_regime$T_lag, kappa), error = function(e) NULL),
    flp     = tryCatch(time_binned_background_recovery(cat_sim_df, fv_flp,     TMIN_USE, TMIN_USE + count_regime$T_lag, kappa), error = function(e) NULL),
    param   = tryCatch(time_binned_background_recovery(cat_sim_df, fv_param,   TMIN_USE, TMIN_USE + count_regime$T_lag, kappa), error = function(e) NULL)
  )

  obj <- list(
    metadata = list(bg_case = scen_id, spatial_case = bg_case, kappa = kappa,
                    count_case = count_regime$count_case, rep = rep_id, seed = seed_rep,
                    n_attempts_classic = fit_classic$n_attempts %||% 1L,
                    n_attempts_flp = fit_flp$n_attempts %||% 1L,
                    n_attempts_param = fit_param$n_attempts %||% 1L),
    calibration = list(mu_cal = mu_cal, k0_cal = k0_cal),
    simulation = sim_result,
    fits = list(classic = fit_classic, flp = fit_flp, param = fit_param),
    metrics = list(
      bg_shape = bg_shape_metrics,
      trigger_bias = trig_bias,
      bg_coef_bias = bg_coef_bias,
      classification = classif,
      time_binned = time_binned
    )
  )

  ns_saveRDS_atomic(obj, file_out)
  ns_append_csv_row(
    data.frame(time = as.character(Sys.time()), bg_case = scen_id, rep = rep_id,
              n_events = nrow(cat_sim_df), n0 = sim_obj$n0, nson = sim_obj$nson,
              classic_ok = isTRUE(fit_classic$ok), flp_ok = isTRUE(fit_flp$ok), param_ok = isTRUE(fit_param$ok),
              flp_n_attempts = fit_flp$n_attempts %||% 1L,
              file = file_out, stringsAsFactors = FALSE),
    file.path(out_root, "logs", "progress_log_nonstat.csv")
  )

  invisible(file_out)
}


## ================================================================
## 8. ORCHESTRAZIONE DEL PIANO COMPLETO
## ================================================================

run_nonstat_plan <- function(bg_cases = NONSTAT_BG_CASES, kappa = NONSTAT_KAPPA,
                             count_case = NONSTAT_COUNT_CASE, nrep = NONSTAT_NREP,
                             out_root = NONSTAT_OUT_ROOT, force = FALSE) {
  count_regime <- count_lookup(count_case)

  calib_file <- file.path(out_root, "calibration_nonstat.rds")
  if (file.exists(calib_file) && !isTRUE(force)) {
    ns_log("Calibrazione gia' presente, la ricarico: %s", calib_file)
    calib_all <- readRDS(calib_file)
  } else {
    calib_all <- list()
    for (bg_case in bg_cases) {
      ns_log("Calibrazione per bg_case=%s, kappa=%.2f ...", bg_case, kappa)
      calib_all[[bg_case]] <- calibrate_cell_nonstat(bg_case, kappa, count_regime, verbose = TRUE)
    }
    ns_saveRDS_atomic(calib_all, calib_file)
  }

  for (bg_case in bg_cases) {
    mu_cal <- calib_all[[bg_case]]$mu_cal
    k0_cal <- calib_all[[bg_case]]$k0_cal
    ns_log("Piano: bg_case=%s | mu_cal=%.5f | k0_cal=%.6f | nrep=%d", bg_case, mu_cal, k0_cal, nrep)
    for (rep_id in seq_len(nrep)) {
      tryCatch(
        run_one_catalog_nonstat(bg_case, kappa, count_regime, rep_id, mu_cal, k0_cal, out_root = out_root, force = force),
        error = function(e) ns_log("Errore fatale rep=%d bg_case=%s: %s", rep_id, bg_case, conditionMessage(e), .level = "ERROR")
      )
      gc(verbose = FALSE)
    }
  }
  invisible(TRUE)
}


## ================================================================
## 9. RACCOLTA RISULTATI IN TABELLE LUNGHE
## ================================================================

collect_nonstat_results <- function(out_root = NONSTAT_OUT_ROOT) {
  files <- list.files(file.path(out_root, "fits"), pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) {
    ns_log("Nessun file trovato in %s", file.path(out_root, "fits"), .level = "WARN")
    return(list(bg_shape = data.frame(), trigger_bias = data.frame(),
               bg_coef_bias = data.frame(), classification = data.frame(), time_binned = data.frame()))
  }

  bg_shape_rows <- list(); trig_rows <- list(); coef_rows <- list()
  classif_rows <- list(); time_rows <- list(); fit_status_rows <- list()

  for (f in files) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj) || is.null(obj$metadata)) next
    meta <- obj$metadata

    fit_status_rows[[length(fit_status_rows) + 1L]] <- data.frame(
      bg_case = meta$bg_case, spatial_case = meta$spatial_case, rep = meta$rep,
      classic_ok = isTRUE(obj$fits$classic$ok),
      flp_ok = isTRUE(obj$fits$flp$ok),
      param_ok = isTRUE(obj$fits$param$ok),
      n_attempts_classic = meta$n_attempts_classic %||% 1L,
      n_attempts_flp = meta$n_attempts_flp %||% 1L,
      n_attempts_param = meta$n_attempts_param %||% 1L,
      stringsAsFactors = FALSE
    )

    if (!is.null(obj$metrics$bg_shape)) {
      d <- obj$metrics$bg_shape
      d$bg_case <- meta$bg_case; d$spatial_case <- meta$spatial_case; d$rep <- meta$rep
      bg_shape_rows[[length(bg_shape_rows) + 1L]] <- d
    }

    for (model_nm in c("classic", "flp", "param")) {
      tb <- obj$metrics$trigger_bias[[model_nm]]
      if (!is.null(tb)) {
        tb$model <- model_nm; tb$bg_case <- meta$bg_case; tb$spatial_case <- meta$spatial_case; tb$rep <- meta$rep
        trig_rows[[length(trig_rows) + 1L]] <- tb
      }
      cl <- obj$metrics$classification[[model_nm]]
      if (!is.null(cl)) {
        cl_df <- as.data.frame(as.list(cl), stringsAsFactors = FALSE)
        cl_df$model <- model_nm; cl_df$bg_case <- meta$bg_case; cl_df$rep <- meta$rep
        classif_rows[[length(classif_rows) + 1L]] <- cl_df
      }
      tbin <- obj$metrics$time_binned[[model_nm]]
      if (!is.null(tbin)) {
        tbin$model <- model_nm; tbin$bg_case <- meta$bg_case; tbin$spatial_case <- meta$spatial_case; tbin$rep <- meta$rep
        time_rows[[length(time_rows) + 1L]] <- tbin
      }
    }

    if (!is.null(obj$metrics$bg_coef_bias)) {
      cb <- obj$metrics$bg_coef_bias
      cb$bg_case <- meta$bg_case; cb$rep <- meta$rep
      coef_rows[[length(coef_rows) + 1L]] <- cb
    }
  }

  rb <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  list(
    bg_shape = rb(bg_shape_rows),
    trigger_bias = rb(trig_rows),
    bg_coef_bias = rb(coef_rows),
    classification = rb(classif_rows),
    time_binned = rb(time_rows),
    fit_status = rb(fit_status_rows)
  )
}

## Tasso di successo e numero medio di tentativi per modello e scenario.
## Da controllare SEMPRE prima di interpretare le tabelle di bias: un
## tasso di successo basso per un modello segnala un campione distorto
## (solo i cataloghi "facili" sopravvivono nel confronto).
summarise_fit_success_rate <- function(results = collect_nonstat_results()) {
  fs <- results$fit_status
  if (!nrow(fs)) { ns_log("Nessun dato di stato dei fit disponibile.", .level = "WARN"); return(invisible(NULL)) }

  out <- do.call(rbind, lapply(c("classic", "flp", "param"), function(nm) {
    ok_col <- fs[[paste0(nm, "_ok")]]
    att_col <- fs[[paste0("n_attempts_", nm)]]
    agg <- aggregate(ok_col ~ bg_case, data = data.frame(bg_case = fs$bg_case, ok_col = ok_col),
                     FUN = function(x) mean(x, na.rm = TRUE))
    agg_att <- aggregate(att_col ~ bg_case, data = data.frame(bg_case = fs$bg_case, att_col = att_col),
                         FUN = function(x) mean(x, na.rm = TRUE))
    data.frame(model = nm, bg_case = agg$bg_case,
              success_rate = agg$ok_col, mean_n_attempts = agg_att$att_col,
              n_reps = as.numeric(table(fs$bg_case)[agg$bg_case]),
              stringsAsFactors = FALSE)
  }))
  out
}


## ================================================================
## 10. TABELLE RIASSUNTIVE (stile Tabelle S1-S8 del supplementary)
## ================================================================

summarise_nonstat_results <- function(results = collect_nonstat_results()) {
  bg_summary <- NULL
  if (nrow(results$bg_shape)) {
    bg_summary <- aggregate(
      rmise_bg ~ bg_case + model, data = results$bg_shape,
      FUN = function(x) c(median = median(x, na.rm = TRUE),
                          q1 = quantile(x, 0.25, na.rm = TRUE),
                          q3 = quantile(x, 0.75, na.rm = TRUE))
    )
  }

  trig_summary <- NULL
  if (nrow(results$trigger_bias)) {
    trig_summary <- aggregate(
      rel_error ~ bg_case + model + parameter, data = results$trigger_bias,
      FUN = function(x) c(median = median(x, na.rm = TRUE),
                          q1 = quantile(x, 0.25, na.rm = TRUE),
                          q3 = quantile(x, 0.75, na.rm = TRUE))
    )
  }

  coef_summary <- NULL
  if (nrow(results$bg_coef_bias)) {
    coef_summary <- aggregate(
      error ~ bg_case + coefficient, data = results$bg_coef_bias,
      FUN = function(x) c(median = median(x, na.rm = TRUE),
                          q1 = quantile(x, 0.25, na.rm = TRUE),
                          q3 = quantile(x, 0.75, na.rm = TRUE))
    )
  }

  list(bg_summary = bg_summary, trig_summary = trig_summary, coef_summary = coef_summary)
}

write_nonstat_summary_csv <- function(results = collect_nonstat_results(), out_root = NONSTAT_OUT_ROOT) {
  if (nrow(results$bg_shape)) utils::write.csv(results$bg_shape, file.path(out_root, "table_bg_shape_rmise.csv"), row.names = FALSE)
  if (nrow(results$trigger_bias)) utils::write.csv(results$trigger_bias, file.path(out_root, "table_trigger_bias.csv"), row.names = FALSE)
  if (nrow(results$bg_coef_bias)) utils::write.csv(results$bg_coef_bias, file.path(out_root, "table_bg_coef_bias.csv"), row.names = FALSE)
  if (nrow(results$classification)) utils::write.csv(results$classification, file.path(out_root, "table_classification.csv"), row.names = FALSE)
  if (nrow(results$time_binned)) utils::write.csv(results$time_binned, file.path(out_root, "table_time_binned.csv"), row.names = FALSE)
  if (nrow(results$fit_status)) utils::write.csv(results$fit_status, file.path(out_root, "table_fit_status.csv"), row.names = FALSE)
  ns_log("Tabelle salvate in %s", out_root)
  invisible(TRUE)
}


## ================================================================
## 11. GRAFICI
## ================================================================

MODEL_LABELS <- c(classic = "ETAS", flp = "ETAS-FLP", param = "ETAS-P")
MODEL_COLORS <- c(classic = "#d7301f", flp = "#238b45", param = "#2c7bb6")
BGCASE_LABELS <- c(constant_trend = "Constant + trend lineare", cov1_trend = "Cov. Z1 + trend lineare")

## Fig. A: boxplot rMISEbg per i tre modelli, faceted per scenario spaziale.
## Analogo alla Fig. 2 del paper: risponde a "il drift finisce nel background spaziale?"
plot_bg_shape_recovery <- function(results = collect_nonstat_results()) {
  d <- results$bg_shape
  if (!nrow(d)) { ns_log("Nessun dato per plot_bg_shape_recovery", .level = "WARN"); return(invisible(NULL)) }
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))
  d$bg_case <- factor(d$bg_case, levels = names(BGCASE_LABELS), labels = BGCASE_LABELS)

  p <- ggplot(d, aes(x = model, y = rmise_bg, fill = model)) +
    geom_boxplot(outlier.alpha = 0.4) +
    facet_wrap(~ bg_case, nrow = 1) +
    scale_x_discrete(labels = MODEL_LABELS) +
    scale_fill_manual(values = unname(MODEL_COLORS), labels = MODEL_LABELS, guide = "none") +
    labs(title = "Recupero della forma spaziale del background sotto trend temporale non modellato",
        subtitle = "rMISE della densita' spaziale stimata vs vera f0(s); ci si aspetta NESSUNA differenza sistematica rispetto al caso stazionario",
        x = NULL, y = expression(rMISE[bg])) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
  p
}

## Fig. B: bias relativo dei parametri di triggering (mu, k0, c, p, d, q, magnitude),
## per i tre modelli, faceted per scenario spaziale. Analogo alla Fig. 4 del paper:
## risponde a "il drift finisce nel triggering?"
plot_trigger_bias <- function(results = collect_nonstat_results()) {
  d <- results$trigger_bias
  if (!nrow(d)) { ns_log("Nessun dato per plot_trigger_bias", .level = "WARN"); return(invisible(NULL)) }
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))
  d$bg_case <- factor(d$bg_case, levels = names(BGCASE_LABELS), labels = BGCASE_LABELS)
  d <- d[is.finite(d$rel_error), , drop = FALSE]

  p <- ggplot(d, aes(x = parameter, y = rel_error, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
    facet_wrap(~ bg_case, nrow = 2) +
    scale_fill_manual(values = unname(MODEL_COLORS), labels = MODEL_LABELS, name = NULL) +
    labs(title = "Bias relativo dei parametri ETAS sotto trend temporale non modellato",
        subtitle = "Se il drift viene erroneamente assorbito nel triggering, ci aspettiamo bias sistematico soprattutto su k0, c, p",
        x = NULL, y = "Errore relativo (stima - vero) / vero") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))
  p
}

## Fig. C: bias del coefficiente spaziale Z1 (solo ETAS-P, solo scenario cov1).
## Analogo diretto alla preoccupazione del revisore sul coefficiente di distanza
## dalla faglia nell'applicazione reale.
plot_bg_coef_bias <- function(results = collect_nonstat_results()) {
  d <- results$bg_coef_bias
  if (!nrow(d)) { ns_log("Nessun dato per plot_bg_coef_bias", .level = "WARN"); return(invisible(NULL)) }

  p <- ggplot(d, aes(x = coefficient, y = error)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(fill = MODEL_COLORS[["param"]], alpha = 0.6, outlier.alpha = 0.3) +
    labs(title = "Bias del coefficiente spaziale Z1 in ETAS-P sotto trend temporale non modellato",
        subtitle = "Analogo diretto del coefficiente di distanza dalla faglia nell'applicazione al catalogo italiano",
        x = NULL, y = "Errore (stima - vero)") +
    theme_minimal(base_size = 12)
  p
}

## Fig. D: la diagnostica CHIAVE. Confronta, per bin temporali, la massa
## vera di mu(t) con la massa di probabilita' posteriore di background
## stimata dai tre modelli, mediata sulle repliche. Se le curve stimate
## sono piatte mentre quella vera ha una pendenza marcata, e' la prova
## diretta che il modello non "vede" la non-stazionarieta' temporale.
plot_time_binned_recovery <- function(results = collect_nonstat_results()) {
  d <- results$time_binned
  if (!nrow(d)) { ns_log("Nessun dato per plot_time_binned_recovery", .level = "WARN"); return(invisible(NULL)) }
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))
  d$bg_case <- factor(d$bg_case, levels = names(BGCASE_LABELS), labels = BGCASE_LABELS)

  ## media sulle repliche, per bin, modello e scenario
  agg <- aggregate(
    cbind(fitted_mass_frac, true_mass_frac) ~ bin + t_mid + model + bg_case,
    data = d, FUN = mean, na.rm = TRUE
  )

  ## la massa vera e' identica per costruzione a parita' di t_mid/bg_case;
  ## per il grafico prendiamo una sola curva "true" per pannello
  true_curve <- unique(agg[agg$model == "classic", c("bin", "t_mid", "bg_case", "true_mass_frac")])

  p <- ggplot() +
    geom_line(data = agg, aes(x = t_mid, y = fitted_mass_frac, color = model), linewidth = 1) +
    geom_point(data = agg, aes(x = t_mid, y = fitted_mass_frac, color = model), size = 2) +
    geom_line(data = true_curve, aes(x = t_mid, y = true_mass_frac),
             color = "black", linewidth = 1, linetype = "dashed") +
    facet_wrap(~ bg_case, nrow = 1) +
    scale_color_manual(values = unname(MODEL_COLORS), labels = MODEL_LABELS, name = "Modello stimato") +
    labs(title = "Recupero temporale del background: massa vera vs massa stimata per bin",
        subtitle = "Linea tratteggiata nera = vera mu(t) (trend lineare); linee colorate = massa di probabilita' posteriore di background stimata (media sulle repliche)",
        x = "Tempo (bin, punto medio)", y = "Frazione di massa di background nel bin") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
  p
}

## Salva tutti i grafici come file, e li ritorna in una lista per ispezione interattiva.
build_all_nonstat_plots <- function(results = collect_nonstat_results(), out_root = NONSTAT_OUT_ROOT) {
  fig_dir <- file.path(out_root, "figures")
  ns_dir_create(fig_dir)

  p1 <- plot_bg_shape_recovery(results)
  p2 <- plot_trigger_bias(results)
  p3 <- plot_bg_coef_bias(results)
  p4 <- plot_time_binned_recovery(results)

  if (!is.null(p1)) ggsave(file.path(fig_dir, "fig_A_bg_shape_recovery.png"), p1, width = 9, height = 5, dpi = 300)
  if (!is.null(p2)) ggsave(file.path(fig_dir, "fig_B_trigger_bias.png"), p2, width = 9, height = 7, dpi = 300)
  if (!is.null(p3)) ggsave(file.path(fig_dir, "fig_C_bg_coef_bias.png"), p3, width = 6, height = 5, dpi = 300)
  if (!is.null(p4)) ggsave(file.path(fig_dir, "fig_D_time_binned_recovery.png"), p4, width = 9, height = 5, dpi = 300)

  ns_log("Grafici salvati in %s", fig_dir)
  list(bg_shape_recovery = p1, trigger_bias = p2, bg_coef_bias = p3, time_binned_recovery = p4)
}


## ================================================================
## 12. AUTORUN SWITCH
## ================================================================

if (identical(NONSTAT_MODE, "test")) {
  ns_log("NONSTAT_MODE='test': eseguo una sola replica per bg_case, per verificare che tutto funzioni.")
  run_nonstat_plan(nrep = 1L, force = TRUE)
  res <- collect_nonstat_results()
  print(str(res, max.level = 1))

} else if (identical(NONSTAT_MODE, "run")) {
  ns_log("NONSTAT_MODE='run': eseguo il piano completo (%d repliche per scenario).", NONSTAT_NREP)
  run_nonstat_plan()
  res <- collect_nonstat_results()
  write_nonstat_summary_csv(res)
  cat("\n--- Tasso di successo dei fit (controllare PRIMA di interpretare il bias) ---\n")
  print(summarise_fit_success_rate(res))
  cat("\n--- Riassunto metriche ---\n")
  print(summarise_nonstat_results(res))
  build_all_nonstat_plots(res)

} else if (identical(NONSTAT_MODE, "collect_and_plot")) {
  ns_log("NONSTAT_MODE='collect_and_plot': raccolgo i risultati esistenti e produco tabelle/grafici.")
  res <- collect_nonstat_results()
  write_nonstat_summary_csv(res)
  cat("\n--- Tasso di successo dei fit (controllare PRIMA di interpretare il bias) ---\n")
  print(summarise_fit_success_rate(res))
  cat("\n--- Riassunto metriche ---\n")
  print(summarise_nonstat_results(res))
  build_all_nonstat_plots(res)

} else if (identical(NONSTAT_MODE, "none")) {
  ns_log("NONSTAT_MODE='none': funzioni caricate, nessuna simulazione lanciata.")

} else {
  stop("NONSTAT_MODE sconosciuto: ", NONSTAT_MODE, ". Valori ammessi: test, run, collect_and_plot, none.")
}
