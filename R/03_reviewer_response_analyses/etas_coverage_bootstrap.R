#######################################################################
# COPERTURA DEGLI INTERVALLI DI CONFIDENZA PER I COEFFICIENTI SPAZIALI
# DI ETAS-P: CONDIZIONALE (NAIVE) vs BOOTSTRAP PARAMETRICO COMPLETO
#
# Risposta al punto 3 di Reviewer 1.
#
# LIVELLO A -- copertura naive (condizionale, costo zero)
#   Riusa le repliche gia' fittate nel piano principale (01_main_plan,
#   scenari cov1 e linear_xy, regime Balanced N=1000). Per ciascuna,
#   legge la stima e il SE del GAM (model.bg$mod_global) e verifica se
#   l'intervallo Wald al 95% copre il vero coefficiente. NON lancia
#   nessuna simulazione: legge solo file .rds gia' esistenti.
#
# LIVELLO B -- copertura bootstrap (nuovo calcolo, contenuto)
#   Per un sottoinsieme di M repliche esterne (tra quelle gia' fittate
#   con successo), tratta l'intero vettore stimato (mu, k0, c, p, d, q,
#   magnitude, coefficiente/i di background) come "verita'", genera B
#   cataloghi sintetici da questa verita', rifitta l'INTERA procedura
#   ETAS-P (non solo il passo GAM finale) su ciascuno, e costruisce un
#   intervallo percentile al 95% dalle B stime bootstrap. La copertura
#   e' la frazione delle M repliche esterne per cui questo intervallo
#   contiene il vero coefficiente FISSO dello scenario.
#
# Il coefficiente vero entra nel simulatore come argomento esplicito
# (bg_scenario$sim_spec$bg_lp_coefs), non come costante globale interna:
# per il bootstrap basta quindi costruire uno scenario con questo campo
# sostituito dalla stima, e riusare simulate_one_full()/etas.par.sim_v4
# esattamente come nel piano principale -- nessun simulatore nuovo.
#
# DIPENDENZE: carica (senza autorun) le sole DEFINIZIONI da
#   etas_parametric_final_plan_no_smooth.R
# Non serve lo script FLP: la copertura riguarda solo ETAS-P.
#
# USO:
#   Sys.setenv(COV_MODE = "test")              # 1 esterna, B piccolo, verifica
#   source("etas_coverage_bootstrap.R")
#
#   Sys.setenv(COV_MODE = "run")               # piano completo
#   source("etas_coverage_bootstrap.R")
#
#   Sys.setenv(COV_MODE = "collect_and_plot")  # solo tabelle sui risultati esistenti
#   source("etas_coverage_bootstrap.R")
#
# Variabili opzionali:
#   COV_PLAN_SCRIPT   percorso di etas_parametric_final_plan_no_smooth.R
#   COV_MAIN_OUT_ROOT default: ~/sim_paper1_marco/etas parametrico/etas_parametric_final_no_smooth_outputs_v1
#   COV_OUT_ROOT      dove salvare i risultati di questo script (default: sotto-cartella "05_coverage_bootstrap" di COV_MAIN_OUT_ROOT)
#   COV_M_OUTER       default 30
#   COV_B_INNER       default 100
#######################################################################


## ================================================================
## 0. UTILITY DI BASE
## ================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

cov_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cov_log <- function(..., .level = "INFO") {
  cat(sprintf("[%s] [%s] ", cov_now(), .level), sprintf(...), "\n", sep = "")
  flush.console()
}
cov_dir_create <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
cov_env_chr <- function(name, default = NULL) { x <- Sys.getenv(name, unset = NA_character_); if (is.na(x) || !nzchar(x)) default else x }
cov_env_int <- function(name, default) { x <- Sys.getenv(name, unset = NA_character_); if (is.na(x) || !nzchar(x)) return(default); out <- suppressWarnings(as.integer(x)); if (is.na(out)) default else out }

cov_saveRDS_atomic <- function(object, file, compress = "gzip") {
  cov_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_cov_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) { file.copy(tmp, file, overwrite = TRUE); unlink(tmp) }
  invisible(file)
}

cov_safe_run <- function(expr_fun) {
  t0 <- Sys.time()
  out <- tryCatch(
    list(ok = TRUE, value = expr_fun(), error_message = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error_message = conditionMessage(e))
  )
  out$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}


## ================================================================
## 1. CARICAMENTO SICURO DELLE DEFINIZIONI DEL PIANO PRINCIPALE
## ================================================================

COV_PLAN_SCRIPT <- cov_env_chr("COV_PLAN_SCRIPT", "etas_parametric_final_plan_no_smooth.R")
if (!file.exists(COV_PLAN_SCRIPT)) {
  stop("Non trovo '", COV_PLAN_SCRIPT, "'. Imposta Sys.setenv(COV_PLAN_SCRIPT = '/percorso/completo/file.R').")
}

cov_load_plan_definitions <- function(plan_script = COV_PLAN_SCRIPT, envir = .GlobalEnv) {
  if (exists("etasclass.par", envir = envir, inherits = TRUE) &&
      exists("fit_etas_parametric_full", envir = envir, inherits = TRUE) &&
      exists("BACKGROUND_SCENARIOS_FULL", envir = envir, inherits = TRUE)) {
    cov_log("Definizioni gia' presenti in ambiente: non ricarico.")
    return(invisible(TRUE))
  }
  cov_log("Carico le definizioni da: %s", plan_script)
  txt <- readLines(plan_script, warn = FALSE)
  idx <- grep("^##[[:space:]]*11\\.[[:space:]]*AUTORUN SWITCH", txt)
  if (!length(idx)) idx <- grep("^print_plan_sizes\\(\\)", txt)
  if (!length(idx)) stop("Non trovo il marcatore di AUTORUN SWITCH: mi rifiuto di sourceare l'intero file.")
  code <- paste(txt[seq_len(idx[1] - 1L)], collapse = "\n")
  eval(parse(text = code), envir = envir)
  invisible(TRUE)
}
cov_load_plan_definitions()

required_objs <- c(
  "etas.par.sim_v4", "make_eqcat_from_sim", "fit_etas_parametric_full", "background_fit_spec",
  "BACKGROUND_SCENARIOS_FULL", "COUNT_MAIN_1000", "count_lookup",
  "TRUE_BASE_FULL", "TRUE_BETACOV_FULL", "TRUE_BETA_Z1_FULL", "TRUE_BETA_XSTD_FULL", "TRUE_BETA_YSTD_FULL",
  "M0_FULL", "B_FULL", "LONG_RANGE", "LAT_RANGE", "TMIN_USE", "MAX_EVENTS"
)
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
if (length(missing_objs)) stop("Mancano questi oggetti dopo il caricamento: ", paste(missing_objs, collapse = ", "))
cov_log("Definizioni caricate correttamente.")


## ================================================================
## 2. CONFIGURAZIONE
## ================================================================

COV_MAIN_OUT_ROOT <- cov_env_chr(
  "COV_MAIN_OUT_ROOT",
  path.expand("~/sim_paper1_marco/etas parametrico/etas_parametric_final_no_smooth_outputs_v1")
)
if (!dir.exists(file.path(COV_MAIN_OUT_ROOT, "01_main_plan", "fits"))) {
  stop("Non trovo '", file.path(COV_MAIN_OUT_ROOT, "01_main_plan", "fits"),
       "'. Verifica COV_MAIN_OUT_ROOT (percorso del main plan gia' eseguito).")
}

COV_OUT_ROOT <- cov_env_chr("COV_OUT_ROOT", file.path(COV_MAIN_OUT_ROOT, "05_coverage_bootstrap"))
COV_M_OUTER  <- cov_env_int("COV_M_OUTER", 30L)
COV_B_INNER  <- cov_env_int("COV_B_INNER", 100L)
COV_N_CORES  <- cov_env_int("COV_N_CORES", max(1L, parallel::detectCores() - 1L))
COV_MODE     <- tolower(cov_env_chr("COV_MODE", "none"))
COV_SEED_BASE <- cov_env_int("COV_SEED_BASE", 91000000L)

COV_SCENARIOS <- list(
  cov1      = list(bg_case = "cov1",      coef_names = "Z1",              true_coefs = c(Z1 = TRUE_BETA_Z1_FULL)),
  linear_xy = list(bg_case = "linear_xy", coef_names = c("x_std","y_std"), true_coefs = c(x_std = TRUE_BETA_XSTD_FULL, y_std = TRUE_BETA_YSTD_FULL))
)
COV_COUNT_CASE <- "N1000_balanced_50_50"

cov_dir_create(COV_OUT_ROOT)
cov_dir_create(file.path(COV_OUT_ROOT, "bootstrap_fits"))
cov_dir_create(file.path(COV_OUT_ROOT, "logs"))

cov_log("Config: main_plan=%s | out=%s | M_outer=%d | B_inner=%d | n_cores=%d | scenarios=%s",
        COV_MAIN_OUT_ROOT, COV_OUT_ROOT, COV_M_OUTER, COV_B_INNER, COV_N_CORES,
        paste(names(COV_SCENARIOS), collapse = ","))

## Log di avanzamento su file condiviso: con l'esecuzione parallela,
## l'output a console dei worker puo' non arrivare (o arrivare
## disordinato) al processo principale, specialmente con backend PSOCK
## su Windows. Ogni job scrive quindi anche una riga qui, indipendente
## dagli altri (append di una singola riga: rischio di interleaving
## trascurabile per un semplice log di avanzamento, non per i dati).
cov_append_progress <- function(..., out_root = COV_OUT_ROOT) {
  line <- sprintf("[%s] %s", cov_now(), sprintf(...))
  cat(line, "\n", file = file.path(out_root, "logs", "progress_bootstrap.log"), append = TRUE)
}


## ================================================================
## 3. LIVELLO A -- COPERTURA NAIVE (CONDIZIONALE, SOLA LETTURA)
## ================================================================

## Elenca i file del main plan per un dato bg_case/count_case, filtrando
## per metadata (robusto rispetto al formato esatto del nome file).
list_main_plan_files <- function(bg_case, count_case, out_root = COV_MAIN_OUT_ROOT) {
  fits_dir <- file.path(out_root, "01_main_plan", "fits")
  all_files <- list.files(fits_dir, pattern = "\\.rds$", full.names = TRUE)
  keep <- logical(length(all_files))
  for (i in seq_along(all_files)) {
    meta <- tryCatch(readRDS(all_files[i])$metadata, error = function(e) NULL)
    keep[i] <- !is.null(meta) && identical(meta$bg_case, bg_case) && identical(meta$count_case, count_case)
  }
  all_files[keep]
}

## Estrae stima, SE (dalla tabella parametrica del GAM) e costruisce la
## verifica di copertura Wald per un singolo file gia' fittato.
naive_ci_from_file <- function(f, coef_names, true_coefs) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj) || !isTRUE(obj$fits$param_correct$ok)) return(NULL)
  fit <- obj$fits$param_correct$value
  mod <- fit$model.bg$mod_global
  if (is.null(mod)) return(NULL)

  ptab <- tryCatch(summary(mod)$p.table, error = function(e) NULL)
  if (is.null(ptab) || !all(coef_names %in% rownames(ptab))) return(NULL)

  do.call(rbind, lapply(coef_names, function(cn) {
    est <- ptab[cn, "Estimate"]
    se  <- ptab[cn, "Std. Error"]
    lo <- est - 1.96 * se; hi <- est + 1.96 * se
    data.frame(
      rep = obj$metadata$rep, coefficient = cn,
      true_value = as.numeric(true_coefs[cn]),
      estimate = est, se = se, ci_lo = lo, ci_hi = hi,
      covered = (as.numeric(true_coefs[cn]) >= lo) & (as.numeric(true_coefs[cn]) <= hi),
      width = hi - lo,
      stringsAsFactors = FALSE
    )
  }))
}

## Copertura naive per uno scenario: legge TUTTI i file gia' fittati con
## successo, nessuna nuova simulazione.
naive_coverage_one_scenario <- function(scenario_key, count_case = COV_COUNT_CASE) {
  spec <- COV_SCENARIOS[[scenario_key]]
  files <- list_main_plan_files(spec$bg_case, count_case)
  cov_log("Livello A [%s]: %d file trovati per bg_case=%s, count_case=%s",
          scenario_key, length(files), spec$bg_case, count_case)
  rows <- lapply(files, naive_ci_from_file, coef_names = spec$coef_names, true_coefs = spec$true_coefs)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$scenario <- scenario_key
  out
}


## ================================================================
## 4. LIVELLO B -- BOOTSTRAP PARAMETRICO DELL'INTERA PROCEDURA
## ================================================================

## Simulatore per il bootstrap: identico a simulate_one_full() del piano
## principale, ma con l'intero vettore di parametri (incluso il/i
## coefficiente/i di background) passato esplicitamente come argomento
## anziche' letto da costanti globali fisse (TRUE_BASE_FULL, ecc.).
## Riusa etas.par.sim_v4 senza modifiche.
simulate_one_bootstrap <- function(bg_case, bg_lp_coefs_hat, count_regime, seed,
                                   mu_hat, k0_hat, c_hat, p_hat, d_hat, q_hat, betacov_hat,
                                   return_bg_info = FALSE, return_support = FALSE) {
  params_now <- c(mu = mu_hat, k0 = k0_hat, c = c_hat, p = p_hat, gamma = 0, d = d_hat, q = q_hat)

  bg_scenario <- BACKGROUND_SCENARIOS_FULL[[bg_case]]
  bg_scenario$sim_spec$bg_lp_coefs <- bg_lp_coefs_hat

  sim_args <- list(
    params = params_now, m0 = M0_FULL, b = B_FULL,
    tmin = TMIN_USE, t.lag = count_regime$T_lag,
    long.range = LONG_RANGE, lat.range = LAT_RANGE, longlat.to.km = TRUE,
    sectoday = FALSE,
    trig_lp_type = "linear", trig_lp_coefs = c(m_rel = betacov_hat),
    n_support_obs = 1500, mult_support = 4, max_events = MAX_EVENTS,
    return_bg_info = return_bg_info, return_support = return_support,
    plot_bg = FALSE, plot_catalog = FALSE, seed = seed
  )
  sim_args <- c(sim_args, bg_scenario$sim_spec)
  do.call(etas.par.sim_v4, sim_args)
}

## Estrae il vettore di parametri stimati "come se fossero la verita'"
## da una replica gia' fittata del main plan. Include un controllo di
## plausibilita': repliche con parametri di triggering troppo lontani
## dai veri valori della simulazione (fit numericamente degeneri, gia'
## osservati occasionalmente nello studio principale) vengono scartate,
## perche' bootstrappare da una "verita'" instabile non e' significativo
## e puo' produrre cataloghi esplosivi enormi da fittare ripetutamente.
extract_theta_hat <- function(fit_obj, coef_names, max_ratio = 5) {
  p <- fit_obj$params.MLtot
  mod <- fit_obj$model.bg$mod_global
  if (is.null(p) || is.null(mod)) return(NULL)
  cf <- stats::coef(mod)
  if (!all(coef_names %in% names(cf))) return(NULL)

  ## Nota: mu e k0 NON hanno un valore "vero" globale con cui confrontarli
  ## -- sono calibrati scenario per scenario (TRUE_BASE_FULL li definisce
  ## come NA per costruzione). Il controllo di plausibilita' si applica
  ## solo ai parametri di triggering fissi in tutti gli scenari.
  check_params <- c("c", "d", "q")
  for (nm in check_params) {
    true_val <- as.numeric(TRUE_BASE_FULL[nm])
    est_val <- as.numeric(p[nm])
    if (!is.finite(true_val) || !is.finite(est_val) || true_val == 0) next
    ratio <- est_val / true_val
    if (!is.finite(ratio) || ratio < 1 / max_ratio || ratio > max_ratio) {
      return(NULL)  # scartata: stima implausibile, probabile fit degenere
    }
  }

  list(
    mu_hat = as.numeric(p["mu"]), k0_hat = as.numeric(p["k0"]),
    c_hat = as.numeric(p["c"]), p_hat = as.numeric(p["p"]),
    d_hat = as.numeric(p["d"]), q_hat = as.numeric(p["q"]),
    betacov_hat = as.numeric(p["magnitude"]),
    bg_lp_coefs_hat = cf[coef_names]
  )
}

## Esegue il bootstrap interno (B repliche) per UNA replica esterna gia'
## fittata, e restituisce le B stime bootstrap dei coefficienti di
## background, piu' l'intervallo percentile al 95% e l'indicatore di
## copertura rispetto al vero coefficiente FISSO dello scenario.
run_bootstrap_for_replicate <- function(scenario_key, theta_hat, count_regime,
                                        B = COV_B_INNER, seed_base, verbose = TRUE,
                                        progress_label = NULL, log_every = 10L) {
  spec <- COV_SCENARIOS[[scenario_key]]
  coef_names <- spec$coef_names

  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(coef_names),
                     dimnames = list(NULL, coef_names))
  n_ok <- 0L
  t_start <- Sys.time()

  starts <- list(
    mu = theta_hat$mu_hat, k0 = theta_hat$k0_hat, c = theta_hat$c_hat, p = theta_hat$p_hat,
    gamma = 0, d = theta_hat$d_hat, q = theta_hat$q_hat, betacov = theta_hat$betacov_hat
  )

  for (b in seq_len(B)) {
    seed_b <- seed_base + b
    sim_res <- cov_safe_run(function() {
      simulate_one_bootstrap(
        bg_case = spec$bg_case, bg_lp_coefs_hat = theta_hat$bg_lp_coefs_hat,
        count_regime = count_regime, seed = seed_b,
        mu_hat = theta_hat$mu_hat, k0_hat = theta_hat$k0_hat,
        c_hat = theta_hat$c_hat, p_hat = theta_hat$p_hat,
        d_hat = theta_hat$d_hat, q_hat = theta_hat$q_hat, betacov_hat = theta_hat$betacov_hat
      )
    })
    if (!isTRUE(sim_res$ok)) next

    cat_boot_df <- tryCatch(make_eqcat_from_sim(sim_res$value), error = function(e) NULL)
    if (is.null(cat_boot_df)) next

    ## Guardia: cataloghi esplosi o anomalmente grandi vengono scartati
    ## subito invece di tentare un fit costosissimo che potrebbe
    ## richiedere molto piu' tempo del previsto. La soglia (5x la
    ## dimensione tipica N=1000 del disegno principale) e' volutamente
    ## larga per non scartare variabilita' campionaria normale.
    if (isTRUE(sim_res$value$exploded) || nrow(cat_boot_df) > 5000L) {
      if (!is.null(progress_label)) {
        cov_append_progress("[%s] iterazione %d/%d SCARTATA: catalogo esploso o troppo grande (n=%d)",
                            progress_label, b, B, nrow(cat_boot_df), out_root = COV_OUT_ROOT)
      }
      next
    }

    fit_res <- cov_safe_run(function() {
      fit_etas_parametric_full(cat_boot_df, bg_case = spec$bg_case, formula_key = "correct", starts = starts)
    })
    if (!isTRUE(fit_res$ok)) next

    cf_boot <- tryCatch(stats::coef(fit_res$value$model.bg$mod_global), error = function(e) NULL)
    if (is.null(cf_boot) || !all(coef_names %in% names(cf_boot))) next

    boot_mat[b, ] <- cf_boot[coef_names]
    n_ok <- n_ok + 1L

    if (!is.null(progress_label) && (b %% log_every == 0L || b == B)) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      rate <- elapsed / b
      eta_min <- rate * (B - b)
      cov_append_progress("[%s] iterazione %d/%d (successi: %d) | %.1f min trascorsi | ~%.1f min al termine di questo job",
                          progress_label, b, B, n_ok, elapsed, eta_min, out_root = COV_OUT_ROOT)
    }
  }

  if (isTRUE(verbose)) cov_log("  bootstrap interno: %d/%d fit riusciti", n_ok, B)

  out <- lapply(coef_names, function(cn) {
    vals <- boot_mat[, cn]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 10) {
      return(data.frame(coefficient = cn, n_boot_ok = length(vals),
                        ci_lo = NA_real_, ci_hi = NA_real_, covered = NA))
    }
    qs <- quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
    true_val <- as.numeric(spec$true_coefs[cn])
    data.frame(
      coefficient = cn, n_boot_ok = length(vals),
      ci_lo = unname(qs[1]), ci_hi = unname(qs[2]),
      covered = (true_val >= qs[1]) && (true_val <= qs[2])
    )
  })
  out <- do.call(rbind, out)
  attr(out, "boot_values") <- boot_mat
  out
}

## Orchestrazione completa del Livello B, con checkpoint per replica
## esterna (resumable: se il file esiste gia', lo salta), eseguita in
## parallelo sui ~60 job (scenario x replica esterna). Il ciclo bootstrap
## interno (B) resta sequenziale DENTRO ciascun job, per evitare
## parallelismo annidato ed eventuale sovra-sottoscrizione dei core.
##
## Backend: parallel::mclapply (fork) su Unix/macOS -- niente da
## esportare, i worker ereditano l'intero ambiente per copy-on-write.
## Su Windows (dove fork non e' disponibile), usa un cluster PSOCK con
## clusterExport dell'intero ambiente globale. Con COV_N_CORES = 1
## resta sequenziale (comportamento originale, utile per debug).
cov_run_parallel <- function(job_list, worker_fun, n_cores = COV_N_CORES) {
  if (n_cores <= 1L || length(job_list) <= 1L) {
    return(lapply(job_list, worker_fun))
  }
  is_windows <- .Platform$OS.type == "windows"
  if (!is_windows) {
    parallel::mclapply(job_list, worker_fun, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, { library(mgcv) })
    parallel::parLapply(cl, job_list, worker_fun)
  }
}

run_coverage_bootstrap_plan <- function(scenario_keys = names(COV_SCENARIOS),
                                        M_outer = COV_M_OUTER, B_inner = COV_B_INNER,
                                        count_case = COV_COUNT_CASE,
                                        out_root = COV_OUT_ROOT, seed_base = COV_SEED_BASE,
                                        n_cores = COV_N_CORES, force = FALSE) {
  count_regime <- count_lookup(count_case)

  ## Costruisce la lista piatta dei job (scenario x replica esterna
  ## scelta), filtrando quelli gia' completati (a meno di force = TRUE).
  job_list <- list()
  for (scenario_key in scenario_keys) {
    spec <- COV_SCENARIOS[[scenario_key]]
    files <- list_main_plan_files(spec$bg_case, count_case)

    ok_files <- Filter(function(f) {
      obj <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(obj) || !isTRUE(obj$fits$param_correct$ok)) return(FALSE)
      fit_val <- obj$fits$param_correct$value
      if (is.null(fit_val$model.bg$mod_global)) return(FALSE)
      !is.null(extract_theta_hat(fit_val, spec$coef_names))  # scarta anche i fit numericamente implausibili
    }, files)

    cov_log("Livello B [%s]: %d/%d repliche disponibili con fit riuscito", scenario_key, length(ok_files), length(files))
    if (length(ok_files) < M_outer) {
      cov_log("  Attenzione: disponibili meno repliche (%d) di M_outer richiesto (%d); uso tutte quelle disponibili.",
              length(ok_files), M_outer, .level = "WARN")
    }

    set.seed(seed_base + match(scenario_key, names(COV_SCENARIOS)))
    chosen <- sample(ok_files, size = min(M_outer, length(ok_files)))

    for (f_outer in chosen) {
      outer_rep_id <- readRDS(f_outer)$metadata$rep
      out_file <- file.path(out_root, "bootstrap_fits",
                            sprintf("%s__outerrep%03d.rds", scenario_key, outer_rep_id))
      if (file.exists(out_file) && !isTRUE(force)) next

      job_list[[length(job_list) + 1L]] <- list(
        scenario_key = scenario_key, f_outer = f_outer, outer_rep_id = outer_rep_id,
        out_file = out_file, count_regime = count_regime, B_inner = B_inner,
        seed_base_job = seed_base + 100000L * match(scenario_key, names(COV_SCENARIOS)) + 1000L * outer_rep_id
      )
    }
  }

  cov_log("Totale job da eseguire: %d (gia' completati: saltati automaticamente)", length(job_list))
  if (!length(job_list)) return(invisible(TRUE))

  worker_fun <- function(job) {
    spec <- COV_SCENARIOS[[job$scenario_key]]
    obj_outer <- readRDS(job$f_outer)

    theta_hat <- extract_theta_hat(obj_outer$fits$param_correct$value, spec$coef_names)
    if (is.null(theta_hat)) {
      cov_append_progress("[%s] outerrep %03d: impossibile estrarre theta_hat, salto.",
                          job$scenario_key, job$outer_rep_id, out_root = COV_OUT_ROOT)
      return(invisible(NULL))
    }

    result <- run_bootstrap_for_replicate(job$scenario_key, theta_hat, job$count_regime,
                                          B = job$B_inner, seed_base = job$seed_base_job, verbose = FALSE,
                                          progress_label = sprintf("%s/outerrep%03d", job$scenario_key, job$outer_rep_id))

    out_obj <- list(
      metadata = list(scenario = job$scenario_key, bg_case = spec$bg_case, count_case = COV_COUNT_CASE,
                      outer_rep = job$outer_rep_id, M_outer = COV_M_OUTER, B_inner = job$B_inner,
                      source_file = job$f_outer),
      theta_hat = theta_hat,
      naive = naive_ci_from_file(job$f_outer, spec$coef_names, spec$true_coefs),
      bootstrap = result,
      boot_values = attr(result, "boot_values")
    )
    cov_saveRDS_atomic(out_obj, job$out_file)
    cov_append_progress("[%s] outerrep %03d: completato (%d/%d bootstrap ok).",
                        job$scenario_key, job$outer_rep_id,
                        sum(!is.na(attr(result, "boot_values")[, 1])), job$B_inner, out_root = COV_OUT_ROOT)
    invisible(NULL)
  }

  cov_run_parallel(job_list, worker_fun, n_cores = n_cores)
  invisible(TRUE)
}


## ================================================================
## 5. RACCOLTA E TABELLA RIASSUNTIVA FINALE (Tabella S12)
## ================================================================

collect_coverage_results <- function(out_root = COV_OUT_ROOT) {
  files <- list.files(file.path(out_root, "bootstrap_fits"), pattern = "\\.rds$", full.names = TRUE)
  naive_rows <- list(); boot_rows <- list()

  for (f in files) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) next
    if (!is.null(obj$naive)) {
      d <- obj$naive; d$scenario <- obj$metadata$scenario
      naive_rows[[length(naive_rows) + 1L]] <- d
    }
    if (!is.null(obj$bootstrap)) {
      d <- obj$bootstrap; d$scenario <- obj$metadata$scenario; d$outer_rep <- obj$metadata$outer_rep
      boot_rows[[length(boot_rows) + 1L]] <- d
    }
  }
  rb <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  list(naive_subset = rb(naive_rows), bootstrap = rb(boot_rows))
}

## Copertura naive calcolata su TUTTE le repliche disponibili nel main
## plan (non solo il sottoinsieme M_outer usato per il bootstrap): e'
## il numero di riferimento piu' solido per il confronto, dato che non
## richiede nuove simulazioni ed e' quindi calcolabile sull'intero
## campione di 100 repliche.
build_table_S12 <- function(out_root_coverage = COV_OUT_ROOT, count_case = COV_COUNT_CASE) {
  boot_res <- collect_coverage_results(out_root_coverage)

  rows <- list()
  for (scenario_key in names(COV_SCENARIOS)) {
    spec <- COV_SCENARIOS[[scenario_key]]
    naive_full <- naive_coverage_one_scenario(scenario_key, count_case)

    boot_sub <- boot_res$bootstrap[boot_res$bootstrap$scenario == scenario_key, ]

    for (cn in spec$coef_names) {
      naive_cn <- naive_full[naive_full$coefficient == cn, ]
      boot_cn  <- boot_sub[boot_sub$coefficient == cn, ]

      rows[[length(rows) + 1L]] <- data.frame(
        scenario = scenario_key, coefficient = cn, true_value = as.numeric(spec$true_coefs[cn]),
        n_naive = nrow(naive_cn),
        naive_coverage = if (nrow(naive_cn)) mean(naive_cn$covered, na.rm = TRUE) else NA_real_,
        naive_mean_width = if (nrow(naive_cn)) mean(naive_cn$width, na.rm = TRUE) else NA_real_,
        n_boot = nrow(boot_cn),
        boot_coverage = if (nrow(boot_cn)) mean(boot_cn$covered, na.rm = TRUE) else NA_real_,
        boot_mean_width = if (nrow(boot_cn)) mean(boot_cn$ci_hi - boot_cn$ci_lo, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


## ================================================================
## 6. AUTORUN SWITCH
## ================================================================

if (identical(COV_MODE, "test")) {
  cov_log("COV_MODE='test': 1 esterna, B=10, seriale (n_cores=1) per verificare che tutto funzioni.")
  run_coverage_bootstrap_plan(M_outer = 1L, B_inner = 10L, n_cores = 1L, force = TRUE)
  print(collect_coverage_results())

} else if (identical(COV_MODE, "run")) {
  cov_log("COV_MODE='run': piano completo (M_outer=%d, B_inner=%d, n_cores=%d).", COV_M_OUTER, COV_B_INNER, COV_N_CORES)
  run_coverage_bootstrap_plan()
  tab <- build_table_S12()
  utils::write.csv(tab, file.path(COV_OUT_ROOT, "table_S12_coverage.csv"), row.names = FALSE)
  print(tab, row.names = FALSE, digits = 3)

} else if (identical(COV_MODE, "collect_and_plot")) {
  cov_log("COV_MODE='collect_and_plot': ricalcolo la tabella sui risultati esistenti.")
  tab <- build_table_S12()
  utils::write.csv(tab, file.path(COV_OUT_ROOT, "table_S12_coverage.csv"), row.names = FALSE)
  print(tab, row.names = FALSE, digits = 3)

} else if (identical(COV_MODE, "none")) {
  cov_log("COV_MODE='none': funzioni caricate, nessun calcolo lanciato.")

} else {
  stop("COV_MODE sconosciuto: ", COV_MODE, ". Valori ammessi: test, run, collect_and_plot, none.")
}
