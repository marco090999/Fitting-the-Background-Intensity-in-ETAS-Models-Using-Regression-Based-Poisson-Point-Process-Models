#######################################################################
# BOOTSTRAP PARAMETRICO DELL'INTERA PROCEDURA ETAS-P SULL'APPLICAZIONE
# (catalogo italiano reale, M0=2.5 -- il fit riportato in Tabella 1/2)
#
# Risposta al punto 3 di Reviewer 1, applicata direttamente al catalogo
# reale (complemento della validazione gia' fatta su dati simulati in
# SI Sec. 2.6).
#
# DISEGNO
# --------
# 1) Riusa il fit ETAS-P gia' ottenuto a M0=2.5 (da
#    m0_sensitivity_raw_fits.rds) come "verita'" per il bootstrap:
#    l'intero vettore stimato (mu, k0, c, p, d, q, alpha, beta_distmin
#    nel triggering, theta0/superficie smooth/theta_distmin nel
#    background) viene trattato come se fosse il vero processo
#    generatore.
# 2) Il background di ETAS-P qui e' una superficie GAM stimata (non una
#    funzione algebrica nota come nel piano di simulazione): le
#    posizioni di background vengono campionate per rifiuto-accettazione
#    valutando la superficie fittata (intercetta + smooth(x,y) +
#    theta_distmin * distmin) su una griglia.
# 3) La covariata distmin (distanza dalla faglia piu' vicina) non ha
#    una forma nota per punti nuovi: il pacchetto etasFLP distribuisce
#    solo la colonna gia' calcolata, non la geometria delle faglie
#    sorgente (verificato dalla documentazione di catalog.withcov, che
#    non cita alcun database di faglie). Costruiamo quindi un
#    interpolatore spaziale liscio (GAM) di distmin sui punti osservati,
#    fittato UNA VOLTA SOLA sui dati reali e riusato per assegnare un
#    valore plausibile di distmin a qualunque posizione simulata. Questo
#    va documentato esplicitamente come approssimazione nel testo SI:
#    la distanza da una faglia e' una funzione spazialmente liscia
#    (tranne esattamente sulle linee di faglia), quindi l'interpolazione
#    dai valori osservati dovrebbe approssimare bene il valore vero per
#    posizioni simulate nella stessa area geografica del catalogo.
# 4) B repliche bootstrap: per ciascuna, si simula un catalogo sintetico
#    dalla "verita'" cosi' definita, e si rifitta l'INTERA procedura
#    ETAS-P (non solo il passo GAM finale) con la stessa identica
#    chiamata usata per il fit originale (stesso M0, stesse impostazioni).
#
# DIPENDENZE: richiede in sessione le definizioni di etasclass.par()
# (dal file sorgente del piano di simulazione) e l'oggetto
# catalog.withcov. Legge m0_sensitivity_raw_fits.rds gia' presente.
#
# USO:
#   Sys.setenv(APPBOOT_MODE = "test")   # 2 repliche, seriale, verifica
#   source("etas_application_bootstrap.R")
#
#   Sys.setenv(APPBOOT_N_CORES = "50")
#   Sys.setenv(APPBOOT_MODE = "run")    # piano completo, parallelo
#   source("etas_application_bootstrap.R")
#######################################################################

suppressPackageStartupMessages({ library(mgcv); library(parallel) })

ab_log <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }
`%||%` <- function(a, b) if (!is.null(a)) a else b
ab_dir_create <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }

ab_saveRDS_atomic <- function(object, file, compress = "gzip") {
  ab_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_ab_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) { file.copy(tmp, file, overwrite = TRUE); unlink(tmp) }
  invisible(file)
}

ab_safe_run <- function(expr_fun) {
  t0 <- Sys.time()
  out <- tryCatch(list(ok = TRUE, value = expr_fun(), error_message = NA_character_),
                  error = function(e) list(ok = FALSE, value = NULL, error_message = conditionMessage(e)))
  out$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}


## ================================================================
## 1. CARICAMENTO SICURO DI etasclass.par() E COSTANTI DI BASE
## ================================================================

AB_PLAN_SCRIPT <- Sys.getenv("APPBOOT_PLAN_SCRIPT", unset = "etas_parametric_final_plan_no_smooth.R")
if (!exists("etasclass.par", mode = "function", inherits = TRUE)) {
  if (!file.exists(AB_PLAN_SCRIPT)) stop("Non trovo '", AB_PLAN_SCRIPT, "'. Imposta Sys.setenv(APPBOOT_PLAN_SCRIPT = '/percorso/completo/file.R').")
  ab_log("Carico le definizioni da: %s", AB_PLAN_SCRIPT)
  txt <- readLines(AB_PLAN_SCRIPT, warn = FALSE)
  idx <- grep("^##[[:space:]]*11\\.[[:space:]]*AUTORUN SWITCH", txt)
  if (!length(idx)) idx <- grep("^print_plan_sizes\\(\\)", txt)
  if (!length(idx)) stop("Non trovo il marcatore di AUTORUN SWITCH.")
  eval(parse(text = paste(txt[seq_len(idx[1] - 1L)], collapse = "\n")), envir = .GlobalEnv)
}
if (!exists("catalog.withcov")) stop("Non trovo 'catalog.withcov' in sessione.")

AB_MAIN_DIR <- Sys.getenv("APPBOOT_MAIN_DIR", unset = "/home/nicolettadangelo/sim_paper1_marco/etas parametrico")
AB_RAW_FITS_FILE <- file.path(AB_MAIN_DIR, "m0_sensitivity_raw_fits.rds")
if (!file.exists(AB_RAW_FITS_FILE)) stop("Non trovo '", AB_RAW_FITS_FILE, "'.")

AB_OUT_DIR <- Sys.getenv("APPBOOT_OUT_DIR", unset = file.path(AB_MAIN_DIR, "application_bootstrap_outputs"))
AB_B <- as.integer(Sys.getenv("APPBOOT_B", unset = "100"))
AB_N_CORES <- as.integer(Sys.getenv("APPBOOT_N_CORES", unset = as.character(max(1L, parallel::detectCores() - 1L))))
AB_MODE <- tolower(Sys.getenv("APPBOOT_MODE", unset = "none"))
AB_SEED_BASE <- as.integer(Sys.getenv("APPBOOT_SEED_BASE", unset = "77300000"))

ab_dir_create(AB_OUT_DIR)
ab_dir_create(file.path(AB_OUT_DIR, "boot_fits"))
ab_dir_create(file.path(AB_OUT_DIR, "logs"))

UNIT_KM <- 6371.3 * pi / 180   # stessa convenzione usata ovunque nel codice per long/lat -> km


## ================================================================
## 2. ESTRAZIONE DI theta_hat DAL FIT GIA' ESISTENTE (M0=2.5)
## ================================================================

ab_load_theta_hat <- function() {
  raw_fits <- readRDS(AB_RAW_FITS_FILE)
  fit <- raw_fits[["2.5"]]$param
  if (is.null(fit)) stop("Non trovo il fit ETAS-P a M0=2.5 dentro '", AB_RAW_FITS_FILE, "'.")

  p <- fit$params.MLtot
  mod <- fit$model.bg$mod_global
  if (is.null(mod)) stop("Il fit non contiene model.bg$mod_global.")
  if (!("distmin" %in% names(p))) stop("params.MLtot non contiene 'distmin': controlla formula1 del fit originale.")

  list(
    fit_param = fit,
    mu_hat = as.numeric(p["mu"]), k0_hat = as.numeric(p["k0"]),
    c_hat = as.numeric(p["c"]), p_hat = as.numeric(p["p"]),
    d_hat = as.numeric(p["d"]), q_hat = as.numeric(p["q"]),
    betacov_hat = as.numeric(p["magnitude"]),
    beta_distmin_hat = as.numeric(p["distmin"]),
    mod_global = mod
  )
}

ab_theta_hat <- ab_load_theta_hat()
ab_log("theta_hat caricato: mu=%.4f k0=%.4f c=%.4f p=%.4f d=%.4f q=%.4f alpha=%.4f beta_distmin=%.4f",
       ab_theta_hat$mu_hat, ab_theta_hat$k0_hat, ab_theta_hat$c_hat, ab_theta_hat$p_hat,
       ab_theta_hat$d_hat, ab_theta_hat$q_hat, ab_theta_hat$betacov_hat, ab_theta_hat$beta_distmin_hat)


## ================================================================
## 3. FINESTRA SPAZIO-TEMPORALE E DISTRIBUZIONE DI MAGNITUDO
## ================================================================
## tmin = tempo del primo evento nel catalogo filtrato per M0 (identico
## alla convenzione interna di etasclass.par, verificata nel sorgente:
## "tmin <- min(times.tot)"). Con M0=2.5 il catalogo filtrato coincide
## con quello intero.

AB_M0 <- 2.5
AB_TMIN <- min(catalog.withcov$time[catalog.withcov$magn1 >= AB_M0])
AB_TMAX <- max(catalog.withcov$time[catalog.withcov$magn1 >= AB_M0])
AB_TLAG <- AB_TMAX - AB_TMIN
AB_LONG_RANGE <- range(catalog.withcov$long)
AB_LAT_RANGE <- range(catalog.withcov$lat)

## b-value generatore delle magnitudo sintetiche: stimatore ML di
## Aki/Utsu alla soglia M0, con la stessa risoluzione dM gia' usata
## nell'analisi di completezza (nuova Sezione 3 del SI).
ab_detect_dm <- function(mags, candidates = c(0.01, 0.02, 0.1, 0.2)) {
  for (dm in sort(candidates, decreasing = TRUE)) {
    if (all(abs(mags / dm - round(mags / dm)) < 1e-6)) return(dm)
  }
  0.1
}
AB_DM <- ab_detect_dm(catalog.withcov$magn1)
ab_m_above <- catalog.withcov$magn1[catalog.withcov$magn1 >= AB_M0]
AB_B_GR <- log10(exp(1)) / (mean(ab_m_above) - (AB_M0 - AB_DM / 2))   # b-value
AB_BETA_GR <- log(10) * AB_B_GR                                      # tasso esponenziale per generare magnitudo
ab_log("Finestra: tmin=%.3f tmax=%.3f (t.lag=%.3f) | dM=%.3f | b=%.4f", AB_TMIN, AB_TMAX, AB_TLAG, AB_DM, AB_B_GR)


## ================================================================
## 4. INTERPOLATORE SPAZIALE DI distmin (fittato una volta sola)
## ================================================================

ab_build_distmin_interpolator <- function(cat_data = catalog.withcov, unit = UNIT_KM) {
  x_km <- cat_data$long * unit
  y_km <- cat_data$lat * unit
  df <- data.frame(x = x_km, y = y_km, distmin = cat_data$distmin)
  mod <- mgcv::gam(distmin ~ s(x, y, k = 30), data = df)
  list(mod = mod, residuals = as.numeric(residuals(mod, type = "response")))
}
AB_DISTMIN_INTERP <- ab_build_distmin_interpolator()
ab_log("Interpolatore di distmin costruito (R^2 adj = %.3f, SD residui = %.3f).",
       summary(AB_DISTMIN_INTERP$mod)$r.sq, sd(AB_DISTMIN_INTERP$residuals))

## IMPORTANTE: la sola previsione smooth spiega solo il 66% circa della
## varianza di distmin nei dati osservati (verificato empiricamente).
## Usare solo la previsione smooth eliminerebbe la componente di
## variazione "indipendente dalla posizione" che nei dati reali e' in
## parte cio' che identifica theta_distmin separatamente dal termine
## smooth s(x,y) del modello di background -- con conseguenza osservata
## empiricamente in un primo test: theta_distmin_hat sistematicamente
## di segno opposto e ampiezza doppia rispetto al vero valore su
## repliche indipendenti (bias strutturale, non rumore campionario).
## Correggiamo reiniettando un residuo ricampionato dai residui
## osservati (bootstrap iid dei residui): approssimazione che restaura
## la stessa VARIANZA "extra-spaziale" presente nei dati veri, pur non
## replicando l'eventuale autocorrelazione spaziale fine dei residui
## (approssimazione ulteriore, dichiarata esplicitamente).
ab_predict_distmin <- function(x_km, y_km, add_noise = TRUE) {
  pred <- as.numeric(predict(AB_DISTMIN_INTERP$mod, newdata = data.frame(x = x_km, y = y_km)))
  if (isTRUE(add_noise)) {
    pred <- pred + sample(AB_DISTMIN_INTERP$residuals, length(pred), replace = TRUE)
  }
  pmax(0, pred)
}


## ================================================================
## 5. SIMULATORE: BACKGROUND DA SUPERFICIE GAM STIMATA + TRIGGERING
## ================================================================

simulate_one_application_bootstrap <- function(theta_hat, tmin, tmax, long_range, lat_range,
                                                m0, beta_gr, unit = UNIT_KM,
                                                n_support_obs = 1500, mult_support = 4,
                                                max_events = 50000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  xmin <- long_range[1] * unit; xmax <- long_range[2] * unit
  ymin <- lat_range[1] * unit; ymax <- lat_range[2] * unit
  area_km2 <- (xmax - xmin) * (ymax - ymin)
  t_lag <- tmax - tmin

  mod <- theta_hat$mod_global

  ## ---- supporto per stimare eta_max (serve al rifiuto-accettazione) ----
  n_sup <- n_support_obs * (1 + mult_support)
  x_sup <- runif(n_sup, xmin, xmax)
  y_sup <- runif(n_sup, ymin, ymax)
  distmin_sup <- ab_predict_distmin(x_sup, y_sup)
  eta_sup <- as.numeric(predict(mod, newdata = data.frame(x = x_sup, y = y_sup, distmin = distmin_sup), type = "link"))
  eta_max <- max(eta_sup, na.rm = TRUE)

  ## ---- numero di eventi di background: Poisson(mu * t_lag) ----
  n0 <- rpois(1, theta_hat$mu_hat * t_lag)
  if (n0 == 0L) return(list(cat.sim = data.frame(), n0 = 0L, nson = 0L, exploded = FALSE))

  t0 <- runif(n0, tmin, tmax)

  ## ---- posizioni di background per rifiuto-accettazione dalla
  ##      superficie GAM completa (intercetta + smooth(x,y) + distmin) ----
  ## IMPORTANTE: il valore di distmin usato per calcolare la probabilita'
  ## di accettazione viene CONSERVATO e riusato come valore finale del
  ## punto accettato (non ridisegnato in seguito) -- con rumore casuale
  ## nell'interpolatore, un secondo campionamento indipendente
  ## produrrebbe un valore diverso da quello che ha determinato
  ## l'accettazione, rompendo la coerenza interna dell'algoritmo.
  sample_bg_xy <- function(n_need) {
    out_x <- numeric(0); out_y <- numeric(0); out_dmin <- numeric(0)
    batch <- max(1000L, ceiling(3 * n_need))
    guard <- 0L
    while (length(out_x) < n_need) {
      guard <- guard + 1L
      if (guard > 20000L) stop("Rifiuto-accettazione per il background non converge.")
      xcand <- runif(batch, xmin, xmax); ycand <- runif(batch, ymin, ymax)
      dmin_cand <- ab_predict_distmin(xcand, ycand)
      etac <- as.numeric(predict(mod, newdata = data.frame(x = xcand, y = ycand, distmin = dmin_cand), type = "link"))
      acc <- runif(batch) <= exp(etac - eta_max)
      if (any(acc)) {
        out_x <- c(out_x, xcand[acc]); out_y <- c(out_y, ycand[acc]); out_dmin <- c(out_dmin, dmin_cand[acc])
      }
      if (length(out_x) < n_need && guard %% 20L == 0L) batch <- min(batch * 2L, 200000L)
    }
    data.frame(x = out_x[seq_len(n_need)], y = out_y[seq_len(n_need)], distmin = out_dmin[seq_len(n_need)])
  }
  xy0 <- sample_bg_xy(n0)
  x0 <- xy0$x; y0 <- xy0$y
  distmin0 <- xy0$distmin

  m0_bg <- m0 + rexp(n0, rate = beta_gr)
  long0 <- x0 / unit; lat0 <- y0 / unit

  cat.new <- data.frame(
    event_id = seq_len(n0), father_id = 0L, lgen = 0L,
    time = t0, lat = lat0, long = long0, z = 0, magn1 = m0_bg,
    x_km = x0, y_km = y0, m_rel = m0_bg - m0, distmin = distmin0,
    stringsAsFactors = FALSE
  )
  cat.new <- cat.new[order(cat.new$time), , drop = FALSE]
  rownames(cat.new) <- NULL
  cat.new$event_id <- seq_len(nrow(cat.new))

  ## ---- branching/triggering: kernel Omori-Utsu temporale + kernel
  ##      isotropo spaziale, con la covariata distmin nel predittore di
  ##      produttivita' (formula1 = "time ~ magnitude + distmin - 1"):
  ##      intensita' attesa di figli = k0 * exp(alpha*m_rel + beta_distmin*distmin) * (integrale Omori) * (integrale kernel spaziale) ----
  ak <- theta_hat$k0_hat * theta_hat$c_hat^(1 - theta_hat$p_hat) / (theta_hat$p_hat - 1)
  sk <- (pi * theta_hat$d_hat^(1 - theta_hat$q_hat)) / (theta_hat$q_hat - 1)

  sample_xy_etas <- function(n, x0c, y0c, d, q) {
    theta <- runif(n, 0, 2 * pi)
    U <- runif(n)
    R <- sqrt(d * (U^(1 / (1 - q)) - 1))
    cbind(x0c + R * cos(theta), y0c + R * sin(theta))
  }

  exploded <- FALSE
  i <- 0L
  while (i < nrow(cat.new)) {
    i <- i + 1L
    n_exp_i <- ak * sk * exp(theta_hat$betacov_hat * cat.new$m_rel[i] + theta_hat$beta_distmin_hat * cat.new$distmin[i])
    if (!is.finite(n_exp_i) || n_exp_i < 0) n_exp_i <- 0
    ni <- rpois(1, n_exp_i)

    if (ni > 0L) {
      t_child <- theta_hat$c_hat * runif(ni)^(-1 / (theta_hat$p_hat - 1)) - theta_hat$c_hat + cat.new$time[i]
      xy_child <- sample_xy_etas(ni, cat.new$x_km[i], cat.new$y_km[i], theta_hat$d_hat, theta_hat$q_hat)
      inside <- (t_child > tmin) & (t_child < tmax) &
        (xy_child[, 1] > xmin) & (xy_child[, 1] < xmax) &
        (xy_child[, 2] > ymin) & (xy_child[, 2] < ymax)

      if (any(inside)) {
        nt <- sum(inside)
        x1 <- xy_child[inside, 1]; y1 <- xy_child[inside, 2]; t1 <- t_child[inside]
        m1 <- m0 + rexp(nt, rate = beta_gr)
        distmin1 <- ab_predict_distmin(x1, y1)
        long1 <- x1 / unit; lat1 <- y1 / unit

        child_df <- data.frame(
          event_id = NA_integer_, father_id = as.integer(cat.new$event_id[i]),
          lgen = as.integer(cat.new$lgen[i] + 1L),
          time = t1, lat = lat1, long = long1, z = 0, magn1 = m1,
          x_km = x1, y_km = y1, m_rel = m1 - m0, distmin = distmin1,
          stringsAsFactors = FALSE
        )
        cat.new <- rbind(cat.new, child_df)
      }
    }
    if (nrow(cat.new) > max_events) { exploded <- TRUE; break }
  }

  cat.new <- cat.new[order(cat.new$time), , drop = FALSE]
  rownames(cat.new) <- NULL
  cat.new$event_id <- seq_len(nrow(cat.new))

  list(cat.sim = cat.new, n0 = n0, nson = nrow(cat.new) - n0, exploded = exploded)
}


## ================================================================
## 6. UNA REPLICA BOOTSTRAP COMPLETA: SIMULA + RIFITTA L'INTERA ETAS-P
## ================================================================

run_one_application_bootstrap <- function(b, theta_hat, seed_base = AB_SEED_BASE, out_dir = AB_OUT_DIR) {
  out_file <- file.path(out_dir, "boot_fits", sprintf("boot_%04d.rds", b))
  if (file.exists(out_file)) return(invisible(out_file))   # resumable: salta se gia' fatto

  seed_b <- seed_base + b
  sim_res <- ab_safe_run(function() {
    simulate_one_application_bootstrap(
      theta_hat, AB_TMIN, AB_TMAX, AB_LONG_RANGE, AB_LAT_RANGE,
      AB_M0, AB_BETA_GR, seed = seed_b
    )
  })

  fit_res <- NULL
  if (isTRUE(sim_res$ok) && nrow(sim_res$value$cat.sim) > 10 && !isTRUE(sim_res$value$exploded)) {
    cat_boot <- sim_res$value$cat.sim
    n <- nrow(cat_boot)   # a M0=2.5 filtrato==non filtrato, quindi w/hvarx/hvary hanno tutti la stessa lunghezza n

    starts <- list(
      mu = theta_hat$mu_hat, k0 = theta_hat$k0_hat, c = theta_hat$c_hat, p = theta_hat$p_hat,
      gamma = 0, d = theta_hat$d_hat, q = theta_hat$q_hat, betacov = theta_hat$betacov_hat
    )

    fit_res <- ab_safe_run(function() {
      etasclass.par(
        cat.orig = cat_boot, time.update = FALSE, magn.threshold = AB_M0, magn.threshold.back = 3.9,
        tmax = max(cat_boot$time), long.range = AB_LONG_RANGE, lat.range = AB_LAT_RANGE,
        mu = starts$mu, k0 = starts$k0, c = starts$c, p = starts$p, gamma = 0,
        d = starts$d, q = starts$q, betacov = starts$betacov,
        params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
        formula1 = "time ~ magnitude + distmin - 1",
        offset = 0, hdef = c(1, 1), w = replicate(n, 1), hvarx = replicate(n, 1), hvary = replicate(n, 1),
        declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
        ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE,
        description = "", cat.back = NULL, back.smooth = 1,
        sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS",
        compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
        formula.bg = ~ s(x, y, k = 25) + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
        type.cov.values.bg = list(distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
        verbose.bg = FALSE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE,
        mark.c.bg = FALSE, seed.bg = 2
      )
    })
  }

  out_obj <- list(
    b = b, seed = seed_b,
    n_events_sim = if (isTRUE(sim_res$ok)) nrow(sim_res$value$cat.sim) else NA_integer_,
    exploded = if (isTRUE(sim_res$ok)) sim_res$value$exploded else NA,
    sim_ok = isTRUE(sim_res$ok), sim_error = sim_res$error_message,
    fit_ok = if (!is.null(fit_res)) isTRUE(fit_res$ok) else FALSE,
    fit_error = if (!is.null(fit_res)) fit_res$error_message else "simulation_or_size_guard_failed",
    theta0_hat = if (!is.null(fit_res) && isTRUE(fit_res$ok)) tryCatch(coef(fit_res$value$model.bg$mod_global)["(Intercept)"], error = function(e) NA) else NA,
    theta_distmin_hat = if (!is.null(fit_res) && isTRUE(fit_res$ok)) tryCatch(coef(fit_res$value$model.bg$mod_global)["distmin"], error = function(e) NA) else NA,
    beta_distmin_trig_hat = if (!is.null(fit_res) && isTRUE(fit_res$ok)) tryCatch(as.numeric(fit_res$value$params.MLtot["distmin"]), error = function(e) NA) else NA,
    params_MLtot = if (!is.null(fit_res) && isTRUE(fit_res$ok)) fit_res$value$params.MLtot else NULL
  )
  ab_saveRDS_atomic(out_obj, out_file)
  cat(sprintf("[%s] boot %04d: sim_ok=%s fit_ok=%s n_sim=%s\n",
             format(Sys.time(), "%H:%M:%S"), b, out_obj$sim_ok, out_obj$fit_ok, out_obj$n_events_sim %||% "NA"),
     file = file.path(out_dir, "logs", "progress_appboot.log"), append = TRUE)
  invisible(out_file)
}


## ================================================================
## 7. ORCHESTRAZIONE PARALLELA (stesso pattern gia' usato altrove)
## ================================================================

ab_run_parallel <- function(job_ids, worker_fun, n_cores = AB_N_CORES) {
  if (n_cores <= 1L || length(job_ids) <= 1L) return(lapply(job_ids, worker_fun))
  if (.Platform$OS.type != "windows") {
    parallel::mclapply(job_ids, worker_fun, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, { library(mgcv) })
    parallel::parLapply(cl, job_ids, worker_fun)
  }
}

run_application_bootstrap_plan <- function(B = AB_B, n_cores = AB_N_CORES, theta_hat = ab_theta_hat) {
  todo <- seq_len(B)
  already <- vapply(todo, function(b) file.exists(file.path(AB_OUT_DIR, "boot_fits", sprintf("boot_%04d.rds", b))), logical(1))
  ab_log("Repliche gia' presenti: %d/%d. Da eseguire: %d.", sum(already), B, sum(!already))
  ab_run_parallel(todo, function(b) run_one_application_bootstrap(b, theta_hat), n_cores = n_cores)
  invisible(TRUE)
}


## ================================================================
## 8. RACCOLTA, INTERVALLO PERCENTILE, TABELLA RIASSUNTIVA
## ================================================================

collect_application_bootstrap <- function(out_dir = AB_OUT_DIR) {
  files <- list.files(file.path(out_dir, "boot_fits"), pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    o <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(o)) return(NULL)
    p <- o$params_MLtot
    data.frame(
      b = o$b, sim_ok = o$sim_ok, fit_ok = o$fit_ok, n_events_sim = o$n_events_sim,
      theta0_hat = as.numeric(o$theta0_hat), theta_distmin_hat = as.numeric(o$theta_distmin_hat),
      beta_distmin_trig_hat = as.numeric(o$beta_distmin_trig_hat),
      mu_hat = if (!is.null(p)) as.numeric(p["mu"]) else NA_real_,
      k0_hat = if (!is.null(p)) as.numeric(p["k0"]) else NA_real_,
      c_hat = if (!is.null(p)) as.numeric(p["c"]) else NA_real_,
      d_hat = if (!is.null(p)) as.numeric(p["d"]) else NA_real_,
      q_hat = if (!is.null(p)) as.numeric(p["q"]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

## Filtro di plausibilita' sulle repliche bootstrap INTERNE: stessa
## logica gia' usata nello studio di copertura sulla simulazione
## (S2.6), qui applicata a ciascuno dei B rifit invece che a repliche
## esterne del main plan. Scarta repliche con k0/c/d/q stimato troppo
## lontano (>5x in entrambe le direzioni) dal valore usato per
## generarle -- segno di ottimizzazione non convergente/degenere su
## quel particolare catalogo sintetico, osservato empiricamente (es.
## replica con beta_distmin_trig_hat di segno opposto e ampiezza >2x).
ab_flag_plausible <- function(df, theta_hat = ab_theta_hat, max_ratio = 5) {
  ratio_ok <- function(est, true_val) {
    if (!is.finite(est) || !is.finite(true_val) || true_val == 0) return(TRUE)
    r <- est / true_val
    is.finite(r) && r >= 1 / max_ratio && r <= max_ratio
  }
  df$plausible <- mapply(function(k0, c, d, q) {
    ratio_ok(k0, theta_hat$k0_hat) && ratio_ok(c, theta_hat$c_hat) &&
      ratio_ok(d, theta_hat$d_hat) && ratio_ok(q, theta_hat$q_hat)
  }, df$k0_hat, df$c_hat, df$d_hat, df$q_hat)
  df
}

summarise_application_bootstrap <- function(df = collect_application_bootstrap(), theta_hat = ab_theta_hat, mad_k = 5) {
  ok <- df[isTRUE(df$fit_ok) | df$fit_ok == TRUE, ]
  ok <- ok[is.finite(ok$theta_distmin_hat), ]
  ok <- ab_flag_plausible(ok, theta_hat)   # filtro 1: k0/c/d/q plausibili (cattura fit degeneri del TRIGGERING)

  n_ok <- nrow(ok)
  n_implausible_trig <- sum(!ok$plausible, na.rm = TRUE)
  ok_stage1 <- ok[ok$plausible, ]

  ## filtro 2: outlier robusti (MAD) sul singolo coefficiente di
  ## interesse -- cattura fit degeneri specifici del passo di
  ## regressione del BACKGROUND, che il filtro 1 non vede (i due passi
  ## dell'algoritmo alternato possono fallire indipendentemente).
  mad_outlier_flag <- function(x, k = mad_k) {
    x_fin <- x[is.finite(x)]
    if (length(x_fin) < 10) return(rep(TRUE, length(x)))
    med <- median(x_fin); mad_x <- mad(x_fin)
    if (!is.finite(mad_x) || mad_x == 0) return(rep(TRUE, length(x)))
    abs(x - med) <= k * mad_x
  }

  ab_log("Repliche con fit riuscito: %d/%d. Flaggate implausibili (k0/c/d/q): %d.", n_ok, nrow(df), n_implausible_trig)

  fit0 <- theta_hat$fit_param
  point_theta_distmin <- tryCatch(coef(fit0$model.bg$mod_global)["distmin"], error = function(e) NA)
  point_beta_distmin <- as.numeric(fit0$params.MLtot["distmin"])

  summarise_one <- function(vals, label, point_estimate) {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 10) return(data.frame(coefficient = label, point_estimate = point_estimate,
                                             n_ok = length(vals), n_mad_excluded = NA, boot_se = NA, ci_lo = NA, ci_hi = NA))
    keep <- mad_outlier_flag(vals)
    vals_robust <- vals[keep]
    qs <- quantile(vals_robust, probs = c(0.025, 0.975))
    data.frame(coefficient = label, point_estimate = point_estimate,
              n_ok = length(vals_robust), n_mad_excluded = sum(!keep),
              boot_se = sd(vals_robust), ci_lo = unname(qs[1]), ci_hi = unname(qs[2]))
  }

  tab_all <- rbind(
    summarise_one(ok$theta_distmin_hat, "theta_distmin (background) -- tutte le repliche riuscite", point_theta_distmin),
    summarise_one(ok$beta_distmin_trig_hat, "beta_distmin (triggering) -- tutte le repliche riuscite", point_beta_distmin)
  )
  tab_stable <- rbind(
    summarise_one(ok_stage1$theta_distmin_hat, "theta_distmin (background) -- k0/c/d/q plausibili + MAD", point_theta_distmin),
    summarise_one(ok_stage1$beta_distmin_trig_hat, "beta_distmin (triggering) -- k0/c/d/q plausibili + MAD", point_beta_distmin)
  )
  list(all_replicates = tab_all, stable_replicates = tab_stable,
      n_ok = n_ok, n_implausible = n_implausible_trig,
      raw_values = ok_stage1[, c("b", "theta_distmin_hat", "beta_distmin_trig_hat")])
}


## ================================================================
## 9. AUTORUN SWITCH
## ================================================================

if (identical(AB_MODE, "test")) {
  ab_log("APPBOOT_MODE='test': 2 repliche seriali per verificare che tutto funzioni.")
  run_application_bootstrap_plan(B = 2L, n_cores = 1L)
  print(collect_application_bootstrap())

} else if (identical(AB_MODE, "run")) {
  ab_log("APPBOOT_MODE='run': B=%d, n_cores=%d.", AB_B, AB_N_CORES)
  run_application_bootstrap_plan()
  res <- summarise_application_bootstrap()
  utils::write.csv(res$all_replicates, file.path(AB_OUT_DIR, "table_application_bootstrap_all.csv"), row.names = FALSE)
  utils::write.csv(res$stable_replicates, file.path(AB_OUT_DIR, "table_application_bootstrap_stable.csv"), row.names = FALSE)
  cat("\n--- Tutte le repliche riuscite ---\n"); print(res$all_replicates, row.names = FALSE, digits = 4)
  cat("\n--- Solo repliche numericamente stabili (n_implausible =", res$n_implausible, "su", res$n_ok, ") ---\n")
  print(res$stable_replicates, row.names = FALSE, digits = 4)

} else if (identical(AB_MODE, "collect")) {
  res <- summarise_application_bootstrap()
  utils::write.csv(res$all_replicates, file.path(AB_OUT_DIR, "table_application_bootstrap_all.csv"), row.names = FALSE)
  utils::write.csv(res$stable_replicates, file.path(AB_OUT_DIR, "table_application_bootstrap_stable.csv"), row.names = FALSE)
  cat("\n--- Tutte le repliche riuscite ---\n"); print(res$all_replicates, row.names = FALSE, digits = 4)
  cat("\n--- Solo repliche numericamente stabili (n_implausible =", res$n_implausible, "su", res$n_ok, ") ---\n")
  print(res$stable_replicates, row.names = FALSE, digits = 4)

} else {
  ab_log("APPBOOT_MODE='%s': funzioni caricate, nessun calcolo lanciato.", AB_MODE)
}
