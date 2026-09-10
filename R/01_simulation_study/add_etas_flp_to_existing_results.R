#######################################################################
# NO-OVERWRITE VERSION
#
# This script NEVER modifies the original output directory by default.
# It reads existing .rds files from:
#   ETAS_SOURCE_OUT_ROOT
# default:
#   <getwd()>/etas_parametric_final_no_smooth_outputs_v1
#
# It copies only the .rds result files into a separate target root:
#   ETAS_FLP_TARGET_ROOT
# default:
#   <ETAS_SOURCE_OUT_ROOT>_with_flp
#
# Then it adds obj$fits$flp only inside the copied .rds files in the
# target root. The original .rds files are left untouched.
#
# Recommended use:
#   setwd("/home/nicolettadangelo/sim_paper1_marco/etas parametrico")
#
#   Sys.setenv(ETAS_FLP_MODE = "check")
#   source("add_etas_flp_to_existing_results_NO_OVERWRITE.R")
#
#   Sys.setenv(ETAS_FLP_MODE = "test")
#   source("add_etas_flp_to_existing_results_NO_OVERWRITE.R")
#
#   Sys.setenv(ETAS_FLP_MODE = "run")
#   source("add_etas_flp_to_existing_results_NO_OVERWRITE.R")
#
# Useful optional variables:
#   ETAS_SOURCE_OUT_ROOT     original output root to read from
#   ETAS_FLP_TARGET_ROOT     new output root to write to
#   ETAS_FLP_REFRESH_TARGET  TRUE/FALSE; if TRUE, recopies source .rds over target .rds
#
# IMPORTANT: the script still fits only ETAS-FLP. It does not refit
# obj$fits$classic, obj$fits$param_correct, or obj$fits$param_misspec.
#######################################################################

flp_no_overwrite_bool <- function(name, default = FALSE) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("1", "true", "t", "yes", "y", "si", "sì")
}

ETAS_SOURCE_OUT_ROOT_NO_OVERWRITE <- Sys.getenv(
  "ETAS_SOURCE_OUT_ROOT",
  unset = file.path(getwd(), "etas_parametric_final_no_smooth_outputs_v1")
)

ETAS_FLP_TARGET_ROOT_NO_OVERWRITE <- Sys.getenv(
  "ETAS_FLP_TARGET_ROOT",
  unset = paste0(ETAS_SOURCE_OUT_ROOT_NO_OVERWRITE, "_with_flp")
)

## Force the downstream post-processing code to work on the copied target root only.
Sys.setenv(ETAS_OUT_ROOT = ETAS_FLP_TARGET_ROOT_NO_OVERWRITE)

flp_prepare_no_overwrite_target <- function(source_root = ETAS_SOURCE_OUT_ROOT_NO_OVERWRITE,
                                            target_root = ETAS_FLP_TARGET_ROOT_NO_OVERWRITE,
                                            blocks = c(
                                              "01_main_plan",
                                              "02_secondary_sample_size",
                                              "03_secondary_misspecification",
                                              "04_secondary_extreme_triggering"
                                            ),
                                            refresh_target = flp_no_overwrite_bool("ETAS_FLP_REFRESH_TARGET", FALSE)) {
  if (!dir.exists(source_root)) {
    stop("Original source output root not found: ", source_root)
  }
  if (normalizePath(source_root, mustWork = TRUE) == normalizePath(target_root, mustWork = FALSE)) {
    stop("Refusing to run: source_root and target_root are the same. Choose a different ETAS_FLP_TARGET_ROOT.")
  }

  if (!dir.exists(target_root)) dir.create(target_root, recursive = TRUE, showWarnings = FALSE)

  manifest <- data.frame()

  for (block in blocks) {
    src_fits <- file.path(source_root, block, "fits")
    tgt_fits <- file.path(target_root, block, "fits")

    if (!dir.exists(src_fits)) {
      warning("Source fits directory not found and will be skipped: ", src_fits)
      next
    }

    if (!dir.exists(tgt_fits)) dir.create(tgt_fits, recursive = TRUE, showWarnings = FALSE)

    src_files <- list.files(src_fits, pattern = "\\.rds$", full.names = TRUE)
    if (!length(src_files)) next

    for (src_file in src_files) {
      tgt_file <- file.path(tgt_fits, basename(src_file))
      action <- "kept_existing_target"
      copied <- FALSE

      if (refresh_target || !file.exists(tgt_file)) {
        copied <- file.copy(src_file, tgt_file, overwrite = refresh_target)
        action <- if (copied) {
          if (refresh_target) "copied_refreshed_from_source" else "copied_from_source"
        } else {
          "copy_failed"
        }
      }

      manifest <- rbind(
        manifest,
        data.frame(
          block = block,
          source_file = src_file,
          target_file = tgt_file,
          action = action,
          copied = copied,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  manifest_file <- file.path(target_root, "flp_no_overwrite_copy_manifest.csv")
  utils::write.csv(manifest, manifest_file, row.names = FALSE)

  cat("[NO-OVERWRITE] Source root: ", source_root, "\n", sep = "")
  cat("[NO-OVERWRITE] Target root: ", target_root, "\n", sep = "")
  cat("[NO-OVERWRITE] Manifest: ", manifest_file, "\n", sep = "")
  cat("[NO-OVERWRITE] Source .rds files copied where needed: ", sum(manifest$action %in% c("copied_from_source", "copied_refreshed_from_source")), "\n", sep = "")
  cat("[NO-OVERWRITE] Existing target .rds files kept: ", sum(manifest$action == "kept_existing_target"), "\n", sep = "")

  invisible(manifest)
}

flp_prepare_no_overwrite_target()

#######################################################################
# ADD ETAS-FLP FITS TO EXISTING PARAMETRIC-ETAS SIMULATION OUTPUTS
#
# Purpose
# -------
# This script is a POST-PROCESSING script. It reads the existing .rds
# files produced by etas_parametric_final_plan_no_smooth.R, reconstructs
# the exact same simulated catalogues already stored in each object, fits
# ETAS with FLP bandwidth selection, and saves the new fit as:
#
#   obj$fits$flp
#
# It then recomputes obj$payload so that metric tables can compare:
#   - obj$fits$classic        : ETAS with kernel background, Silverman h
#   - obj$fits$param_correct  : proposed parametric-background ETAS
#   - obj$fits$flp            : ETAS with kernel background, FLP-selected h
#   - obj$fits$param_misspec  : optional misspecified parametric model
#
# IMPORTANT
# ---------
# Do NOT source the original simulation file directly for this task,
# because the current version contains a forced RUN_MODE = "main" near
# the end of the file. This script safely loads only the function and
# constant definitions before the AUTORUN SWITCH section.
#
# Recommended use from your working directory:
#   setwd("/home/nicolettadangelo/sim_paper1_marco/etas parametrico")
#
#   ## 1) First check existing files, without fitting anything
#   Sys.setenv(ETAS_FLP_MODE = "check")
#   source("add_etas_flp_to_existing_results.R")
#
#   ## 2) Test on one catalogue only
#   Sys.setenv(ETAS_FLP_MODE = "test")
#   source("add_etas_flp_to_existing_results.R")
#
#   ## 3) Run all remaining FLP fits sequentially
#   Sys.setenv(ETAS_FLP_MODE = "run")
#   source("add_etas_flp_to_existing_results.R")
#
# Useful optional environment variables:
#   ETAS_BASE_DIR      default: getwd()
#   ETAS_PLAN_SCRIPT   default: etas_parametric_final_plan_no_smooth.R in BASE_DIR
#   ETAS_OUT_ROOT      default: BASE_DIR/etas_parametric_final_no_smooth_outputs_v1
#   ETAS_FLP_BLOCKS    comma-separated block names; default: all four blocks
#   ETAS_FLP_FORCE     TRUE/FALSE; refit FLP even if obj$fits$flp already exists
#   ETAS_FLP_MAX_FILES integer; useful for testing
#   ETAS_FLP_BACKUP    TRUE/FALSE; copy each .rds before modifying it
#   ETAS_FLP_CAPTURE_OUTPUT TRUE/FALSE; write verbose etasclass output to logs_flp/
#######################################################################

## ================================================================
## 0. Small independent utilities
## ================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

flp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

flp_log <- function(..., .level = "INFO") {
  cat(sprintf("[%s] [%s] ", flp_now(), .level), sprintf(...), "\n", sep = "")
  flush.console()
}

flp_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

flp_env_bool <- function(name, default = FALSE) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("1", "true", "t", "yes", "y", "si", "sì")
}

flp_env_int <- function(name, default = NA_integer_) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) default else out
}

flp_env_chr <- function(name, default = NULL) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) default else x
}

flp_split_env <- function(name, default) {
  x <- flp_env_chr(name, NULL)
  if (is.null(x)) return(default)
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

flp_true <- function(x) {
  !is.na(x) & as.logical(x)
}

flp_append_csv_row <- function(row, file) {
  flp_dir_create(dirname(file))
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  utils::write.table(
    row,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file)
  )
  invisible(file)
}

flp_saveRDS_atomic <- function(object, file, compress = "gzip") {
  flp_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_flp_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) {
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(file)
}

flp_safe_run <- function(expr_fun) {
  warnings <- character(0)
  t0 <- Sys.time()
  res <- withCallingHandlers(
    tryCatch(
      expr_fun(),
      error = function(e) structure(list(error = conditionMessage(e)), class = "flp_safe_error")
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  t1 <- Sys.time()
  if (inherits(res, "flp_safe_error")) {
    list(
      ok = FALSE,
      value = NULL,
      error_message = res$error,
      warnings = warnings,
      elapsed_sec = as.numeric(difftime(t1, t0, units = "secs"))
    )
  } else {
    list(
      ok = TRUE,
      value = res,
      error_message = NA_character_,
      warnings = warnings,
      elapsed_sec = as.numeric(difftime(t1, t0, units = "secs"))
    )
  }
}

flp_rbind_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && is.data.frame(z) && nrow(z) > 0, x)
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(df) {
    miss <- setdiff(all_names, names(df))
    for (m in miss) df[[m]] <- NA
    df <- df[, all_names, drop = FALSE]
    rownames(df) <- NULL
    df
  })
  do.call(rbind, x)
}

flp_pluck <- function(x, name, default = NULL) {
  if (is.list(x) && name %in% names(x) && !is.null(x[[name]])) x[[name]] else default
}

flp_as_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(default)
  x <- x[1]
  if (is.null(x) || is.na(x)) default else x
}

## ================================================================
## 1. Configuration
## ================================================================

FLP_BASE_DIR <- flp_env_chr("ETAS_BASE_DIR", getwd())

flp_guess_plan_script <- function(base_dir) {
  explicit <- flp_env_chr("ETAS_PLAN_SCRIPT", NULL)
  candidates <- c(
    explicit,
    file.path(base_dir, "etas_parametric_final_plan_no_smooth.R"),
    file.path(base_dir, "etas_parametric_final_plan_no_smooth(1).R"),
    file.path(base_dir, "etas_parametric_final_plan_no_smooth_no_smooth.R")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop(
      "Could not find the original simulation-plan script. Tried:\n  ",
      paste(candidates, collapse = "\n  "),
      "\nSet Sys.setenv(ETAS_PLAN_SCRIPT = '/full/path/to/etas_parametric_final_plan_no_smooth.R')."
    )
  }
  normalizePath(hit[1], mustWork = TRUE)
}

FLP_PLAN_SCRIPT <- flp_guess_plan_script(FLP_BASE_DIR)
FLP_OUT_ROOT <- flp_env_chr(
  "ETAS_OUT_ROOT",
  file.path(FLP_BASE_DIR, "etas_parametric_final_no_smooth_outputs_v1")
)

FLP_DEFAULT_BLOCKS <- c(
  "01_main_plan",
  "02_secondary_sample_size",
  "03_secondary_misspecification",
  "04_secondary_extreme_triggering"
)

FLP_BLOCKS <- flp_split_env("ETAS_FLP_BLOCKS", FLP_DEFAULT_BLOCKS)
FLP_MODE <- tolower(flp_env_chr("ETAS_FLP_MODE", "check"))
FLP_FORCE <- flp_env_bool("ETAS_FLP_FORCE", FALSE)
FLP_BACKUP <- flp_env_bool("ETAS_FLP_BACKUP", FALSE)
FLP_CAPTURE_OUTPUT <- flp_env_bool("ETAS_FLP_CAPTURE_OUTPUT", TRUE)
FLP_MAX_FILES <- flp_env_int("ETAS_FLP_MAX_FILES", NA_integer_)

## Optional fitting controls. Defaults intentionally match the final plan.
FLP_NDECLUST <- flp_env_int("ETAS_FLP_NDECLUST", NA_integer_)
FLP_ITERLIM <- flp_env_int("ETAS_FLP_ITERLIM", NA_integer_)
FLP_NITERWEIGHT <- flp_env_int("ETAS_FLP_NITERWEIGHT", 1L)
FLP_NTHETA <- flp_env_int("ETAS_FLP_NTHETA", 36L)
FLP_EPSMAX <- as.numeric(flp_env_chr("ETAS_FLP_EPSMAX", "1e-4"))

## ================================================================
## 2. Safely load the original simulation functions, without autorun
## ================================================================

load_etas_plan_definitions_no_autorun <- function(plan_script = FLP_PLAN_SCRIPT,
                                                  envir = .GlobalEnv) {
  if (!file.exists(plan_script)) stop("plan_script not found: ", plan_script)

  flp_log("Loading definitions from: %s", plan_script)
  txt <- readLines(plan_script, warn = FALSE)

  idx <- grep("^##[[:space:]]*11\\.[[:space:]]*AUTORUN SWITCH", txt)
  if (!length(idx)) {
    idx <- grep("^print_plan_sizes\\(\\)", txt)
  }
  if (!length(idx)) {
    stop(
      "Could not locate the AUTORUN SWITCH in the original script. ",
      "I refuse to source the whole file because it may launch simulations."
    )
  }

  keep <- seq_len(idx[1] - 1L)
  code <- paste(txt[keep], collapse = "\n")
  eval(parse(text = code), envir = envir)

  ## The copied/local etasclass() in the simulation script calls an
  ## internal etasFLP helper, flp1.etas.nlmNEW(). In the original package
  ## this helper is visible from the etasFLP namespace, but after sourcing
  ## a local copy of etasclass() its enclosing environment is no longer the
  ## package namespace. Therefore we explicitly import the helper into the
  ## working environment. This is needed only for flp = TRUE; the previous
  ## simulations used flp = FALSE and therefore never triggered this call.
  if (!requireNamespace("etasFLP", quietly = TRUE)) {
    stop("Package 'etasFLP' is required for ETAS-FLP fitting but is not installed/available.")
  }

  if (!exists("flp1.etas.nlmNEW", envir = envir, inherits = FALSE)) {
    flp_fun <- tryCatch(
      getFromNamespace("flp1.etas.nlmNEW", "etasFLP"),
      error = function(e) NULL
    )

    if (!is.function(flp_fun)) {
      ns_flp <- tryCatch(
        ls(getNamespace("etasFLP"), pattern = "flp", all.names = TRUE),
        error = function(e) character(0)
      )
      stop(
        "Could not import internal helper 'flp1.etas.nlmNEW' from package etasFLP. ",
        "Available namespace objects matching 'flp': ",
        paste(ns_flp, collapse = ", "),
        ". Check the installed etasFLP version."
      )
    }

    assign("flp1.etas.nlmNEW", flp_fun, envir = envir)
    flp_log("Imported internal etasFLP helper: flp1.etas.nlmNEW")
  }

  required <- c(
    "etasclass", "flp1.etas.nlmNEW", "make_eqcat_from_sim", "get_cat_sim", "make_metric_payload",
    "make_grid_eval", "TRUE_BASE_FULL", "TRUE_BETACOV_FULL", "M0_FULL",
    "LONG_RANGE", "LAT_RANGE", "NDECLUST_FIT_NEW", "ITERLIM_FIT_NEW"
  )
  missing_required <- required[!vapply(required, exists, logical(1), envir = envir, inherits = TRUE)]
  if (length(missing_required)) {
    stop("Missing required objects after loading definitions: ", paste(missing_required, collapse = ", "))
  }

  flp_log("Definitions loaded. No simulation runner has been launched.")
  invisible(TRUE)
}

## ================================================================
## 3. ETAS-FLP fitting wrapper
## ================================================================

make_flp_starts_from_obj <- function(obj) {
  cal <- flp_pluck(obj, "calibration", list())

  list(
    mu = as.numeric(flp_pluck(cal, "mu_cal", 0.4)),
    k0 = as.numeric(flp_pluck(cal, "k0_cal", 0.01)),
    c = as.numeric(TRUE_BASE_FULL["c"]),
    p = as.numeric(TRUE_BASE_FULL["p"]),
    gamma = 0,
    d = as.numeric(TRUE_BASE_FULL["d"]),
    q = as.numeric(TRUE_BASE_FULL["q"]),
    betacov = as.numeric(TRUE_BETACOV_FULL)
  )
}

fit_etas_flp_full <- function(cat_sim_df,
                              starts = NULL,
                              ndeclust = if (is.na(FLP_NDECLUST)) NDECLUST_FIT_NEW else FLP_NDECLUST,
                              iterlim = if (is.na(FLP_ITERLIM)) ITERLIM_FIT_NEW else FLP_ITERLIM,
                              n.iterweight = FLP_NITERWEIGHT,
                              epsmax = FLP_EPSMAX,
                              ntheta = FLP_NTHETA,
                              m1 = NULL) {
  if (is.null(starts)) {
    starts <- list(
      mu = 0.4,
      k0 = 0.01,
      c = as.numeric(TRUE_BASE_FULL["c"]),
      p = as.numeric(TRUE_BASE_FULL["p"]),
      gamma = 0,
      d = as.numeric(TRUE_BASE_FULL["d"]),
      q = as.numeric(TRUE_BASE_FULL["q"]),
      betacov = as.numeric(TRUE_BETACOV_FULL)
    )
  }

  etasclass(
    cat.orig = cat_sim_df,
    magn.threshold = M0_FULL,
    magn.threshold.back = M0_FULL,
    tmax = max(cat_sim_df$time),
    long.range = LONG_RANGE,
    lat.range = LAT_RANGE,
    mu = starts$mu,
    k0 = starts$k0,
    c = starts$c,
    p = starts$p,
    gamma = 0,
    d = starts$d,
    q = max(starts$q, 1.01),
    betacov = starts$betacov,
    params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    w = rep(1, nrow(cat_sim_df)),
    hvarx = rep(1, nrow(cat_sim_df)),
    hvary = rep(1, nrow(cat_sim_df)),
    formula1 = "time ~ magnitude - 1",
    declustering = TRUE,
    thinning = FALSE,
    flp = TRUE,
    m1 = m1,
    ndeclust = ndeclust,
    n.iterweight = n.iterweight,
    onlytime = FALSE,
    is.backconstant = FALSE,
    description = "ETAS-FLP added as post-processing fit",
    sectoday = FALSE,
    longlat.to.km = TRUE,
    usenlm = TRUE,
    compsqm = TRUE,
    epsmax = epsmax,
    iterlim = iterlim,
    ntheta = ntheta
  )
}

## ================================================================
## 4. Status and validation helpers
## ================================================================

flp_block_from_file <- function(file) {
  ## Expected: OUT_ROOT/<block>/fits/<file>.rds
  basename(dirname(dirname(file)))
}

flp_status_one_file <- function(file, obj = NULL) {
  if (is.null(obj)) {
    obj <- tryCatch(readRDS(file), error = function(e) {
      structure(list(.read_error = conditionMessage(e)), class = "flp_read_error")
    })
  }

  if (inherits(obj, "flp_read_error")) {
    return(data.frame(
      file = file,
      block = flp_block_from_file(file),
      bg_case = NA_character_,
      count_case = NA_character_,
      rep = NA_integer_,
      read_ok = FALSE,
      simulation_ok = NA,
      nobs = NA_integer_,
      has_classic = NA,
      classic_ok = NA,
      has_param_correct = NA,
      param_correct_ok = NA,
      has_param_misspec = NA,
      param_misspec_ok = NA,
      has_flp = NA,
      flp_ok = NA,
      payload_has_flp = NA,
      read_error = obj$.read_error,
      stringsAsFactors = FALSE
    ))
  }

  meta <- flp_pluck(obj, "metadata", list())
  fits <- flp_pluck(obj, "fits", list())
  sim <- flp_pluck(obj, "simulation", NULL)
  payload <- flp_pluck(obj, "payload", list())

  sim_ok <- !is.null(sim) && isTRUE(sim$ok)
  nobs <- NA_integer_
  if (sim_ok) {
    nobs <- tryCatch(nrow(get_cat_sim(sim$value)), error = function(e) NA_integer_)
  }

  data.frame(
    file = file,
    block = flp_block_from_file(file),
    bg_case = as.character(flp_pluck(meta, "bg_case", NA_character_)),
    count_case = as.character(flp_pluck(meta, "count_case", NA_character_)),
    rep = as.integer(flp_pluck(meta, "rep", NA_integer_)),
    read_ok = TRUE,
    simulation_ok = sim_ok,
    nobs = as.integer(nobs),
    has_classic = !is.null(fits$classic),
    classic_ok = if (!is.null(fits$classic)) isTRUE(fits$classic$ok) else NA,
    has_param_correct = !is.null(fits$param_correct),
    param_correct_ok = if (!is.null(fits$param_correct)) isTRUE(fits$param_correct$ok) else NA,
    has_param_misspec = !is.null(fits$param_misspec),
    param_misspec_ok = if (!is.null(fits$param_misspec)) isTRUE(fits$param_misspec$ok) else NA,
    has_flp = !is.null(fits$flp),
    flp_ok = if (!is.null(fits$flp)) isTRUE(fits$flp$ok) else NA,
    payload_has_flp = !is.null(payload$fit_summaries) && "flp" %in% names(payload$fit_summaries),
    read_error = NA_character_,
    stringsAsFactors = FALSE
  )
}

flp_simulation_signature <- function(sim_value) {
  cat_df <- get_cat_sim(sim_value)
  data.frame(
    nobs = nrow(cat_df),
    min_time = suppressWarnings(min(cat_df$time, na.rm = TRUE)),
    max_time = suppressWarnings(max(cat_df$time, na.rm = TRUE)),
    sum_time = suppressWarnings(sum(cat_df$time, na.rm = TRUE)),
    sum_long = suppressWarnings(sum(cat_df$long, na.rm = TRUE)),
    sum_lat = suppressWarnings(sum(cat_df$lat, na.rm = TRUE)),
    sum_magn1 = suppressWarnings(sum(cat_df$magn1, na.rm = TRUE)),
    n_background_true = if ("father_id" %in% names(cat_df)) sum(cat_df$father_id == 0, na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

check_flp_fit <- function(fit_result, n_events) {
  if (is.null(fit_result)) {
    return(data.frame(check = "fit_result_exists", ok = FALSE, value = NA_character_))
  }

  if (!isTRUE(fit_result$ok)) {
    return(data.frame(
      check = c("safe_run_ok", "error_message"),
      ok = c(FALSE, FALSE),
      value = c("FALSE", as.character(fit_result$error_message %||% NA_character_)),
      stringsAsFactors = FALSE
    ))
  }

  fit <- fit_result$value
  nstep_flp <- fit[["nstep.flp"]] %||% NA_integer_
  nstep_kde <- fit[["nstep.kde"]] %||% NA_integer_
  hdef <- fit[["hdef"]] %||% NA_real_
  rho <- fit[["rho.weights"]] %||% NULL
  lambda_vec <- fit[["l"]] %||% NULL
  logl <- fit[["logl"]] %||% NA_real_

  data.frame(
    check = c(
      "safe_run_ok",
      "fit_flp_flag_TRUE",
      "nstep_flp_positive",
      "nstep_kde_positive",
      "hdef_length_2_finite",
      "rho_weights_length_matches_events",
      "lambda_length_matches_events",
      "loglik_finite"
    ),
    ok = c(
      TRUE,
      isTRUE(fit[["flp"]]),
      is.finite(nstep_flp) && nstep_flp > 0,
      is.finite(nstep_kde) && nstep_kde > 0,
      length(hdef) == 2 && all(is.finite(hdef)),
      !is.null(rho) && length(rho) == n_events && all(is.finite(rho)),
      !is.null(lambda_vec) && length(lambda_vec) == n_events && all(is.finite(lambda_vec)),
      is.finite(as.numeric(logl)[1])
    ),
    value = c(
      "TRUE",
      as.character(fit[["flp"]] %||% NA),
      as.character(nstep_flp),
      as.character(nstep_kde),
      paste(as.numeric(hdef), collapse = ";"),
      as.character(length(rho %||% numeric(0))),
      as.character(length(lambda_vec %||% numeric(0))),
      as.character(as.numeric(logl)[1])
    ),
    stringsAsFactors = FALSE
  )
}

list_flp_result_files <- function(out_root = FLP_OUT_ROOT,
                                  blocks = FLP_BLOCKS) {
  if (!dir.exists(out_root)) stop("Output root not found: ", out_root)

  files <- unlist(lapply(blocks, function(block) {
    fits_dir <- file.path(out_root, block, "fits")
    if (!dir.exists(fits_dir)) {
      warning("fits directory not found and will be skipped: ", fits_dir)
      return(character(0))
    }
    list.files(fits_dir, pattern = "\\.rds$", full.names = TRUE)
  }), use.names = FALSE)

  sort(unique(files))
}

scan_flp_status <- function(files,
                            out_csv = file.path(FLP_OUT_ROOT, "flp_status_scan.csv"),
                            write_csv = TRUE) {
  if (!length(files)) {
    out <- data.frame()
  } else {
    out <- flp_rbind_fill(lapply(files, flp_status_one_file))
  }
  if (isTRUE(write_csv)) {
    flp_dir_create(dirname(out_csv))
    utils::write.csv(out, out_csv, row.names = FALSE)
    flp_log("Status scan written to: %s", out_csv)
  }
  out
}

print_flp_status_summary <- function(status_df) {
  if (!nrow(status_df)) {
    flp_log("No result files found.", .level = "WARN")
    return(invisible(status_df))
  }

  flp_log("Files scanned: %d", nrow(status_df))
  flp_log("Readable files: %d", sum(flp_true(status_df$read_ok), na.rm = TRUE))
  flp_log("Simulation OK: %d", sum(flp_true(status_df$simulation_ok), na.rm = TRUE))
  flp_log("Existing FLP fits: %d", sum(flp_true(status_df$has_flp), na.rm = TRUE))
  flp_log("Successful FLP fits: %d", sum(flp_true(status_df$flp_ok), na.rm = TRUE))

  by_block <- stats::aggregate(
    cbind(read_ok, simulation_ok, has_flp, flp_ok) ~ block,
    data = transform(
      status_df,
      read_ok = as.integer(flp_true(read_ok)),
      simulation_ok = as.integer(flp_true(simulation_ok)),
      has_flp = as.integer(flp_true(has_flp)),
      flp_ok = as.integer(flp_true(flp_ok))
    ),
    FUN = sum,
    na.rm = TRUE
  )

  print(by_block)
  invisible(status_df)
}

## ================================================================
## 5. Add FLP to one existing .rds object
## ================================================================

add_flp_to_one_file <- function(file,
                                force = FLP_FORCE,
                                backup = FLP_BACKUP,
                                capture_output = FLP_CAPTURE_OUTPUT,
                                recompute_payload = TRUE,
                                dry_run = FALSE) {
  block <- flp_block_from_file(file)
  block_dir <- dirname(dirname(file))
  block_log <- file.path(block_dir, "progress_log_flp.csv")
  global_log <- file.path(FLP_OUT_ROOT, "progress_log_flp_all_blocks.csv")

  obj <- tryCatch(readRDS(file), error = function(e) {
    stop("Could not read RDS file: ", file, "\n", conditionMessage(e))
  })

  obj$fits <- obj$fits %||% list()
  meta <- flp_pluck(obj, "metadata", list())
  bg_case <- as.character(flp_pluck(meta, "bg_case", NA_character_))
  count_case <- as.character(flp_pluck(meta, "count_case", NA_character_))
  rep_id <- as.integer(flp_pluck(meta, "rep", NA_integer_))

  if (is.null(obj$simulation) || !isTRUE(obj$simulation$ok)) {
    row <- data.frame(
      time = as.character(Sys.time()),
      action = "skip_simulation_failed_or_missing",
      block = block,
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      nobs = NA_integer_,
      flp_ok = NA,
      flp_elapsed = NA_real_,
      nstep_flp = NA_integer_,
      hdef_x = NA_real_,
      hdef_y = NA_real_,
      error_message = as.character(flp_pluck(obj$simulation, "error_message", NA_character_)),
      stringsAsFactors = FALSE
    )
    flp_append_csv_row(row, block_log)
    flp_append_csv_row(row, global_log)
    flp_log("Skipping failed/missing simulation: %s", basename(file), .level = "WARN")
    return(invisible(row))
  }

  existing_flp_ok <- !is.null(obj$fits$flp) && isTRUE(obj$fits$flp$ok)

  ## Skip only successful existing FLP fits. If obj$fits$flp exists but
  ## ok = FALSE (for example from a previous failed test run), refit it.
  if (!isTRUE(force) && existing_flp_ok) {
    if (isTRUE(recompute_payload) && (is.null(obj$payload$fit_summaries) || !("flp" %in% names(obj$payload$fit_summaries)))) {
      flp_log("Successful FLP exists but payload lacks FLP; recomputing payload: %s", basename(file))
      bg_case_metric <- if (is.na(bg_case) || !nzchar(bg_case)) NULL else bg_case
      obj$payload <- make_metric_payload(obj$simulation$value, obj$fits, bg_case = bg_case_metric, calibration = obj$calibration)
      obj$grid_truth <- tryCatch(if (!is.null(bg_case_metric)) make_grid_eval(bg_case_metric, grid_n = GRID_N) else obj$grid_truth %||% NULL, error = function(e) obj$grid_truth %||% NULL)
      if (!isTRUE(dry_run)) flp_saveRDS_atomic(obj, file)
    }

    row <- data.frame(
      time = as.character(Sys.time()),
      action = "skip_existing_successful_flp",
      block = block,
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      nobs = tryCatch(nrow(get_cat_sim(obj$simulation$value)), error = function(e) NA_integer_),
      flp_ok = TRUE,
      flp_elapsed = obj$fits$flp$elapsed_sec,
      nstep_flp = obj$fits$flp$value[["nstep.flp"]],
      hdef_x = flp_as_scalar(obj$fits$flp$value$hdef[1]),
      hdef_y = flp_as_scalar(obj$fits$flp$value$hdef[2]),
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    flp_append_csv_row(row, block_log)
    flp_append_csv_row(row, global_log)
    flp_log("Skipping existing successful FLP: %s", basename(file))
    return(invisible(row))
  }

  if (!isTRUE(force) && !is.null(obj$fits$flp) && !existing_flp_ok) {
    flp_log("Existing FLP found but not successful; refitting: %s", basename(file), .level = "WARN")
  }

  sig_before <- flp_simulation_signature(obj$simulation$value)
  cat_sim_df <- make_eqcat_from_sim(obj$simulation$value)

  if ("mark_cat" %in% names(get_cat_sim(obj$simulation$value))) {
    cat_sim_df$mark_cat <- factor(get_cat_sim(obj$simulation$value)$mark_cat, levels = c("A", "B"))
  }

  if (!all(c("time", "lat", "long", "z", "magn1") %in% names(cat_sim_df))) {
    stop("Reconstructed catalogue lacks required eqcat columns: ", file)
  }

  if (!is.na(sig_before$nobs) && nrow(cat_sim_df) != sig_before$nobs) {
    stop("Internal check failed: reconstructed catalogue has a different number of events.")
  }

  starts <- make_flp_starts_from_obj(obj)

  if (isTRUE(dry_run)) {
    row <- data.frame(
      time = as.character(Sys.time()),
      action = "dry_run_would_fit_flp",
      block = block,
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      nobs = nrow(cat_sim_df),
      flp_ok = NA,
      flp_elapsed = NA_real_,
      nstep_flp = NA_integer_,
      hdef_x = NA_real_,
      hdef_y = NA_real_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    flp_log("DRY RUN: would fit FLP on %s", basename(file))
    return(invisible(row))
  }

  if (isTRUE(backup)) {
    bak <- paste0(file, ".before_flp.bak")
    if (!file.exists(bak)) {
      ok_bak <- file.copy(file, bak, overwrite = FALSE)
      if (!isTRUE(ok_bak)) warning("Could not create backup: ", bak)
    }
  }

  fit_log_file <- file.path(block_dir, "logs_flp", paste0(tools::file_path_sans_ext(basename(file)), "__flp_fit.log"))

  flp_log("Fitting ETAS-FLP: block=%s | bg=%s | count=%s | rep=%03d | n=%d",
          block, bg_case, count_case, rep_id, nrow(cat_sim_df))

  obj$fits$flp <- flp_safe_run(function() {
    if (isTRUE(capture_output)) {
      flp_dir_create(dirname(fit_log_file))
      con <- file(fit_log_file, open = "wt")
      sink(con, split = FALSE)
      on.exit({
        sink()
        close(con)
      }, add = TRUE)
    }

    fit_etas_flp_full(
      cat_sim_df = cat_sim_df,
      starts = starts,
      ndeclust = if (is.na(FLP_NDECLUST)) NDECLUST_FIT_NEW else FLP_NDECLUST,
      iterlim = if (is.na(FLP_ITERLIM)) ITERLIM_FIT_NEW else FLP_ITERLIM,
      n.iterweight = FLP_NITERWEIGHT,
      epsmax = FLP_EPSMAX,
      ntheta = FLP_NTHETA,
      m1 = NULL
    )
  })

  fit_checks <- check_flp_fit(obj$fits$flp, n_events = nrow(cat_sim_df))
  all_checks_ok <- all(fit_checks$ok, na.rm = FALSE)

  obj$checks <- obj$checks %||% list()
  obj$checks$flp <- list(
    created_at = Sys.time(),
    source_file = file,
    fit_log_file = if (isTRUE(capture_output)) fit_log_file else NA_character_,
    simulation_signature_before = sig_before,
    n_events_input = nrow(cat_sim_df),
    fit_checks = fit_checks,
    all_checks_ok = all_checks_ok
  )

  if (isTRUE(recompute_payload)) {
    bg_case_metric <- if (is.na(bg_case) || !nzchar(bg_case)) NULL else bg_case
    obj$payload <- make_metric_payload(obj$simulation$value, obj$fits, bg_case = bg_case_metric, calibration = obj$calibration)
    obj$grid_truth <- tryCatch(if (!is.null(bg_case_metric)) make_grid_eval(bg_case_metric, grid_n = GRID_N) else obj$grid_truth %||% NULL, error = function(e) obj$grid_truth %||% NULL)
  }

  sig_after <- flp_simulation_signature(obj$simulation$value)
  same_signature <- isTRUE(all.equal(sig_before, sig_after, check.attributes = FALSE))
  obj$checks$flp$simulation_signature_after <- sig_after
  obj$checks$flp$simulation_signature_unchanged <- same_signature

  flp_saveRDS_atomic(obj, file)

  ## Re-read the file as a final check that the saved object is usable.
  reloaded <- tryCatch(readRDS(file), error = function(e) NULL)
  saved_ok <- !is.null(reloaded) && !is.null(reloaded$fits$flp)

  fit_value <- if (isTRUE(obj$fits$flp$ok)) obj$fits$flp$value else NULL
  hdef <- if (!is.null(fit_value)) fit_value$hdef else c(NA_real_, NA_real_)

  row <- data.frame(
    time = as.character(Sys.time()),
    action = "fit_flp",
    block = block,
    file = file,
    bg_case = bg_case,
    count_case = count_case,
    rep = rep_id,
    nobs = nrow(cat_sim_df),
    flp_ok = isTRUE(obj$fits$flp$ok),
    flp_elapsed = as.numeric(obj$fits$flp$elapsed_sec),
    nstep_flp = if (!is.null(fit_value)) as.integer(fit_value[["nstep.flp"]] %||% NA_integer_) else NA_integer_,
    hdef_x = flp_as_scalar(hdef[1]),
    hdef_y = flp_as_scalar(hdef[2]),
    all_checks_ok = isTRUE(all_checks_ok),
    simulation_signature_unchanged = isTRUE(same_signature),
    saved_ok = isTRUE(saved_ok),
    fit_log_file = if (isTRUE(capture_output)) fit_log_file else NA_character_,
    error_message = as.character(obj$fits$flp$error_message %||% NA_character_),
    warnings = paste(unique(obj$fits$flp$warnings %||% character(0)), collapse = " | "),
    stringsAsFactors = FALSE
  )

  flp_append_csv_row(row, block_log)
  flp_append_csv_row(row, global_log)

  if (isTRUE(obj$fits$flp$ok)) {
    flp_log("ETAS-FLP done: ok=TRUE | elapsed=%.2f sec | checks_ok=%s | %s",
            obj$fits$flp$elapsed_sec, all_checks_ok, basename(file))
  } else {
    flp_log("ETAS-FLP failed: %s | %s", basename(file), obj$fits$flp$error_message, .level = "ERROR")
  }

  invisible(row)
}

## ================================================================
## 6. Process many files
## ================================================================

add_flp_to_existing_results <- function(out_root = FLP_OUT_ROOT,
                                        blocks = FLP_BLOCKS,
                                        force = FLP_FORCE,
                                        max_files = FLP_MAX_FILES,
                                        backup = FLP_BACKUP,
                                        capture_output = FLP_CAPTURE_OUTPUT,
                                        dry_run = FALSE) {
  files <- list_flp_result_files(out_root = out_root, blocks = blocks)

  if (!length(files)) {
    flp_log("No .rds result files found under selected blocks.", .level = "WARN")
    return(invisible(data.frame()))
  }

  status <- scan_flp_status(files, write_csv = TRUE)
  print_flp_status_summary(status)

  if (!isTRUE(force)) {
    ## Refit files with no successful FLP fit. This includes files with
    ## obj$fits$flp present but ok = FALSE from a previous failed test.
    todo <- status$file[flp_true(status$read_ok) & flp_true(status$simulation_ok) & !flp_true(status$flp_ok)]
  } else {
    todo <- status$file[flp_true(status$read_ok) & flp_true(status$simulation_ok)]
  }

  if (!length(todo)) {
    flp_log("No files need FLP fitting under the current force=%s setting.", force)
    return(invisible(status))
  }

  if (!is.na(max_files) && max_files > 0L) {
    todo <- head(todo, max_files)
  }

  flp_log("Files selected for FLP fitting: %d", length(todo))

  rows <- vector("list", length(todo))
  for (ii in seq_along(todo)) {
    flp_log("FLP file %d/%d", ii, length(todo))
    rows[[ii]] <- tryCatch(
      add_flp_to_one_file(
        file = todo[ii],
        force = force,
        backup = backup,
        capture_output = capture_output,
        dry_run = dry_run
      ),
      error = function(e) {
        row <- data.frame(
          time = as.character(Sys.time()),
          action = "fatal_error",
          block = flp_block_from_file(todo[ii]),
          file = todo[ii],
          bg_case = NA_character_,
          count_case = NA_character_,
          rep = NA_integer_,
          nobs = NA_integer_,
          flp_ok = FALSE,
          flp_elapsed = NA_real_,
          nstep_flp = NA_integer_,
          hdef_x = NA_real_,
          hdef_y = NA_real_,
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
        flp_append_csv_row(row, file.path(out_root, "progress_log_flp_all_blocks.csv"))
        flp_log("Fatal error on file: %s | %s", basename(todo[ii]), conditionMessage(e), .level = "ERROR")
        row
      }
    )
    gc(verbose = FALSE)
  }

  out <- flp_rbind_fill(rows)
  summary_csv <- file.path(out_root, "flp_run_last_summary.csv")
  utils::write.csv(out, summary_csv, row.names = FALSE)
  flp_log("Last run summary written to: %s", summary_csv)

  final_status <- scan_flp_status(files, out_csv = file.path(out_root, "flp_status_scan_after_last_run.csv"), write_csv = TRUE)
  print_flp_status_summary(final_status)

  invisible(out)
}

## ================================================================
## 7. Preflight checks and autorun modes
## ================================================================

run_flp_preflight_checks <- function() {
  flp_log("Base directory: %s", FLP_BASE_DIR)
  flp_log("Plan script: %s", FLP_PLAN_SCRIPT)
  flp_log("Output root: %s", FLP_OUT_ROOT)
  flp_log("Blocks: %s", paste(FLP_BLOCKS, collapse = ", "))
  flp_log("Mode: %s | force=%s | max_files=%s | backup=%s | capture_output=%s",
          FLP_MODE, FLP_FORCE, as.character(FLP_MAX_FILES), FLP_BACKUP, FLP_CAPTURE_OUTPUT)
  flp_log("FLP controls: ndeclust=%s | iterlim=%s | n.iterweight=%s | ntheta=%s | epsmax=%s",
          as.character(if (is.na(FLP_NDECLUST)) NDECLUST_FIT_NEW else FLP_NDECLUST),
          as.character(if (is.na(FLP_ITERLIM)) ITERLIM_FIT_NEW else FLP_ITERLIM),
          as.character(FLP_NITERWEIGHT),
          as.character(FLP_NTHETA),
          as.character(FLP_EPSMAX))

  if (!dir.exists(FLP_OUT_ROOT)) stop("Output root does not exist: ", FLP_OUT_ROOT)

  block_dirs <- file.path(FLP_OUT_ROOT, FLP_BLOCKS)
  missing_blocks <- block_dirs[!dir.exists(block_dirs)]
  if (length(missing_blocks)) {
    warning("Some block directories do not exist and will be skipped:\n  ", paste(missing_blocks, collapse = "\n  "))
  }

  files <- list_flp_result_files(FLP_OUT_ROOT, FLP_BLOCKS)
  status <- scan_flp_status(files, write_csv = TRUE)
  print_flp_status_summary(status)
  invisible(status)
}

import_etasflp_internals <- function() {
  ns <- asNamespace("etasFLP")
  
  needed <- c(
    "flp1.etas.nlmNEW",
    "flpkspace"
  )
  
  for (nm in needed) {
    assign(nm, get(nm, envir = ns), envir = .GlobalEnv)
  }
  
  invisible(TRUE)
}

import_etasflp_internals()

## Load definitions before any mode action.
load_etas_plan_definitions_no_autorun()

FLP_MAX_FILES <- NA
FLP_MODE <- "run"

if (identical(FLP_MODE, "check")) {
  flp_log("ETAS_FLP_MODE='check': scanning existing files only; no FLP fit will be run.")
  run_flp_preflight_checks()

} else if (identical(FLP_MODE, "test")) {
  flp_log("ETAS_FLP_MODE='test': fitting FLP on one selected catalogue only.")
  if (is.na(FLP_MAX_FILES) || FLP_MAX_FILES < 1L) FLP_MAX_FILES <- 1L
  run_flp_preflight_checks()
  add_flp_to_existing_results(max_files = FLP_MAX_FILES, dry_run = FALSE)

} else if (identical(FLP_MODE, "run")) {
  flp_log("ETAS_FLP_MODE='run': fitting FLP on all selected catalogues without a successful obj$fits$flp.")
  run_flp_preflight_checks()
  add_flp_to_existing_results(dry_run = FALSE)

} else if (identical(FLP_MODE, "dry_run")) {
  flp_log("ETAS_FLP_MODE='dry_run': selecting files but not fitting/saving.")
  run_flp_preflight_checks()
  add_flp_to_existing_results(dry_run = TRUE)

} else if (identical(FLP_MODE, "summary")) {
  flp_log("ETAS_FLP_MODE='summary': scanning status only.")
  files <- list_flp_result_files(FLP_OUT_ROOT, FLP_BLOCKS)
  status <- scan_flp_status(files, write_csv = TRUE)
  print_flp_status_summary(status)

} else if (identical(FLP_MODE, "none")) {
  flp_log("ETAS_FLP_MODE='none': functions loaded; no scan and no fit launched.")

} else {
  stop("Unknown ETAS_FLP_MODE: ", FLP_MODE,
       ". Allowed values: check, test, run, dry_run, summary, none.")
}
