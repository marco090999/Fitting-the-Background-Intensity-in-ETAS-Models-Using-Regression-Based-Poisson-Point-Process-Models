

#######################################################################
# STANDALONE FINAL SIMULATION PLAN: PARAMETRIC-BACKGROUND ETAS vs CLASSICAL ETAS
#
# This file is self-contained: it does NOT source/import the previous
# poster_simulation_plan_etas_parametric_full_functions_v4.R file.
# All core functions needed for simulation, classical ETAS fitting,
# parametric-background ETAS fitting, calibration, checkpointing, and
# result collection are included below.
#
# Recommended execution:
#   Sys.setenv(ETAS_RUN_MODE = "main");    source("etas_parametric_final_plan_no_smooth.R")
#   Sys.setenv(ETAS_RUN_MODE = "sample");  source("etas_parametric_final_plan_no_smooth.R")
#   Sys.setenv(ETAS_RUN_MODE = "misspec"); source("etas_parametric_final_plan_no_smooth.R")
#   Sys.setenv(ETAS_RUN_MODE = "extreme"); source("etas_parametric_final_plan_no_smooth.R")
#   Sys.setenv(ETAS_RUN_MODE = "none");    source("etas_parametric_final_plan_no_smooth.R")
#######################################################################


#######################################################################
# LEGACY CORE FUNCTIONS COPIED FROM THE PREVIOUS FILE
#######################################################################

#######################################################################
# POSTER SIMULATION PLAN: PARAMETRIC-BACKGROUND ETAS vs CLASSICAL ETAS
#
# Versione autosufficiente: tutte le funzioni sono definite direttamente nel file; nessuna patch a runtime.
# Le funzioni etasclass, etasclass.par, WPI_bkgd.fit, build_points_covs,
# etas.starting, etas.mod2NEW e etas.par.sim_v4 sono riportate direttamente
# in forma modificata.
#
# Piano finale no-smooth:
#   - background truth: constant, linear_xy, cov1, mark_cat, hotspot.
#   - main plan: N1000 bgdom/balanced/trigdom, 50 repliche.
#   - secondary plans: sample-size, misspecification, extreme-triggering, 30 repliche.
#   - checkpoint dopo ogni simulazione, fit classico, fit parametrico e metriche.
#   - il caso smooth è escluso dai piani finali.
#######################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(spatstat)
  library(stopp)
  library(stpp)
  library(mgcv)
  library(sf)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(plotly)
  library(Rcpp)
  library(etasFLP)
  library(tidyr)
  library(grid)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

cartesian.product.catmark <- function(pp_obj, win_obj, list_levels_marks) {
  
  if (!length(list_levels_marks)) {
    stop("'list_levels_marks' is empty: no categorical mark levels supplied.")
  }
  
  levels_df <- do.call(
    expand.grid,
    c(
      list_levels_marks,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  
  ## preserva i nomi delle marks
  if (!is.null(names(list_levels_marks))) {
    names(levels_df) <- names(list_levels_marks)
  }
  
  ## forza i livelli come factor coerenti
  for (nm in names(levels_df)) {
    levels_df[[nm]] <- factor(
      levels_df[[nm]],
      levels = list_levels_marks[[nm]]
    )
  }
  
  n_comb <- nrow(levels_df)
  
  ############################################################
  ## Caso 2D: ppp
  ############################################################
  
  if (inherits(pp_obj, "ppp")) {
    
    n_pts <- pp_obj$n
    
    x <- rep(pp_obj$x, each = n_comb)
    y <- rep(pp_obj$y, each = n_comb)
    
    marks_rep <- levels_df[
      rep.int(seq_len(n_comb), times = n_pts),
      ,
      drop = FALSE
    ]
    
    out <- spatstat.geom::ppp(
      x = x,
      y = y,
      window = win_obj,
      marks = marks_rep,
      check = FALSE
    )
    
    return(out)
  }
  
  ############################################################
  ## Caso 3D: pp3
  ############################################################
  
  if (inherits(pp_obj, "pp3")) {
    
    df <- as.data.frame(pp_obj$data)
    
    n_pts <- nrow(df)
    
    x <- rep(df$x, each = n_comb)
    y <- rep(df$y, each = n_comb)
    z <- rep(df$z, each = n_comb)
    
    marks_rep <- levels_df[
      rep.int(seq_len(n_comb), times = n_pts),
      ,
      drop = FALSE
    ]
    
    out <- spatstat.geom::pp3(
      x = x,
      y = y,
      z = z,
      box = win_obj,
      marks = marks_rep
    )
    
    return(out)
  }
  
  stop("'pp_obj' must be either a 'ppp' or a 'pp3' object.")
}



catmark_as_df <- function(obj) {
  
  ############################################################
  ## Caso 2D: ppp
  ############################################################
  
  if (inherits(obj, "ppp")) {
    
    df <- data.frame(
      x = obj$x,
      y = obj$y,
      check.names = FALSE
    )
    
    m <- spatstat.geom::marks(obj)
    
    if (!is.null(m)) {
      
      if (is.data.frame(m)) {
        mdf <- m
      } else if (is.matrix(m)) {
        mdf <- as.data.frame(m, stringsAsFactors = FALSE)
      } else {
        mdf <- data.frame(mark = m, check.names = FALSE)
      }
      
      df <- cbind(df, mdf)
    }
    
    return(df)
  }
  
  ############################################################
  ## Caso 3D: pp3
  ############################################################
  
  if (inherits(obj, "pp3")) {
    
    df <- as.data.frame(obj$data, stringsAsFactors = FALSE)
    
    m <- spatstat.geom::marks(obj)
    
    if (!is.null(m)) {
      
      if (is.data.frame(m)) {
        mdf <- m
      } else if (is.matrix(m)) {
        mdf <- as.data.frame(m, stringsAsFactors = FALSE)
      } else {
        mdf <- data.frame(mark = m, check.names = FALSE)
      }
      
      ## Se le marks non sono già dentro obj$data, le aggiungo.
      new_cols <- setdiff(names(mdf), names(df))
      
      if (length(new_cols)) {
        df <- cbind(df, mdf[, new_cols, drop = FALSE])
      }
    }
    
    return(df)
  }
  
  if (!is.null(obj$data)) {
    return(as.data.frame(obj$data, stringsAsFactors = FALSE))
  }
  
  as.data.frame(obj, stringsAsFactors = FALSE)
}



extract_ppp_marks_for_WPI <- function(X0, formula = NULL) {
  
  m <- spatstat.geom::marks(X0)
  
  if (is.null(m)) {
    return(NULL)
  }
  
  if (is.data.frame(m)) {
    mdf <- m
  } else if (is.matrix(m)) {
    mdf <- as.data.frame(m, stringsAsFactors = FALSE)
  } else {
    mdf <- data.frame(.mark = m)
  }
  
  ## Se c'è una sola mark senza nome informativo, uso il nome dalla formula
  if (!is.null(formula) && ncol(mdf) == 1) {
    vars_formula <- all.vars(formula)
    mark_vars <- setdiff(vars_formula, c("x", "y", "t", "z"))
    
    if (length(mark_vars) == 1) {
      names(mdf) <- mark_vars
    }
  }
  
  mdf
}



rbind_fill_base <- function(...) {
  dots <- list(...)
  dots <- dots[!vapply(dots, is.null, logical(1))]
  if (!length(dots)) return(data.frame())
  dots <- lapply(dots, as.data.frame, stringsAsFactors = FALSE)
  dots <- dots[vapply(dots, nrow, integer(1)) > 0L]
  if (!length(dots)) return(data.frame())
  all_names <- unique(unlist(lapply(dots, names), use.names = FALSE))
  dots <- lapply(dots, function(d) {
    miss <- setdiff(all_names, names(d))
    if (length(miss)) {
      for (m in miss) d[[m]] <- NA
    }
    d[, all_names, drop = FALSE]
  })
  rownames_out <- NULL
  out <- do.call(base::rbind, dots)
  rownames(out) <- rownames_out
  out
}

.safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

saveRDS_atomic <- function(object, file, compress = "gzip") {
  .safe_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) {
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(file)
}

append_csv_row <- function(row, file) {
  .safe_dir_create(dirname(file))
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

#######################################################################
# 1) FUNZIONI ETAS CLASSICHE RISCRITTE/PATCHATE
#######################################################################

etasclass <- function (cat.orig, time.update = FALSE, magn.threshold = 2.5, 
                       magn.threshold.back = magn.threshold + 2, tmax = max(cat.orig$time), 
                       long.range = range(cat.orig$long), lat.range = range(cat.orig$lat), 
                       mu = 1, k0 = 1, c = 0.5, p = 1.01, gamma = 0.5, d = 1, q = 1.5, 
                       betacov = 0.7, params.ind = replicate(7, TRUE), formula1 = "time~magnitude-1", 
                       offset = 0, hdef = c(1, 1), w = replicate(nrow(cat.orig), 1), hvarx = replicate(nrow(cat.orig), 1), 
                       hvary = replicate(nrow(cat.orig),1), declustering = TRUE, thinning = FALSE, flp = TRUE, 
                       m1 = NULL, ndeclust = 5, n.iterweight = 1, onlytime = FALSE, 
                       is.backconstant = FALSE, description = "", cat.back = NULL, 
                       back.smooth = 1, sectoday = FALSE, longlat.to.km = TRUE, 
                       usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, 
                       iterlim = 50, ntheta = 36) 
{
  iprint <- FALSE
  fastML <- FALSE
  parallel <- FALSE
  fast.eps <- 0.001
  params.lim = c(0, 0, 0, 1, 0, 0, 1)
  this.call <- match.call()
  flag <- eqcat(cat.orig)
  if (!flag$ok) 
    stop("WRONG EARTHQUAKE CATALOG DEFINITION")
  cat.orig <- flag$cat
  iter <- 0
  AIC.iter <- numeric(0)
  AIC.iter2 <- numeric(0)
  params.iter <- numeric(0)
  sqm.iter <- numeric(0)
  rho.weights.iter <- numeric(0)
  hdef.iter <- numeric(0)
  wmat <- numeric(0)
  fl <- 0
  fl.iter <- numeric(0)
  AIC.decrease <- TRUE
  trace <- TRUE
  eps <- 2 * epsmax
  eps.par <- 2 * epsmax
  if (onlytime) {
    is.backconstant <- TRUE
    declustering <- FALSE
    params.ind[5:7] <- c(0, 0, 0)
    gamma <- 0
    d <- 0
    q <- 0
  }
  if (is.backconstant) {
    declustering <- FALSE
  }
  if (!declustering) {
    ndeclust <- 1
    thinning <- FALSE
    flp <- FALSE
  }
  if (missing(tmax)) 
    is.na(tmax) <- TRUE
  if (missing(long.range)) 
    is.na(long.range) <- TRUE
  if (missing(lat.range)) 
    is.na(lat.range) <- TRUE
  etas.current <- list(parallel = parallel, this.call = match.call(), 
                       nstep.flp = 0, nstep.kde = 0, nstep.par = 0, description = description, 
                       time.start = Sys.time(), magn.threshold = magn.threshold, 
                       magn.threshold.back = magn.threshold.back, onlytime = onlytime, 
                       tmax = tmax, lat.range = lat.range, long.range = long.range, 
                       hvarx = hvarx, hvary = hvary, is.backconstant = is.backconstant, 
                       usenlm = usenlm, method = method, cat.orig = cat.orig, 
                       declustering = declustering, thinning = thinning, flp = flp, 
                       back.smooth = back.smooth, ndeclust = ndeclust, eps = eps, 
                       longlat.to.km = longlat.to.km, sectoday = sectoday, ntheta = ntheta)
  class(etas.current) <- "etasclass"
  etas.current <- cat.select(etas.current, longlat.to.km, sectoday = sectoday)
  etas.current$cat <- data.frame(etas.current$cat, hvarx, hvary)
  cat <- etas.current$cat[etas.current$cat$ord, ][etas.current$cat$ind[etas.current$cat$ord], 
  ]
  etas.current$cat <- cat
  n <- nrow(cat)
  missingw <- missing(w)
  if (missingw) 
    w = replicate(n, 1)
  if (time.update) {
    n.old <- length(w)
    w <- c(w, replicate(n - n.old, 0.5))
    hvarx <- c(hvarx, replicate(n - n.old, 1))
    hvary <- c(hvary, replicate(n - n.old, 1))
  }
  else {
    hvarx <- cat$hvarx/(prod(cat$hvarx))^(1/n)
    hvary <- cat$hvary/(prod(cat$hvary))^(1/n)
    etas.current$w <- w
  }
  ycat.work <- cat$ycat.work
  xcat.work <- cat$xcat.work
  if (missing(cat.back) || is.null(cat.back)) 
    ind.back <- (cat$magn1 >= magn.threshold.back)
  else stop("cat.back argument no more allowed, computed internally")
  xback.work <- cat$xcat.work[ind.back]
  yback.work <- cat$ycat.work[ind.back]
  if (!missingw) {
    xback.work <- xcat.work
    yback.work <- ycat.work
  }
  if (missingw) 
    w = replicate(length(xback.work), 1)
  starting <- etas.starting(cat, magn.threshold = magn.threshold, 
                            longlat.to.km = longlat.to.km, sectoday = sectoday, p.start = p, 
                            gamma.start = gamma, q.start = q, betacov.start = betacov[1], 
                            onlytime = onlytime)
  if (missing(m1) || is.null(m1)) 
    m1 <- as.integer(nrow(cat)/2)
  if (missing(mu) || is.null(mu)) 
    mu <- starting$mu.start
  if (missing(k0) || is.null(k0)) 
    k0 <- starting$k0.start
  if (missing(c) || is.null(c)) 
    c <- starting$c.start
  if (missing(d) || is.null(d)) 
    d <- starting$d.start
  if (missing(hdef) || is.null(hdef)) 
    hdef <- c(bwd.nrd(xback.work, w), bwd.nrd(yback.work, 
                                              w))
  print("Initial ETAS params estimates after checking: ")
  vecpar.init = c(mu, k0, c, p, gamma, d, q)
  names(vecpar.init) = c("mu", "k0", "c", "p", "gamma", "d", 
                         "q")
  print(round(vecpar.init, 4))
  etas.current$xback.work = xback.work
  etas.current$yback.work = yback.work
  etas.current$mu.start = mu
  etas.current$k0.start = k0
  etas.current$c.start = c
  etas.current$d.start = d
  etas.current$p.start = p
  etas.current$q.start = q
  etas.current$gamma.start = gamma
  etas.current$betacov.start = betacov
  etas.current$hdef.start = hdef
  etas.current$formula1 = formula1
  formula1 = as.formula(formula1)
  formula1 = update(formula1, . ~ . - 1)
  cov.matrix = model.matrix(formula1, data = cat)
  offset = model.offset(model.frame(formula1, data = cat))
  if (is.null(offset)) 
    offset = replicate(nrow(cat), 0)
  if (!is.element("matrix", class(cov.matrix))) 
    stop("WRONG FORMULA DEFINITION")
  ncov = ncol(cov.matrix)
  if (length(betacov) != ncov) {
    cat("Wrong number of starting values for parametrs linear predictor of covariates. Correct number of zero starting values inserted", 
        "\n")
    betacov = replicate(ncov, 0)
  }
  if (length(params.ind) != 7) 
    stop("Wrong number of elements of params.ind")
  params.ind <- as.numeric(params.ind)
  if (sum(abs(params.ind - 0.5) == 0.5) != 7) 
    stop("WRONG params.ind DEFINITION: ONLY FALSE/TRUE ALLOWED SEE HELP")
  namespar = c("mu", "k0", "c", "p", "gamma", "d", "q")
  params.fix <- c(mu, k0, c, p, gamma, d, q)
  nparams.etas <- sum(params.ind)
  nparams <- nparams.etas + ncov
  if (onlytime) {
    if (!(prod(params.fix[1:4] >= params.lim[1:4]))) 
      stop("WRONG starting values for parameters: mu, k0, c must  be positive, p must be equal or greater than 1")
    rho.s2 = 0
    region = NULL
  }
  else {
    if (!(prod(params.fix >= params.lim))) 
      stop("WRONG starting values for parameters: mu, k0, c, gamma, d, q must  be positive, p must be equal or greater than 1")
    region = embedding.rect.cat.epsNEW(cat)
    rho.s2 = matrix(0, ntheta, n)
    for (i in (1:n)) rho.s2[, i] = region.P(region, c(xcat.work[i], 
                                                      ycat.work[i]), k = ntheta)$rho
    rho.s2 = t(rho.s2)
  }
  etas.current$cov.matrix = cov.matrix
  etas.current$region = region
  etas.current$rho.s2 = rho.s2
  etas.current$offset = offset
  if (is.backconstant) {
    if (!onlytime) 
      back.dens <- array(1, n)/(diff(range(xcat.work)) * 
                                  diff(range(ycat.work)))
    else back.dens = 1
    back.integral <- 1
  }
  else {
    back.tot <- kde2dnew.fortran(xback.work, yback.work, 
                                 xcat.work, ycat.work, h = hdef, factor.xy = back.smooth, 
                                 eps = 1/n, w = w, hvarx = hvarx, hvary = hvary)
    etas.current$nstep.kde = etas.current$nstep.kde + 1
    back.dens <- back.tot$z
    back.integral <- back.tot$integral
  }
  etas.current$back.integral <- back.integral
  etas.current$back.dens <- back.dens
  while ((iter < ndeclust) & ((eps > epsmax) || (eps.par > 
                                                 epsmax))) {
    params <- log((params.fix - params.lim)[params.ind == 1])
    etas.current$params <- params
    etas.current$params.lim <- params.lim
    etas.current$params.ind <- as.logical(params.ind)
    etas.current$params.fix <- params.fix
    etas.current$betacov <- betacov
    etas.current$nparams <- nparams
    etas.current$nparams.etas <- nparams.etas
    cat("Start ML step number: ", iter + 1, "\n")
    etas.current$fastML <- fastML
    etas.current$fast.eps <- fast.eps
    if (fastML) {
      fast <- (fastML.init(etas.current))
      etas.current$ind <- fast$ind
      etas.current$index.tot <- fast$index.tot
    }
    risult.opt <- generaloptimizationNEW(etas.current, hessian = TRUE, 
                                         iterlim = iterlim, iprint = iprint, trace = trace)
    etas.current <- risult.opt$etas.obj
    params.optim <- risult.opt$params.optim[1:nparams.etas]
    betacov <- risult.opt$params.optim[(nparams.etas + 1):nparams]
    cat("\n")
    for (i in 1:n.iterweight) {
      cat("weighting step n. ", i, "\n")
      l.optim <- risult.opt$l.optim
      risult <- risult.opt$risult
      l <- etas.mod2NEW(params = c(params.optim, betacov), 
                        etas.obj = etas.current, trace = FALSE)
      det.check <- det(risult$hessian)
      check = is.na(det.check) || (abs(det.check) < 1e-20)
      if (check) 
        compsqm <- FALSE
      sqm <- 0
      if (compsqm) 
        sqm <- sqrt(diag(solve(risult$hessian))) * c(exp(params.optim), 
                                                     replicate(ncov, 1))
      sqm.etas <- sqm[1:nparams.etas]
      sqm.cov <- sqm[(nparams.etas + 1):nparams]
      params <- params.fix
      params[params.ind == 1] <- exp(params.optim) + params.lim[params.ind == 
                                                                  1]
      sqm.tot <- array(0, nparams.etas)
      sqm.tot[params.ind == 1] <- sqm.etas
      names(sqm.tot) <- namespar
      params.MLtot <- c(params, betacov)
      sqm.MLtot <- c(sqm.tot, sqm.cov)
      cat("found optimum; end  ML step  ")
      cat(iter + 1, "\n")
      mu <- params[1]
      k0 <- params[2]
      c <- params[3]
      p <- params[4]
      gamma <- params[5]
      d <- params[6]
      q <- params[7]
      predictor <- as.matrix(etas.current$cov.matrix) %*% 
        as.vector(betacov) + as.vector(etas.current$offset)
      rho.weights <- params.MLtot[1] * back.dens/attr(l, 
                                                      "lambda.vec")
      w <- rho.weights
      xback.work <- xcat.work
      yback.work <- ycat.work
      if (flp) {
        etas.l <- attr(l, "etas.vec")
        ris.flp <- flp1.etas.nlmNEW(cat, h.init = hdef, 
                                    etas.params = params.MLtot, etas.l = etas.l, 
                                    w = rho.weights, m1 = m1, m2 = as.integer(nrow(cat) - 
                                                                                1), mh = 1)
        etas.current$nstep.flp = etas.current$nstep.flp + 
          1
        hdef <- ris.flp$hdef
        fl <- ris.flp$fl
      }
      if (thinning) {
        back.ind <- runif(n) < rho.weights
        xback.work <- xcat.work[back.ind]
        yback.work <- ycat.work[back.ind]
        w <- replicate(length(xback.work), 1)
      }
      if (flp) 
        back.tot <- kde2dnew.fortran(xback.work, yback.work, 
                                     xcat.work, ycat.work, eps = 1/n, h = hdef, 
                                     w = w, hvarx = hvarx, hvary = hvary)
      else back.tot <- kde2dnew.fortran(xback.work, yback.work, 
                                        xcat.work, ycat.work, eps = 1/n, w = w, hvarx = hvarx, 
                                        hvary = hvary)
      etas.current$nstep.kde = etas.current$nstep.kde + 
        1
      back.dens <- back.tot$z
      back.integral <- back.tot$integral
      wmat <- back.tot$wmat
      etas.current$rho.weights <- rho.weights
      etas.current$back.integral <- back.integral
      etas.current$back.dens <- back.dens
      l <- etas.mod2NEW(params = c(params.optim, betacov), 
                        etas.obj = etas.current, trace = FALSE)
      AIC2 <- 2 * l + 2 * nparams
    }
    time.res <- array(0, n)
    if (onlytime) {
      times.tot <- cat$time
      magnitudes.tot <- cat$magn1 - magn.threshold
      tmin <- min(times.tot)
      etas.t <- 0
      for (i in 1:n) {
        tmax.i <- times.tot[i]
        times <- subset(times.tot, times.tot < tmax.i)
        magnitudes <- subset(magnitudes.tot, times.tot < 
                               tmax.i)
        if (p == 1) {
          it <- log(c + tmax.i - times) - log(c)
        }
        else {
          it <- ((c + tmax.i - times)^(1 - p) - c^(1 - 
                                                     p))/(1 - p)
        }
        time.res[i] <- (tmax.i - tmin) * mu + k0 * sum(exp(predictor) * 
                                                         it)
      }
    }
    names(params.fix) = namespar
    names(params.ind) = namespar
    names(betacov) = colnames(cov.matrix)
    names(params.MLtot) = c(namespar, names(betacov))
    names(sqm.MLtot) = c(namespar, names(betacov))
    timenow <- Sys.time()
    time.elapsed <- difftime(timenow, etas.current$time.start, 
                             units = "secs")
    AIC.temp <- 2 * l + 2 * nparams
    AIC.decrease <- ((iter == 0) || (AIC.temp <= min(AIC.iter)))
    if (AIC.decrease) {
      iter <- iter + 1
      AIC.iter[iter] <- AIC.temp
      AIC.iter2[iter] <- AIC2
      params.iter <- rbind(params.iter, params.MLtot)
      sqm.iter <- rbind(sqm.iter, sqm.tot)
      rho.weights.iter <- rbind(rho.weights.iter, rho.weights)
      hdef.iter <- rbind(hdef.iter, hdef)
      fl.iter <- c(fl.iter, fl)
      cat("\n", "---ITERATION n.", iter, " ---; AIC = ", 
          round(AIC.iter[iter], 3), "; elapsed time: ", 
          time.elapsed, " time: ")
      print(Sys.time())
      cat("\n")
      cat("Current estimates of parameters: ", "\n")
      cat(round(params.MLtot, 5))
      cat("\n")
    }
    if (iter > 1) {
      eps <- max(abs(back.dens - etas.ris$back.dens))
      eps.par <- max(abs(params.iter[iter] - params.iter[iter - 
                                                           1]))
    }
    rownames(params.iter) <- 1:nrow(params.iter)
    predictor <- as.matrix(etas.current$cov.matrix) %*% as.vector(betacov) + 
      as.vector(etas.current$offset)
    names(params.MLtot) <- c(namespar, names(betacov))
    names(sqm.MLtot) <- names(params.MLtot)
    etas.current$time.elapsed <- time.elapsed
    etas.current$time.end <- timenow
    etas.current$risult <- risult
    etas.current$betacov <- betacov
    etas.current$params.MLtot <- params.MLtot
    etas.current$params <- params.MLtot[1:7]
    etas.current$predictor <- predictor
    etas.current$sqm <- sqm.MLtot
    etas.current$time.res <- time.res
    etas.current$AIC.iter <- AIC.iter
    etas.current$AIC.iter2 <- AIC.iter2
    etas.current$AIC.decrease <- AIC.decrease
    etas.current$params.iter <- params.iter
    etas.current$sqm.iter <- sqm.iter
    etas.current$rho.weights.iter <- rho.weights.iter
    etas.current$rho.weights <- rho.weights
    etas.current$n.iterweight <- n.iterweight
    etas.current$hdef <- hdef
    etas.current$hdef.iter <- hdef.iter
    etas.current$back.integral <- back.integral
    etas.current$back.dens <- back.dens
    etas.current$wmat <- wmat
    etas.current$fl.iter <- fl.iter
    etas.current$iter <- iter
    etas.current$usenlm <- usenlm
    etas.current$method <- method
    etas.current$epsmax <- epsmax
    etas.current$iterlim <- iterlim
    etas.current$compsqm <- compsqm
    etas.current$ntheta <- ntheta
    etas.current$logl <- l
    etas.current$l <- attr(l, "lambda.vec")
    etas.current$integral <- attr(l, "integraltot")
    etas.ris <- etas.current
    if (!AIC.decrease) {
      print("INCREASING AIC")
      return(etas.ris)
    }
  }
  if (!time.update) {
    cat("\n")
    cat("\n")
    cat("--------- FINAL SUMMARY ---------", "\n")
    summary(etas.ris)
  }
  return(etas.ris)
}



etas.starting <- function (cat.orig, magn.threshold = 2.5, p.start = 1, gamma.start = 0.5, 
          q.start = 2, betacov.start = 0.7, longlat.to.km = TRUE, sectoday = FALSE, 
          onlytime = FALSE) 
{
  cat = cat.orig[cat.orig$magn1 > magn.threshold, ]
  if (sectoday) 
    cat$time = cat$time/86400
  t = cat$time
  a = cat[order(t), ]
  t = a$time
  if (longlat.to.km) {
    radius = 6371.3
    ycat.km = radius * a$lat * pi/180
    xcat.km = radius * a$long * pi/180
  }
  else {
    ycat.km = a$lat
    xcat.km = a$long
  }
  n = nrow(a)
  dt = diff(t)
  ds = diff(ycat.km)^2 + diff(xcat.km)^2
  c.start = as.numeric(quantile(dt, 0.25))
  d.start = as.numeric(quantile(ds, 0.05))
  mu.start = n * 0.5/diff(range(t))
  if (onlytime) {
    gamma.start = 0
    q.start = 0
    d.start = 0
  }
  tmax = max(t)
  it = log((c.start + tmax - t)/c.start)
  em = exp((betacov.start[1]) * (a$magn1 - magn.threshold))
  is = d.start^(1 - q.start) * exp(gamma.start * (a$magn1 - 
                                                    magn.threshold)) * pi/(q.start - 1)
  if (onlytime) 
    k0.start = n * 0.5/sum(em * it)
  else k0.start = n * 0.5/sum(em * it * is)
  return(list(mu.start = mu.start, k0.start = k0.start, c.start = c.start, 
              p.start = p.start, gamma.start = gamma.start, d.start = d.start, 
              q.start = q.start, betacov.start = betacov.start, longlat.to.km = longlat.to.km, 
              sectoday = sectoday))
}



kde2dnew.fortran <- function (xkern, ykern, gx, gy, h, factor.xy = 1, eps = 0, w = replicate(length(xkern), 1), 
                              hvarx = replicate(length(xkern), 1), hvary = replicate(length(xkern), 1)) 
{
  n <- length(gx)
  nkern <- length(xkern)
  if (length(ykern) != nkern) 
    stop("Data vectors must be the same length")
  if (missing(h)) 
    h <- c(bwd.nrd(xkern, w), bwd.nrd(ykern, w))
  h <- h * factor.xy
  h <- ifelse(h == 0, max(h), h)
  ris2d <- .Fortran("density2serial", x = as.double(gx), y = as.double(gy), 
                    m = as.integer(n), xkern = as.double(xkern), ykern = as.double(ykern), 
                    nkern = as.integer(nkern), h = as.double(h), w = as.double(w), 
                    hvarx = as.double(hvarx), hvary = as.double(hvary), dens = as.double(gx))
  integral = kde2d.integral(xkern, ykern, gx, gy, factor.xy = factor.xy, 
                            eps = eps, w = w, h = h)
  return(list(x = gx, y = gy, z = ris2d$dens, h = h, integral = integral))
}


kde2d.integral <- function (xkern, ykern, gx = xkern, gy = ykern, eps = 0, factor.xy = 1, 
                            h = c(bwd.nrd(xkern, w), bwd.nrd(ykern, w)), w = replicate(length(xkern), 1), 
                            hvarx = replicate(length(xkern), 1), hvary = replicate(length(xkern), 1)) 
{
  eps.x = diff(range(gx)) * eps/2
  eps.y = diff(range(gy)) * eps/2
  rx = c(min(gx) - eps.x, max(gx) + eps.x)
  ry = c(min(gy) - eps.y, max(gy) + eps.y)
  hx <- factor.xy * h[1]
  hy <- factor.xy * h[2]
  nkern = length(xkern)
  ix = pnorm(rx[2], xkern, hx * hvarx) - pnorm(rx[1], xkern, 
                                               hx * hvarx)
  iy = pnorm(ry[2], ykern, hy * hvary) - pnorm(ry[1], ykern, 
                                               hy * hvary)
  integral = sum(w * ix * iy)/sum(w)
  return(integral)
}


generaloptimizationNEW <- function (etas.obj, hessian, iterlim, iprint, trace) 
{
  n = nrow(etas.obj$cat)
  params = c(etas.obj$params, etas.obj$betacov)
  etas.obj$nstep.par = etas.obj$nstep.par + 1
  if (iprint) {
    print("params in general optimization")
    print(params)
    print(etas.obj$nparams)
    print(etas.obj$nparams.etas)
  }
  etas.obj$dettrasf = replicate(n, 1)
  if (etas.obj$usenlm) {
    risult = nlm(etas.mod2NEW, params, typsize = abs(params), 
                 hessian = hessian, iterlim = iterlim, params.lim = etas.obj$params.lim, 
                 etas.obj = etas.obj)
    params.optim = risult$estimate
    l.optim = risult$minimum
  }
  else {
    risult = optim(params, etas.mod2NEW, method = etas.obj$method, 
                   params.lim = etas.obj$params.lim, hessian = hessian, 
                   control = list(trace = 2, maxit = iterlim, fnscale = n/diff(range(etas.obj$cat$time.work)), 
                                  parscale = sqrt(exp(params))), etas.obj = etas.obj)
    params.optim = risult$par
    l.optim = risult$value
  }
  return(list(params.optim = params.optim, l.optim = l.optim, 
              risult = risult, etas.obj = etas.obj))
}


etas.mod2NEW <- function (params = c(1, 1, 1, 1, 1, 1, 1), etas.obj, params.lim = c(0, 0, 0, 1, 0, 0, 1), trace = TRUE, iprint = FALSE) 
{
  params.ind = etas.obj$params.ind
  params.lim = etas.obj$params.lim
  params.fix = etas.obj$params.fix
  tmax = etas.obj$tmax
  cat = etas.obj$cat
  magn.threshold = etas.obj$magn.threshold
  back.dens = etas.obj$back.dens
  back.integral = etas.obj$back.integral
  onlytime = etas.obj$onlytime
  rho.s2 = etas.obj$rho.s2
  nparams.etas = etas.obj$nparams.etas
  nparams = etas.obj$nparams
  ntheta = etas.obj$ntheta
  params.etas = params[1:nparams.etas]
  betacov = params[(nparams.etas + 1):nparams]
  params.e = params.fix
  params.e[params.ind == 1] = exp(params.etas) + params.lim[params.ind == 
                                                              1]
  lambda = params.e[1]
  k0 = params.e[2]
  c = params.e[3]
  p = params.e[4]
  gamma = params.e[5]
  d = params.e[6]
  q = params.e[7]
  x = cat$xcat.work
  y = cat$ycat.work
  magnitudes = cat$magnitude
  predictor = as.matrix(etas.obj$cov.matrix) %*% as.vector(betacov) + 
    as.vector(etas.obj$offset)
  times = cat$time.work
  range.t = diff(range(times))
  time.init <- Sys.time()
  n = length(times)
  etas.comp = as.double(array(0, n))
  if (etas.obj$fastML) {
    ris = .Fortran("etasfull8fast", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), ind = as.integer(etas.obj$ind), 
                   nindex = as.integer(length(etas.obj$index.tot)), 
                   index = as.integer(etas.obj$index.tot), l = etas.comp)
  }
  else if (etas.obj$parallel) {
    ris = .Fortran("etasfull8newparallel", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), l = etas.comp)
  }
  else {
    ris = .Fortran("etasfull8newserial", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), l = etas.comp)
  }
  timenow <- Sys.time()
  timeelapsed <- difftime(timenow, time.init, units = "secs")
  etas.comp = ris$l
  if (iprint) 
    cat(timeelapsed, "\n")
  ci = lambda * back.dens + etas.comp
  logL.l = sum(log(ci))
  {
    if (p == 1) {
      it = log(c + tmax - times) - log(c)
    }
    else {
      it = ((c + tmax - times)^(1 - p) - c^(1 - p))/(1 - 
                                                       p)
    }
    if (onlytime) {
      integralNEW = k0 * sum(exp(predictor) * it)
    }
    else {
      m1NEW = as.vector(exp(predictor))
      m2 = as.vector(exp(gamma * magnitudes))
      time.init <- Sys.time()
      ci = lambda * back.dens + etas.comp
      etasNEW = rowSums(m1NEW * m2 * ((rho.s2 * rho.s2 * 
                                         etas.obj$dettrasf/m2 + d)^(1 - q) - d^(1 - q)))
      spaceNEW = (pi/((1 - q) * ntheta)) * etasNEW
      integralNEW = sum(k0 * it * spaceNEW)
      if (iprint) {
        cat("integral polar", "\n")
        cat(integral, "\n")
      }
    }
    integral = integralNEW
    timenow <- Sys.time()
    timeelapsed <- difftime(timenow, time.init, units = "secs")
    integraltot = integral + lambda * range.t * back.integral
  }
  logL = -logL.l + integraltot
  if (iprint) {
    cat("whole Integral computation", "\n")
    cat(timeelapsed, "\n")
    cat(params, "\n")
    cat(exp(params), "\n")
  }
  attr(logL, "etas.vec") = etas.comp
  attr(logL, "lambda.vec") = ci
  attr(logL, "integraltot") = integraltot
  if (iprint) {
    cat("ML step: likelihoods and integral", (c(logL, logL.l, 
                                                integraltot)), "\n")
  }
  if (trace) 
    cat("*")
  return(logL)
}




#### FUNZIONI UTILI ANCHE NELL'ETAS PARAMETRICO ####

etas.starting <- function (cat.orig, magn.threshold = 2.5, p.start = 1, gamma.start = 0.5, 
                           q.start = 2, coef.cov.trig.start = 0.7, longlat.to.km = TRUE, sectoday = FALSE, 
                           onlytime = FALSE) 
{
  cat = cat.orig[cat.orig$magn1 > magn.threshold, ]
  if (sectoday) 
    cat$time = cat$time/86400
  t = cat$time
  a = cat[order(t), ]
  t = a$time
  if (longlat.to.km) {
    radius = 6371.3
    ycat.km = radius * a$lat * pi/180
    xcat.km = radius * a$long * pi/180
  }
  else {
    ycat.km = a$lat
    xcat.km = a$long
  }
  n = nrow(a)
  dt = diff(t)
  ds = diff(ycat.km)^2 + diff(xcat.km)^2
  c.start = as.numeric(quantile(dt, 0.25))
  d.start = as.numeric(quantile(ds, 0.05))
  mu.start = n * 0.5/diff(range(t))
  if (onlytime) {
    gamma.start = 0
    q.start = 0
    d.start = 0
  }
  tmax = max(t)
  it = log((c.start + tmax - t)/c.start)
  em = exp((coef.cov.trig.start[1]) * (a$magn1 - magn.threshold))
  is = d.start^(1 - q.start) * exp(gamma.start * (a$magn1 - 
                                                    magn.threshold)) * pi/(q.start - 1)
  if (onlytime) 
    k0.start = n * 0.5/sum(em * it)
  else k0.start = n * 0.5/sum(em * it * is)
  return(list(mu.start = mu.start, k0.start = k0.start, c.start = c.start, 
              p.start = p.start, gamma.start = gamma.start, d.start = d.start, 
              q.start = q.start, coef.cov.trig.start = coef.cov.trig.start, longlat.to.km = longlat.to.km, 
              sectoday = sectoday))
}


embedding.rect.cat.epsNEW <- function (cat, cycle = FALSE, eps = 1/nrow(cat)) 
  embedding.rect.eps(cat$xcat.work, cat$ycat.work, cycle = cycle, eps = eps)


embedding.rect.eps <- function (x, y, cycle = FALSE, eps = 0) 
{
  eps.x = diff(range(x)) * eps/2
  eps.y = diff(range(y)) * eps/2
  rx = c(min(x) - eps.x, max(x) + eps.x)
  ry = c(min(y) - eps.y, max(y) + eps.y)
  rect = c(rx[2], rx[1], rx[1], rx[2], rx[2], ry[2], ry[2], 
           ry[1], ry[1], ry[2])
  rect = matrix(rect, 5, 2)
  if (!cycle) 
    rect = rect[1:4, ]
  return(rect)
}

eqcat <- function (x) 
{
  a <- c("time", "lat", "long", "z", "magn1")
  n <- 5 - length(intersect(a, names(x)))
  if (n > 0) {
    print("WRONG EARTHQUAKE CATALOG DEFINITION")
    print(c(n, " wrong variable names"))
    print("an object of class eqcat (earthquake catalog) must contains at least the five names")
    print(a)
    return(list(cat = x, ok = FALSE))
  }
  else {
    class(x) <- c("eqcat", "data.frame")
    return(list(cat = x, ok = TRUE))
  }
}


cat.select <- function (etas.obj, longlat.to.km, sectoday) 
{
  time.work = etas.obj$cat.orig$time
  if (sectoday) 
    time.work = time.work/86400
  unit = 1
  if (longlat.to.km) 
    unit = 6371.3 * pi/180
  ycat.work = etas.obj$cat.orig$lat * unit
  xcat.work = etas.obj$cat.orig$long * unit
  magnitude = etas.obj$cat.orig$magn1 - etas.obj$magn.threshold
  ind.magn1 = (etas.obj$cat.orig$magn1 >= etas.obj$magn.threshold)
  if (is.na(etas.obj$tmax)) 
    etas.obj$tmax = max(etas.obj$cat.orig$time[ind.magn1])
  if (is.na(etas.obj$long.range[1])) 
    etas.obj$long.range = range(etas.obj$cat.orig$long[ind.magn1])
  if (is.na(etas.obj$lat.range[1])) 
    etas.obj$lat.range = range(etas.obj$cat.orig$lat[ind.magn1])
  ind.time = (etas.obj$cat.orig$time <= etas.obj$tmax)
  ind.long = (etas.obj$cat.orig$long <= etas.obj$long.range[2]) & 
    (etas.obj$cat.orig$long >= etas.obj$long.range[1])
  ind.lat = (etas.obj$cat.orig$lat <= etas.obj$lat.range[2]) & 
    (etas.obj$cat.orig$lat >= etas.obj$lat.range[1])
  ind = as.logical(ind.magn1 * ind.time * ind.long * ind.lat)
  ord <- order(time.work)
  cat = data.frame(etas.obj$cat.orig, time.work, xcat.work, 
                   ycat.work, magnitude, ind, ord)
  etas.obj$cat = cat
  return(etas.obj)
}


bwd.nrd <- function (x, w = replicate(length(x), 1), d = 2) 
{
  if (length(x) < 2L) 
    stop("need at least 2 data points")
  m <- weighted.mean(x, w)
  return(sqrt(weighted.mean((x - m)^2, w)) * (length(x) * (d + 
                                                             2)/4)^(-1/(d + 4)))
}


region.P <- function (region, P, k = 100) 
{
  vert0 = cartesian2polar2D(region, P)
  vert0.theta = vert0[, 2]
  ind.theta = sort.list(vert0.theta)
  vert0 = vert0[ind.theta, ]
  n.vert0 = nrow(vert0)
  vert = matrix(0, length(vert0.theta) + 1, 2)
  vert[1:n.vert0, ] = vert0
  vert[n.vert0 + 1, ] = vert0[1, ]
  vert.theta = vert[, 2]
  vert.rho = vert[, 1]
  n.vert = nrow(vert)
  vert.seq = 1:(n.vert - 1)
  vert.gamma = angle.positive(diff(vert.theta))
  vert.alpha = triangle.abgamma.alpha(vert.rho[vert.seq + 1], 
                                      vert.rho[vert.seq], vert.gamma)
  theta = angle.positive(seq(0, 2 * pi, length.out = k + 1) + 
                           vert.theta[1] - 2 * pi)
  theta = theta[2:(k + 1)]
  ind1 = findInterval(theta, vert.theta[vert.seq])
  ind1 = ifelse(ind1 == 0, max(vert.seq), ind1)
  alpha.inf = vert.alpha[ind1]
  theta.inf = angle.positive(theta - vert.theta[ind1])
  rho.inf = vert.rho[ind1]
  rho = triangle.alphabetac.b(theta.inf, alpha.inf, rho.inf)
  return(list(vert.theta = vert.theta, vert.rho = vert.rho, 
              vert.gamma = vert.gamma, vert.alpha = vert.alpha, theta = theta, 
              vert.seq = vert.seq, theta.inf = theta.inf, rho.inf = rho.inf, 
              alpha.inf = alpha.inf, rho = rho, ind1 = ind1))
}


cartesian2polar2D <- function (x, orig = c(0, 0)) 
{
  x = x - outer(rep(1, dim(x)[1]), orig)
  y = x
  y[, 2] = angle.positive(atan2(x[, 2], x[, 1]))
  y[, 1] = sqrt(x[, 2] * x[, 2] + x[, 1] * x[, 1])
  return(y)
}

angle.positive <- function (theta) 
  ifelse(theta < 0, theta + 2 * pi, theta)


triangle.abgamma.alpha <- function (a, b, gamma) 
{
  c = triangle.abgamma.c(a, b, gamma)
  triangle.abc.gamma(b, c, a)
}

triangle.abgamma.c <- function (a, b, gamma) 
  sqrt(a * a + b * b - 2 * a * b * cos(gamma))

triangle.abc.gamma <- function (a, b, c) 
  acos((a * a + b * b - c * c)/(2 * a * b))

triangle.alphabetac.b <- function (alpha, beta, c) 
  c * sin(beta)/sin(alpha + beta)

#######################################################################
# 2) FUNZIONI ETAS PARAMETRICHE E SIMULAZIONE RISCRITTE/PATCHATE
#######################################################################

WPI_bkgd.fit <- function(
    X,
    formula,
    covs = NULL,
    marked = FALSE,
    spatial.cov = FALSE,
    verbose = TRUE,
    mult = 4,
    seed = NULL,
    ncube = NULL,
    grid = FALSE,
    mark.c = FALSE,
    mark.c.mode = c("idw", "mc"),
    process.type = NULL,
    type.cov.values = NULL,
    Klocal = FALSE,
    interp.p = 81,
    pb = NULL
) {
  
  if (!inherits(X, c("ppp", "stp", "stpm", "s3dp", "s3dpm"))) {
    stop("X should be one of classes: 'ppp','stp','stpm','s3dp','s3dpm'")
  }
  
  if (is.null(process.type)) {
    process.type <- ifelse(
      inherits(X, c("s3dp", "s3dpm")),
      "s3d",
      ifelse(inherits(X, "ppp"), "s2d", "st")
    )
  } else {
    process.type <- match.arg(tolower(process.type), c("s2d", "st", "s3d"))
  }
  
  third_name <- ifelse(
    process.type == "st",
    "t",
    ifelse(process.type == "s3d", "z", NA)
  )
  
  time1 <- Sys.time()
  
  if (verbose) {
    cat("Starting stppm.WPI (process.type = ", process.type, ")", "\n\n")
  }
  
  if (!is.numeric(mult) || mult <= 0) {
    stop("'mult' must be a positive numeric value")
  }
  
  if (!is.null(ncube)) {
    if (!is.numeric(ncube) || ncube <= 0) {
      stop("'ncube' must be a positive numeric value")
    }
  }
  
  X0 <- X
  
  if (process.type == "s2d") {
    
    if (inherits(X0, "ppp")) {
      
      X <- data.frame(
        x = X0$x,
        y = X0$y,
        check.names = FALSE
      )
      
      ## IMPORTANTISSIMO: recupera eventuali marks dal ppp
      marks_df <- extract_ppp_marks_for_WPI(X0, formula = formula)
      
      if (!is.null(marks_df)) {
        if (nrow(marks_df) != nrow(X)) {
          stop(
            "Internal error: ppp marks have wrong length. ",
            "nrow(marks_df) = ", nrow(marks_df),
            ", nrow(X) = ", nrow(X)
          )
        }
        
        X <- cbind(X, marks_df)
      }
      
    } else if (inherits(X0, c("stp", "stpm", "s3dp", "s3dpm"))) {
      
      Xdf <- as.data.frame(X0$df, check.names = FALSE)
      
      if (ncol(Xdf) < 2) {
        stop("For process.type = 's2d', X0$df must contain at least two coordinate columns.")
      }
      
      names(Xdf)[1:2] <- c("x", "y")
      X <- Xdf
      
    } else {
      
      stop("For process.type = 's2d', X must be 'ppp' or include at least x,y in $df")
    }
    
    ## Forza a factor le marks categoriali presenti nella formula
    vars_formula <- all.vars(formula)
    possible_mark_vars <- setdiff(vars_formula, c("x", "y", "t", "z"))
    possible_mark_vars <- intersect(possible_mark_vars, names(X))
    
    for (nm in possible_mark_vars) {
      if (is.character(X[[nm]]) || is.factor(X[[nm]])) {
        X[[nm]] <- factor(X[[nm]])
      }
    }
    
    nX <- nrow(X)
    x <- X[, 1]
    y <- X[, 2]
    
    if (inherits(X0, "ppp") && !is.null(X0$window)) {
      s.region <- matrix(
        c(
          X0$window$xrange[1], X0$window$yrange[1],
          X0$window$xrange[2], X0$window$yrange[1],
          X0$window$xrange[1], X0$window$yrange[2],
          X0$window$xrange[2], X0$window$yrange[2]
        ),
        ncol = 2,
        byrow = TRUE
      )
    } else {
      s.region <- splancs::sbox(cbind(x, y), xfrac = 0, yfrac = 0)
    }
    
    if (verbose) {
      cat("Built 2D spatial window (s.region)\n\n")
    }
    
  }
  
  if (verbose) {
    cat("Observed points: nX = ", nX, "\n\n")
  }
  
  if (!is.null(pb)) {
    pb <- as.numeric(pb)
    if (length(pb) != nX) {
      stop(
        "Length of 'pb' must match the number of observed points. ",
        "length(pb) = ", length(pb), ", nX = ", nX
      )
    }
    if (any(!is.finite(pb))) {
      stop("'pb' contains non-finite values.")
    }
  }
  
  HomLambda <- nX
  rho <- mult * HomLambda
  
  if (grid) {
    
    if (verbose) {
      cat("Generating dummy points on a grid...", "\n\n")
    }
    
    if (process.type == "s2d") {
      
      ff <- max(1L, floor(rho^(1 / 2)))
      x0 <- y0 <- seq_len(ff)
      
      x0 <- scale_to_range(x0, s.region[1, 1], s.region[2, 1])
      y0 <- scale_to_range(y0, s.region[1, 2], s.region[3, 2])
      
      df0 <- expand.grid(x0, y0)
      colnames(df0) <- c("x", "y")
      
      dummy_points <- df0
      
    } else {
      
      ff <- max(1L, floor(rho^(1 / 3)))
      x0 <- y0 <- d30 <- seq_len(ff)
      
      x0 <- scale_to_range(x0, s.region[1, 1], s.region[2, 1])
      y0 <- scale_to_range(y0, s.region[1, 2], s.region[3, 2])
      d30 <- scale_to_range(d30, third.region[1], third.region[2])
      
      df0 <- expand.grid(x0, y0, d30)
      colnames(df0) <- c("x", "y", third_name)
      
      dummy_points <- if (process.type == "st") {
        stp(cbind(df0$x, df0$y, df0[[third_name]]))$df
      } else {
        s3dp(cbind(df0$x, df0$y, df0[[third_name]]))$df
      }
    }
    
    if (verbose) {
      cat("Dummy points (grid): ", nrow(dummy_points), "\n\n")
    }
    
  } else {
    
    if (verbose) {
      cat("Generating dummy points at random...", "\n\n")
    }
    
    if (!is.null(seed)) {
      set.seed(seed)
      if (verbose) {
        cat("Random seed = ", as.character(seed), "\n")
      }
    }
    
    if (process.type == "s2d") {
      
      dummy_points <- rstpp.2D(
        lambda = rho,
        nsim = 1,
        verbose = FALSE,
        minX = s.region[1, 1],
        maxX = s.region[2, 1],
        minY = s.region[1, 2],
        maxY = s.region[3, 2]
      )
      
      dummy_points <- data.frame(
        x = dummy_points$x,
        y = dummy_points$y
      )
      
    } else {
      
      dummy_points <- r_st_s3dpp(
        lambda = rho,
        nsim = 1,
        verbose = FALSE,
        process.type = process.type,
        minX = s.region[1, 1],
        maxX = s.region[2, 1],
        minY = s.region[1, 2],
        maxY = s.region[3, 2],
        minT = third.region[1],
        maxT = third.region[2]
      )$points
    }
    
    if (verbose) {
      cat("Dummy points (random): ", nrow(dummy_points), "\n\n")
    }
  }
  
  if (process.type == "s2d") {
    
    quad_p <- rbind(
      as.matrix(X[, 1:2]),
      as.matrix(dummy_points[, 1:2])
    )
    
    colnames(quad_p) <- c("x", "y")
    
    if (verbose) {
      cat(
        "Quadrature points assembled (2D): total (data + dummy) = ",
        nrow(quad_p), "\n\n"
      )
    }
    
    xx <- quad_p[, 1]
    xy <- quad_p[, 2]
    
    win <- spatstat.geom::owin(
      xrange = c(s.region[1, 1], s.region[2, 1]),
      yrange = c(s.region[1, 2], s.region[3, 2])
    )
    
    if (is.null(ncube)) {
      ncube <- default.ncube.2D(quad_p)
    }
    
    ncube <- rep.int(ncube, 2)
    nx <- ncube[1]
    ny <- ncube[2]
    
    nxy <- nx * ny
    
    cubearea <- spatstat.geom::area.owin(win) / nxy
    volumes <- rep.int(cubearea, nxy)
    
    id <- grid.index.2D(
      xx,
      xy,
      win$xrange,
      win$yrange,
      nx,
      ny
    )$index
    
    w <- counting.weights(id, volumes)
    
  } else {
    
    quad_p <- rbind(
      as.matrix(X[, 1:3]),
      as.matrix(dummy_points[, 1:3])
    )
    
    colnames(quad_p) <- c("x", "y", third_name)
    
    if (verbose) {
      cat(
        "Quadrature points assembled (3D): total (data + dummy) = ",
        nrow(quad_p), "\n\n"
      )
    }
    
    xx <- quad_p[, 1]
    xy <- quad_p[, 2]
    xd3 <- quad_p[, 3]
    
    win <- spatstat.geom::box3(
      xrange = range(xx, na.rm = TRUE),
      yrange = range(xy, na.rm = TRUE),
      zrange = range(xd3, na.rm = TRUE)
    )
    
    if (is.null(ncube)) {
      ncube <- default.ncube(quad_p)
    }
    
    ncube <- rep.int(ncube, 3)
    nx <- ncube[1]
    ny <- ncube[2]
    nt <- ncube[3]
    
    nxyt <- nx * ny * nt
    
    cubevolume <- spatstat.geom::volume(win) / nxyt
    volumes <- rep.int(cubevolume, nxyt)
    
    id <- grid.index(
      xx,
      xy,
      xd3,
      win$xrange,
      win$yrange,
      win$zrange,
      nx,
      ny,
      nt
    )$index
    
    w <- counting.weights(id, volumes)
  }
  
  ndata <- nrow(X)
  ndummy <- nrow(dummy_points)
  
  Wdat <- w[seq_len(ndata)]
  Wdum <- w[(ndata + 1):(ndata + ndummy)]
  
  weights_gam <- w
  
  if (marked == TRUE && spatial.cov == TRUE) {
    
    if (verbose) {
      cat("Case: marked = TRUE, spatial.cov = TRUE. Building external covariates on quadrature...", "\n\n")
    }
    
    points.covs <- build_points_covs(
      quad_p = quad_p,
      covs = covs,
      formula = formula,
      type.cov.values = type.cov.values,
      interp.p = interp.p,
      process.type = process.type
    )
    
    coord_cols <- if (process.type == "s2d") 1:2 else 1:3
    
    if (verbose) {
      cat(
        "External covariates ready: ",
        max(0, ncol(points.covs) - length(coord_cols)),
        "\n\n"
      )
    }
    
    if (verbose) {
      cat("Building replicated quadrature for categorical marks...", "\n\n")
    }
    
    marked.process <- replicated.cubature.dummy(
      X = X,
      formula = formula,
      dummy_points = dummy_points,
      Wdum = Wdum,
      Wdat = Wdat,
      ndata = ndata,
      ndummy = ndummy,
      process.type = process.type
    )
    
    df_dumb <- as.data.frame(marked.process$dumb)
    n_dumb <- nrow(df_dumb)
    n_comb <- marked.process$n_comb_levels
    
    if (verbose) {
      cat(
        "Replicated quadrature: combinations = ",
        n_comb,
        " replicated dummy rows = ",
        n_dumb,
        "\n\n"
      )
    }
    
    z <- if (is.null(pb)) {
      c(rep(1, ndata), rep(0, length(marked.process$Wdumb)))
    } else {
      c(pb, rep(0, length(marked.process$Wdumb)))
    }
    
    w_final <- c(w[seq_len(ndata)], marked.process$Wdumb)
    y_resp <- z / w_final
    
    if (process.type == "s2d") {
      
      dati.modello <- data.frame(
        y_resp = y_resp,
        w = w_final,
        x = c(X$x, df_dumb$x),
        y = c(X$y, df_dumb$y),
        check.names = FALSE
      )
      
    } else {
      
      third_name <- if (process.type == "st") "t" else "z"
      
      dati.modello <- data.frame(
        y_resp = y_resp,
        w = w_final,
        x = c(X$x, df_dumb$x),
        y = c(X$y, df_dumb$y),
        third = c(X[[third_name]], df_dumb$z),
        check.names = FALSE
      )
      
      names(dati.modello)[names(dati.modello) == "third"] <- third_name
    }
    
    dati.modello <- cbind(
      dati.modello,
      rbind(marked.process$df_marks, marked.process$total_dummy_marks)
    )
    
    stopifnot(nrow(dati.modello) == (ndata + ndummy) * n_comb)
    
    if (verbose) {
      cat(
        "Model frame (with replicated marks) assembled: ",
        nrow(dati.modello),
        " rows.",
        "\n\n"
      )
    }
    
    cov_names_only <- names(points.covs)[-coord_cols]
    
    if (length(cov_names_only)) {
      
      if (verbose) {
        cat("Replicating external covariates across mark combinations...", "\n\n")
      }
      
      dati.interpolati.rep <- as.data.frame(
        matrix(
          NA_real_,
          nrow = nrow(dati.modello),
          ncol = length(cov_names_only)
        ),
        check.names = FALSE
      )
      
      colnames(dati.interpolati.rep) <- cov_names_only
      
      stopifnot(nrow(points.covs) == (ndata + ndummy))
      
      for (nm in cov_names_only) {
        v_data <- points.covs[seq_len(ndata), nm]
        v_dummy <- points.covs[(ndata + 1):(ndata + ndummy), nm]
        
        dati.interpolati.rep[[nm]] <- c(
          v_data,
          rep(v_dummy, each = n_comb),
          if (n_comb > 1) rep(v_data, each = n_comb - 1) else NULL
        )
      }
      
      dati.cov.marks <- cbind(dati.modello, dati.interpolati.rep)
      
    } else {
      dati.cov.marks <- dati.modello
    }
    
    weights_gam <- w_final
    
    if (verbose) {
      cat(
        "Design matrix ready (rows: ",
        nrow(dati.cov.marks),
        ", cols: ",
        ncol(dati.cov.marks),
        ")\n\n"
      )
    }
    
  } else if (marked == FALSE && spatial.cov == TRUE) {
    
    if (verbose) {
      cat("Case: marked = FALSE, spatial.cov = TRUE. Building external covariates on quadrature...", "\n\n")
    }
    
    points.covs <- build_points_covs(
      quad_p = quad_p,
      covs = covs,
      formula = formula,
      type.cov.values = type.cov.values,
      interp.p = interp.p,
      process.type = process.type
    )
    
    z <- if (is.null(pb)) {
      c(rep(1, ndata), rep(0, ndummy))
    } else {
      c(pb, rep(0, ndummy))
    }
    
    y_resp <- z / w
    
    dati.cov.marks <- cbind(
      y_resp = y_resp,
      w = w,
      points.covs
    )
    
    weights_gam <- w
    
    if (verbose) {
      cat(
        "Design matrix ready (rows: ",
        nrow(dati.cov.marks),
        ", cols: ",
        ncol(dati.cov.marks),
        ")\n\n"
      )
    }
    
  } else if (marked == FALSE && spatial.cov == FALSE) {
    
    if (verbose) {
      cat("Case: marked = FALSE, spatial.cov = FALSE. Using only coordinates.\n")
    }
    
    if (process.type == "s2d") {
      colnames(quad_p) <- c("x", "y")
    } else {
      third_name <- if (process.type == "st") "t" else "z"
      colnames(quad_p) <- c("x", "y", third_name)
    }
    
    z <- if (is.null(pb)) {
      c(rep(1, ndata), rep(0, ndummy))
    } else {
      c(pb, rep(0, ndummy))
    }
    
    y_resp <- z / w
    
    dati.cov.marks <- cbind(
      y_resp = y_resp,
      w = w,
      quad_p
    )
    
    weights_gam <- w
    
    if (verbose) {
      cat(
        "Design matrix ready (rows: ",
        nrow(dati.cov.marks),
        ", cols: ",
        ncol(dati.cov.marks),
        ")\n\n"
      )
    }
    
  } else if (marked == TRUE && spatial.cov == FALSE) {
    
    if (verbose) {
      cat("Case: marked = TRUE, spatial.cov = FALSE. Building replicated quadrature for categorical marks...\n\n")
    }
    
    mark_vars_formula <- setdiff(all.vars(formula), c("x", "y", "t", "z"))
    mark_vars_formula <- intersect(mark_vars_formula, all.vars(formula))
    
    if (!length(mark_vars_formula)) {
      stop("marked = TRUE but no mark variable was found in the formula.")
    }
    
    missing_marks <- setdiff(mark_vars_formula, names(X))
    
    if (length(missing_marks)) {
      stop(
        "marked = TRUE but the following mark variables are missing from X: ",
        paste(missing_marks, collapse = ", "),
        ". This usually means that marks were lost when converting the ppp object."
      )
    }
    
    for (nm in mark_vars_formula) {
      X[[nm]] <- factor(X[[nm]])
    }
    
    if (verbose) {
      cat("Categorical mark variables found in X: ",
          paste(mark_vars_formula, collapse = ", "), "\n\n")
      print(lapply(X[mark_vars_formula], table))
    }
    
    marked.process <- replicated.cubature.dummy(
      X = X,
      formula = formula,
      dummy_points = dummy_points,
      Wdum = Wdum,
      Wdat = Wdat,
      ndata = ndata,
      ndummy = ndummy,
      process.type = process.type
    )
    
    df_dumb <- as.data.frame(marked.process$dumb)
    n_dumb <- nrow(df_dumb)
    
    z <- if (is.null(pb)) {
      c(rep(1, ndata), rep(0, length(marked.process$Wdumb)))
    } else {
      c(pb, rep(0, length(marked.process$Wdumb)))
    }
    
    w_final <- c(w[seq_len(ndata)], marked.process$Wdumb)
    y_resp <- z / w_final
    
    if (process.type == "s2d") {
      
      dati.modello <- data.frame(
        y_resp = y_resp,
        w = w_final,
        x = c(X$x, df_dumb$x),
        y = c(X$y, df_dumb$y),
        check.names = FALSE
      )
      
    } else {
      
      third_name <- if (process.type == "st") "t" else "z"
      
      dati.modello <- data.frame(
        y_resp = y_resp,
        w = w_final,
        x = c(X$x, df_dumb$x),
        y = c(X$y, df_dumb$y),
        third = c(X[[third_name]], df_dumb$z),
        check.names = FALSE
      )
      
      names(dati.modello)[names(dati.modello) == "third"] <- third_name
    }
    
    if (ncol(marked.process$df_marks) > 0) {
      dati.modello <- cbind(
        dati.modello,
        rbind(marked.process$df_marks, marked.process$total_dummy_marks)
      )
    }
    
    dati.cov.marks <- dati.modello
    weights_gam <- w_final
    
    if (verbose) {
      cat(
        "Design matrix ready (rows: ",
        nrow(dati.cov.marks),
        ", cols: ",
        ncol(dati.cov.marks),
        ")\n\n"
      )
    }
  }
  
  if (mark.c) {
    
    if (is.null(mark.c.mode)) {
      mark.c.mode <- "idw"
    } else {
      mark.c.mode <- match.arg(mark.c.mode, c("idw", "mc"))
    }
    
    if (isTRUE(verbose)) {
      cat("Handling continuous marks in the model formula (mode = ", mark.c.mode, ")...\n\n")
    }
    
    if (process.type == "s2d") {
      coord_cols <- c("x", "y")
    } else {
      third_name <- if (process.type == "st") "t" else "z"
      coord_cols <- c("x", "y", third_name)
    }
    
    form_vars <- all.vars(formula)
    mark_candidates <- setdiff(colnames(X), coord_cols)
    
    numeric_marks <- if (length(mark_candidates)) {
      mark_candidates[
        vapply(X[, mark_candidates, drop = FALSE], is.numeric, logical(1))
      ]
    } else {
      character(0)
    }
    
    cont_mark_names <- intersect(form_vars, numeric_marks)
    
    if (!length(cont_mark_names)) {
      
      if (isTRUE(verbose)) {
        cat("No continuous marks found in the formula; skipping.\n\n")
      }
      
    } else {
      
      n_base <- ndata + ndummy
      n_total <- nrow(dati.cov.marks)
      
      if (n_total %% n_base != 0) {
        stop("Quadrature replication mismatch: nrow(dati.cov.marks) is not a multiple of (ndata + ndummy).")
      }
      
      n_comb <- as.integer(n_total / n_base)
      
      if (identical(mark.c.mode, "idw")) {
        
        df_contmarks <- as.data.frame(
          matrix(
            NA_real_,
            nrow = n_total,
            ncol = length(cont_mark_names)
          ),
          check.names = FALSE
        )
        
        colnames(df_contmarks) <- cont_mark_names
        
        for (nm in cont_mark_names) {
          
          covs_mark <- X[, c(coord_cols, nm), drop = FALSE]
          
          pred_dummy <- idw_interp_fast(
            points = dummy_points[, coord_cols, drop = FALSE],
            covs = covs_mark,
            p = interp.p,
            d = length(coord_cols)
          )
          
          if (n_comb > 1) {
            vec_extended <- c(
              X[[nm]],
              rep(pred_dummy, each = n_comb),
              rep(X[[nm]], each = n_comb - 1)
            )
          } else {
            vec_extended <- c(X[[nm]], pred_dummy)
          }
          
          df_contmarks[[nm]] <- vec_extended
        }
        
        dati.cov.marks <- cbind(dati.cov.marks, df_contmarks)
        
        if (isTRUE(verbose)) {
          cat("Continuous marks appended (IDW). cols added:", ncol(df_contmarks), "\n\n")
        }
        
      } else {
        
        if (isTRUE(verbose)) {
          cat("Replicating dummy rows and sampling continuous marks (Monte Carlo)...\n\n")
        }
        
        R <- if (!is.null(getOption("mark_mc_reps"))) {
          getOption("mark_mc_reps")
        } else {
          100L
        }
        
        if (!is.numeric(R) || R < 1) {
          R <- 100L
        }
        
        data_block <- dati.cov.marks[seq_len(ndata), , drop = FALSE]
        dummy_block <- dati.cov.marks[(ndata + 1):nrow(dati.cov.marks), , drop = FALSE]
        nD <- nrow(dummy_block)
        
        if (nD != ndummy * n_comb) {
          stop("Internal error: dummy block length doesn't match ndummy * n_comb.")
        }
        
        idx_rep <- rep(seq_len(nD), each = R)
        dummy_rep <- dummy_block[idx_rep, , drop = FALSE]
        
        sampler_for <- function(nm) {
          rlo <- min(X[[nm]], na.rm = TRUE)
          rhi <- max(X[[nm]], na.rm = TRUE)
          
          if (!is.finite(rlo) || !is.finite(rhi) || rlo >= rhi) {
            stop(sprintf(
              "Cannot build sampler for mark '%s' (range not finite or degenerate).",
              nm
            ))
          }
          
          function(n) stats::runif(n, rlo, rhi)
        }
        
        M_rep <- as.data.frame(
          setNames(
            lapply(cont_mark_names, function(nm) sampler_for(nm)(nrow(dummy_rep))),
            cont_mark_names
          ),
          check.names = FALSE
        )
        
        if (!("w" %in% names(dummy_rep))) {
          stop("Internal error: 'w' column not found in dummy block.")
        }
        
        dummy_rep$w <- dummy_rep$w / R
        
        dummy_block_mc <- cbind(dummy_rep, M_rep)
        data_block_mc <- cbind(data_block, X[, cont_mark_names, drop = FALSE])
        
        dati.cov.marks <- rbind(data_block_mc, dummy_block_mc)
        
        z <- if (is.null(pb)) {
          c(rep(1, ndata), rep(0, nrow(dummy_block_mc)))
        } else {
          c(pb, rep(0, nrow(dummy_block_mc)))
        }
        
        weights_gam <- c(data_block_mc$w, dummy_block_mc$w)
        y_resp <- z / weights_gam
        
        dati.cov.marks$w <- weights_gam
        dati.cov.marks$y_resp <- y_resp
        
        if (isTRUE(verbose)) {
          cat(
            "Monte Carlo reps =",
            R,
            " -> dummy rows:",
            nD,
            "->",
            nrow(dummy_block_mc),
            "\n\n"
          )
        }
      }
    }
  }
  
  if (isTRUE(Klocal)) {
    
    if (verbose) {
      cat("Computing local K-function offset...", "\n\n")
    }
    
    third_name <- if (process.type == "st") "t" else "z"
    
    if (process.type == "st") {
      
      obj_st <- stp(as.matrix(X[, c("x", "y", "t")]))
      kout <- Khat_st(obj_st, correction = "translate", Klocal = TRUE)
      
      Kloc <- kout$Khat
      Ktheo <- kout$Ktheo
      
      npt <- dim(Kloc)[3]
      
      w_obs <- vapply(seq_len(npt), function(i) {
        sum((Kloc[, , i] - Ktheo) / Ktheo, na.rm = TRUE)
      }, numeric(1))
      
    } else {
      
      obj_s3d <- s3dp(as.matrix(X[, c("x", "y", "z")]))
      kout <- Khat_s3d(obj_s3d, correction = "translate", Klocal = TRUE)
      
      Kloc <- kout$Khat
      Ktheo <- kout$Ktheo
      
      rel <- sweep(Kloc, 1, Ktheo, FUN = function(a, b) (a - b) / b)
      w_obs <- colSums(rel, na.rm = TRUE)
    }
    
    X_Klocal <- X[, c("x", "y", third_name), drop = FALSE]
    X_Klocal$w_klocal <- w_obs
    
    kpred <- idw_interp_fast(
      points = dummy_points[, c("x", "y", third_name), drop = FALSE],
      covs = X_Klocal[, c("x", "y", third_name, "w_klocal"), drop = FALSE],
      p = interp.p,
      d = 3
    )
    
    if (isTRUE(marked) && exists("marked.process") && !is.null(marked.process$n_comb_levels)) {
      m <- marked.process$n_comb_levels
      
      offset_vec <- c(
        w_obs,
        rep(kpred, m),
        rep(w_obs, m - 1L)
      )
      
    } else {
      offset_vec <- c(w_obs, kpred)
    }
    
    offset_vec <- offset_vec - mean(offset_vec, na.rm = TRUE)
    dati.cov.marks$str_Klocal <- offset_vec
    
    if (verbose) {
      cat("Local K offset computed and appended.", "\n\n")
    }
  }
  
  dati.cov.marks <- as.data.frame(dati.cov.marks)
  
  if (length(weights_gam) != nrow(dati.cov.marks)) {
    stop(
      "Internal error: length(weights_gam) does not match nrow(dati.cov.marks). ",
      "length(weights_gam) = ", length(weights_gam),
      ", nrow(dati.cov.marks) = ", nrow(dati.cov.marks)
    )
  }
  
  dati.cov.marks$w <- as.numeric(weights_gam)
  weights_gam <- dati.cov.marks$w
  
  if (!("y_resp" %in% names(dati.cov.marks))) {
    stop("Internal error: response variable 'y_resp' not found.")
  }
  
  y_resp <- dati.cov.marks$y_resp
  
  if (length(y_resp) != nrow(dati.cov.marks)) {
    stop("Internal error: length(y_resp) does not match model frame.")
  }
  
  if (verbose) {
    print(summary(dati.cov.marks))
  }
  
  add_off <- isTRUE(Klocal) && "str_Klocal" %in% names(dati.cov.marks)
  form_full <- build_poisson_formula(formula, add_offset = add_off)
  
  suppressWarnings(
    mod_global <- try(
      mgcv::gam(
        formula = form_full,
        family = poisson(link = "log"),
        data = dati.cov.marks,
        weights = w
      ),
      silent = TRUE
    )
  )
  
  if (inherits(mod_global, "try-error")) {
    stop("Model fit failed: ", as.character(mod_global))
  }
  
  if (verbose) {
    cat("Model fit complete.", "\n\n")
  }
  
  pred_global <- predict(
    mod_global,
    newdata = dati.cov.marks[seq_len(nX), , drop = FALSE],
    type = "response"
  )
  
  pred_global_dummy <- predict(
    mod_global,
    newdata = dati.cov.marks[(nX + 1):nrow(dati.cov.marks), , drop = FALSE],
    type = "response"
  )
  
  res_global <- coef(mod_global)
  
  time2 <- Sys.time()
  elapsed <- round(as.numeric(difftime(time2, time1, units = "sec")), 3)
  
  if (verbose) {
    cat("Done in ", elapsed, " secs.", "\n\n")
  }
  
  list.obj <- list(
    IntCoefs = res_global,
    X = X0,
    nX = ndata,
    I = z,
    y_resp = y_resp,
    formula = formula,
    l = as.vector(pred_global),
    l_dummy = as.vector(pred_global_dummy),
    mod_global = mod_global,
    quad.data = dati.cov.marks,
    newdata = dati.cov.marks[seq_len(ndata), , drop = FALSE],
    dummy.data = dati.cov.marks[(ndata + 1):nrow(dati.cov.marks), , drop = FALSE],
    ncube = ncube,
    weights = weights_gam,
    time = paste0(elapsed, " sec")
  )
  
  class(list.obj) <- switch(
    process.type,
    s2d = "s2dppm",
    st = "stppm",
    s3d = "s3dppm"
  )
  
  return(list.obj)
}









replicated.cubature.dummy <- function(
    X,
    formula,
    dummy_points,
    Wdum,
    Wdat,
    ndata,
    ndummy,
    process.type = c("s2d", "st", "s3d")
) {
  
  process.type <- match.arg(tolower(process.type), c("s2d", "st", "s3d"))
  
  ############################################################
  ## 0. Input coercion and coordinate names
  ############################################################
  
  X <- as.data.frame(X, check.names = FALSE)
  dummy_points <- as.data.frame(dummy_points, check.names = FALSE)
  
  if (process.type == "s2d") {
    
    coord_cols <- c("x", "y")
    
    if (!all(coord_cols %in% names(X))) {
      if (ncol(X) < 2) {
        stop("For process.type = 's2d', X must contain at least two coordinate columns.")
      }
      names(X)[1:2] <- coord_cols
    }
    
    if (!all(coord_cols %in% names(dummy_points))) {
      if (ncol(dummy_points) < 2) {
        stop("For process.type = 's2d', dummy_points must contain at least two coordinate columns.")
      }
      names(dummy_points)[1:2] <- coord_cols
    }
    
  } else {
    
    third_name <- if (process.type == "st") "t" else "z"
    coord_cols <- c("x", "y", third_name)
    
    if (!all(coord_cols %in% names(X))) {
      if (ncol(X) < 3) {
        stop("For process.type = 'st' or 's3d', X must contain at least three coordinate columns.")
      }
      names(X)[1:3] <- coord_cols
    }
    
    if (!all(coord_cols %in% names(dummy_points))) {
      if (ncol(dummy_points) < 3) {
        stop("For process.type = 'st' or 's3d', dummy_points must contain at least three coordinate columns.")
      }
      names(dummy_points)[1:3] <- coord_cols
    }
  }
  
  if (nrow(X) != ndata) {
    stop(
      "ndata does not match nrow(X). ",
      "ndata = ", ndata, ", nrow(X) = ", nrow(X)
    )
  }
  
  if (nrow(dummy_points) != ndummy) {
    stop(
      "ndummy does not match nrow(dummy_points). ",
      "ndummy = ", ndummy, ", nrow(dummy_points) = ", nrow(dummy_points)
    )
  }
  
  if (length(Wdat) != ndata) {
    stop(
      "Length of Wdat must equal ndata. ",
      "length(Wdat) = ", length(Wdat), ", ndata = ", ndata
    )
  }
  
  if (length(Wdum) != ndummy) {
    stop(
      "Length of Wdum must equal ndummy. ",
      "length(Wdum) = ", length(Wdum), ", ndummy = ", ndummy
    )
  }
  
  ############################################################
  ## 1. Identify categorical marks in the formula
  ############################################################
  
  vars_formula <- all.vars(formula)
  
  marks.name.formula <- setdiff(
    vars_formula,
    c("x", "y", "t", "z")
  )
  
  if (!length(marks.name.formula)) {
    stop("No mark variable found in the formula.")
  }
  
  missing_marks <- setdiff(marks.name.formula, names(X))
  
  if (length(missing_marks)) {
    stop(
      "The following mark variables are in the formula but not in X: ",
      paste(missing_marks, collapse = ", "),
      ". Available columns are: ",
      paste(names(X), collapse = ", ")
    )
  }
  
  ############################################################
  ## 2. Force categorical marks to factor and collect levels
  ############################################################
  
  df_marks <- X[, marks.name.formula, drop = FALSE]
  
  for (nm in marks.name.formula) {
    
    if (is.factor(df_marks[[nm]])) {
      lev_nm <- levels(df_marks[[nm]])
    } else {
      lev_nm <- sort(unique(as.character(df_marks[[nm]])))
      df_marks[[nm]] <- factor(as.character(df_marks[[nm]]), levels = lev_nm)
    }
    
    X[[nm]] <- factor(as.character(X[[nm]]), levels = lev_nm)
    df_marks[[nm]] <- factor(as.character(df_marks[[nm]]), levels = lev_nm)
  }
  
  list_levels_marks <- lapply(df_marks, levels)
  names(list_levels_marks) <- marks.name.formula
  
  if (any(vapply(list_levels_marks, length, integer(1)) < 1)) {
    stop("At least one categorical mark has no levels.")
  }
  
  ############################################################
  ## 3. Build all combinations of categorical mark levels
  ############################################################
  
  levels_df <- do.call(
    expand.grid,
    c(
      list_levels_marks,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  
  names(levels_df) <- marks.name.formula
  
  for (nm in marks.name.formula) {
    levels_df[[nm]] <- factor(
      levels_df[[nm]],
      levels = list_levels_marks[[nm]]
    )
  }
  
  n_comb_levels <- nrow(levels_df)
  
  ############################################################
  ## 4. Helper: product between spatial points and mark levels
  ############################################################
  
  make_product_df <- function(points_df, weights_vec) {
    
    points_df <- as.data.frame(points_df, check.names = FALSE)
    n_pts <- nrow(points_df)
    
    coords_rep <- points_df[
      rep(seq_len(n_pts), each = n_comb_levels),
      coord_cols,
      drop = FALSE
    ]
    
    rownames(coords_rep) <- NULL
    
    marks_rep <- levels_df[
      rep(seq_len(n_comb_levels), times = n_pts),
      ,
      drop = FALSE
    ]
    
    rownames(marks_rep) <- NULL
    
    weights_rep <- rep(weights_vec, each = n_comb_levels)
    
    out <- cbind(coords_rep, marks_rep)
    
    list(
      df = out,
      weights = weights_rep
    )
  }
  
  ############################################################
  ## 5. Dummy spatial points × all mark levels
  ############################################################
  
  prod_dummy <- make_product_df(
    points_df = dummy_points[, coord_cols, drop = FALSE],
    weights_vec = Wdum
  )
  
  df_dumdum <- prod_dummy$df
  Wdumdum <- prod_dummy$weights
  
  ############################################################
  ## 6. Observed spatial locations × all mark levels
  ##    Then remove the actually observed combinations
  ############################################################
  
  prod_data <- make_product_df(
    points_df = X[, coord_cols, drop = FALSE],
    weights_vec = Wdat
  )
  
  df_dumdat_all <- prod_data$df
  Wdumdat_all <- prod_data$weights
  
  observed_product_df <- cbind(
    X[, coord_cols, drop = FALSE],
    df_marks
  )
  
  key_fun <- function(df, cols) {
    
    tmp <- df[, cols, drop = FALSE]
    
    for (cc in cols) {
      if (is.numeric(tmp[[cc]])) {
        tmp[[cc]] <- format(tmp[[cc]], digits = 17, scientific = FALSE)
      } else {
        tmp[[cc]] <- as.character(tmp[[cc]])
      }
    }
    
    do.call(paste, c(tmp, sep = "___"))
  }
  
  key_cols <- c(coord_cols, marks.name.formula)
  
  key_dumdat <- key_fun(df_dumdat_all, key_cols)
  key_obs <- key_fun(observed_product_df, key_cols)
  
  ## Remove exactly one replicated row for each observed point-mark combination.
  ## This avoids problems when two observed events share same coordinates and mark.
  remove_idx <- logical(length(key_dumdat))
  
  obs_tab <- table(key_obs)
  
  for (kk in names(obs_tab)) {
    
    idx_kk <- which(key_dumdat == kk)
    
    if (length(idx_kk)) {
      remove_idx[idx_kk[seq_len(min(length(idx_kk), obs_tab[[kk]]))]] <- TRUE
    }
  }
  
  df_dumdat <- df_dumdat_all[!remove_idx, , drop = FALSE]
  Wdumdat <- Wdumdat_all[!remove_idx]
  
  ############################################################
  ## 7. Combine dummy blocks
  ############################################################
  
  total_dummy_df <- rbind(
    df_dumdum,
    df_dumdat
  )
  
  Wdumb <- c(
    Wdumdum,
    Wdumdat
  )
  
  total_dummy_marks <- total_dummy_df[, marks.name.formula, drop = FALSE]
  
  for (nm in marks.name.formula) {
    total_dummy_marks[[nm]] <- factor(
      as.character(total_dummy_marks[[nm]]),
      levels = list_levels_marks[[nm]]
    )
  }
  
  ############################################################
  ## 8. Object returned as 'dumb'
  ##    It is a data.frame, which is enough because WPI_bkgd.fit
  ##    uses as.data.frame(marked.process$dumb).
  ############################################################
  
  dumb <- total_dummy_df[, c(coord_cols, marks.name.formula), drop = FALSE]
  
  ############################################################
  ## 9. Checks
  ############################################################
  
  if (nrow(dumb) != length(Wdumb)) {
    stop(
      "Internal error: nrow(dumb) != length(Wdumb). ",
      "nrow(dumb) = ", nrow(dumb),
      ", length(Wdumb) = ", length(Wdumb)
    )
  }
  
  if (nrow(total_dummy_marks) != length(Wdumb)) {
    stop(
      "Internal error: nrow(total_dummy_marks) != length(Wdumb). ",
      "nrow(total_dummy_marks) = ", nrow(total_dummy_marks),
      ", length(Wdumb) = ", length(Wdumb)
    )
  }
  
  if (nrow(df_marks) != ndata) {
    stop(
      "Internal error: nrow(df_marks) != ndata. ",
      "nrow(df_marks) = ", nrow(df_marks),
      ", ndata = ", ndata
    )
  }
  
  expected_dummy_rows <- ndummy * n_comb_levels + ndata * (n_comb_levels - 1)
  
  if (nrow(dumb) != expected_dummy_rows) {
    warning(
      "Unexpected number of replicated dummy rows. ",
      "Observed = ", nrow(dumb),
      ", expected = ", expected_dummy_rows,
      ". This can happen if duplicated observed coordinate-mark combinations are present."
    )
  }
  
  ############################################################
  ## 10. Return
  ############################################################
  
  list(
    dumb = dumb,
    df_marks = df_marks,
    total_dummy_marks = total_dummy_marks,
    Wdumb = Wdumb,
    n_comb_levels = n_comb_levels,
    levels_df = levels_df,
    list_levels_marks = list_levels_marks,
    marks.name.formula = marks.name.formula,
    coord_cols = coord_cols
  )
}





build_points_covs <- function(quad_p, covs, formula, type.cov.values,
                              interp.p = 81, process.type = c("s2d","st","s3d")) {
  points.covs <- as.data.frame(quad_p)
  
  ## se process.type non è passato, prova ad inferirlo dai nomi delle colonne
  if (missing(process.type) || is.null(process.type)) {
    nms <- names(points.covs)
    if ("t" %in% nms) {
      process.type <- "st"
    } else if ("z" %in% nms) {
      process.type <- "s3d"
    } else {
      process.type <- "s2d"
    }
  } else {
    process.type <- match.arg(tolower(process.type), c("s2d","st","s3d"))
  }
  
  ## coordinate ammesse in base al tipo di processo
  allowed_coords <- switch(
    process.type,
    "s2d" = c("x","y"),
    "st"  = c("x","y","t"),
    "s3d" = c("x","y","z")
  )
  
  cov_names <- intersect(names(covs), all.vars(formula))
  if (!length(cov_names)) {
    stop("No external covariates found in 'covs' following the model formula.")
  }
  
  allowed_types <- c("min","max","interp","exact")
  
  coord_order <- function(nms) intersect(allowed_coords, nms)
  
  coerce_time_numeric <- function(df, cols) {
    if ("t" %in% cols && inherits(df$t, c("POSIXct","POSIXt","Date"))) {
      df$t <- as.numeric(df$t)
    }
    df
  }
  
  for (nm in cov_names) {
    choice <- type.cov.values[[nm]]
    if (is.null(choice) || !(choice %in% allowed_types)) {
      stop("type.cov.values for covariate '", nm, "' must be one of: ",
           paste(allowed_types, collapse = ", "))
    }
    
    if (choice %in% c("min","max")) {
      ## caso distanza minima/massima da covariata esterna
      dfc0 <- covs[[nm]]
      if (is.null(dfc0)) {
        stop("Covariate '", nm, "' not found in 'covs'.")
      }
      present <- coord_order(names(dfc0))
      if (length(present) == 0L || length(present) > length(allowed_coords)) {
        stop(sprintf(
          "Cov '%s': 1-%d coordinate columns needed in external covariate dataset with colnames in %s.",
          nm, length(allowed_coords), paste(allowed_coords, collapse = ", ")
        ))
      }
      if (!all(present %in% names(as.data.frame(quad_p)))) {
        miss <- setdiff(present, names(as.data.frame(quad_p)))
        stop(sprintf("Cov '%s': columns missing in data points: %s",
                     nm, paste(miss, collapse = ", ")))
      }
      dfp <- quad_p[, present, drop = FALSE]
      dfc <- dfc0[, present, drop = FALSE]
      dfp <- coerce_time_numeric(dfp, present)
      dfc <- coerce_time_numeric(dfc, present)
      
      D <- distN_fast(dfp, dfc, d = length(present), coord_names = present)
      agg_fun <- if (choice == "min") min else max
      vals <- apply(D, 1L, agg_fun, na.rm = TRUE)
      points.covs[[nm]] <- vals
      
    } else if (choice == "interp") {
      ## caso interpolazione (IDW)
      dfc0 <- covs[[nm]]
      if (is.null(dfc0)) {
        stop("Covariate '", nm, "' not found in 'covs'.")
      }
      
      coord_p <- coord_order(names(as.data.frame(quad_p)))
      if (!length(coord_p)) {
        stop("Missing coordinates ", paste(allowed_coords, collapse = "/"),
             " in data points for covariate '", nm, "'.")
      }
      req <- c(coord_p, nm)
      if (!all(req %in% names(dfc0))) {
        miss <- setdiff(req, names(dfc0))
        stop(sprintf("Cov '%s' (interp): missing columns in covariate dataset: %s",
                     nm, paste(miss, collapse = ", ")))
      }
      
      pts <- quad_p[, coord_p, drop = FALSE]
      dfc <- dfc0[, req, drop = FALSE]
      pts <- coerce_time_numeric(pts, coord_p)
      dfc <- coerce_time_numeric(dfc, coord_p)
      
      cov.pred <- idw_interp_fast(
        points = pts,
        covs   = dfc,
        p      = interp.p,
        d      = length(coord_p)
      )
      points.covs[[nm]] <- cov.pred
      
    } else if (choice == "exact") {
      f_exact <- covs[[nm]]
      if (!is.function(f_exact)) {
        stop("For type 'exact', covs[['", nm, "']] must be a function of a data.frame with x,y coordinates.")
      }
      vals <- f_exact(as.data.frame(quad_p))
      if (!is.numeric(vals)) stop("Exact covariate function for '", nm, "' must return numeric values.")
      if (length(vals) == 1L) vals <- rep(vals, nrow(points.covs))
      if (length(vals) != nrow(points.covs)) stop("Exact covariate function for '", nm, "' returned the wrong length.")
      points.covs[[nm]] <- as.numeric(vals)

    } else {
      stop("No managed type for covariate '", nm, "': ", choice)
    }
  }
  
  points.covs
}



.select_coords <- function(df, d = NULL, coord_names = c("x","y","z","t"), require_all = FALSE) {
  if (!is.null(d)) {
    if (ncol(df) < d) stop("Argument 'd' exceeds the number of available columns.")
    return(as.matrix(df[, seq_len(d), drop = FALSE]))
  }
  if (!is.null(colnames(df))) {
    if (require_all) {
      if (!all(coord_names %in% colnames(df))) stop("Missing required coordinates: ", paste(coord_names, collapse = ", "))
      return(as.matrix(df[, coord_names, drop = FALSE]))
    } else {
      cn <- intersect(coord_names, colnames(df))
      if (length(cn) == 0L) stop("Unable to infer coordinate columns. Provide 'd' or name columns among {x,y,z,t}.")
      return(as.matrix(df[, cn, drop = FALSE]))
    }
  }
  stop("Provide 'd' or column names to infer coordinates.")
}




distN_fast <- function(A, B, d = NULL, coord_names = c("x","y","z","t"), require_all = FALSE) {
  if (is.null(d)) {
    A <- .select_coords(A, coord_names = coord_names, require_all = require_all)
    B <- .select_coords(B, coord_names = colnames(A), require_all = TRUE)
  } else {
    A <- .select_coords(A, d = d, coord_names = coord_names, require_all = require_all)
    B <- .select_coords(B, d = d, coord_names = coord_names, require_all = require_all)
  }
  if (!is.null(colnames(A)) && "t" %in% colnames(A) && !is.numeric(A[, "t"])) A[, "t"] <- as.numeric(A[, "t"])
  if (!is.null(colnames(B)) && "t" %in% colnames(B) && !is.numeric(B[, "t"])) B[, "t"] <- as.numeric(B[, "t"])
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  A2 <- rowSums(A^2)
  B2 <- rowSums(B^2)
  D2 <- outer(A2, B2, "+") - 2 * (A %*% t(B))
  D2[D2 < 0] <- 0
  sqrt(D2)
}



idw_interp_fast <- function(points, covs, p = 2, d = NULL,
                            coord_names = c("x","y","z","t"),
                            eps = 1e-10, exact_on_zero = TRUE, na.rm = TRUE,
                            require_all = FALSE, rescale01 = FALSE) {
  if (ncol(covs) < 2) stop("'covs' must have at least coordinates plus one value column.")
  value <- covs[[ncol(covs)]]
  if (!is.numeric(value)) stop("The last column of 'covs' must be numeric (values to interpolate).")
  if (na.rm) {
    keep <- !is.na(value)
    covs  <- covs[keep, , drop = FALSE]
    value <- value[keep]
  }
  if (nrow(covs) == 0) stop("No observations available in 'covs' after removing NAs.")
  
  if (is.null(d)) {
    P <- .select_coords(points, coord_names = coord_names, require_all = require_all)
    C <- .select_coords(covs,   coord_names = colnames(P), require_all = TRUE)
    d_use <- ncol(P)
  } else {
    P <- .select_coords(points, d = d, coord_names = coord_names, require_all = require_all)
    C <- .select_coords(covs,   d = d, coord_names = coord_names, require_all = require_all)
    d_use <- d
  }
  
  if (!is.null(colnames(P)) && "t" %in% colnames(P) && !is.numeric(P[, "t"])) P[, "t"] <- as.numeric(P[, "t"])
  if (!is.null(colnames(C)) && "t" %in% colnames(C) && !is.numeric(C[, "t"])) C[, "t"] <- as.numeric(C[, "t"])
  
  if (rescale01) {
    for (j in seq_len(ncol(P))) {
      combo <- c(P[, j], C[, j])
      combo_sc <- scale_to_range(combo, 0, 1)
      nP <- nrow(P)
      P[, j] <- combo_sc[seq_len(nP)]
      C[, j] <- combo_sc[(nP + 1L):length(combo_sc)]
    }
  }
  
  storage.mode(P) <- "double"
  storage.mode(C) <- "double"
  
  D <- distN_fast(P, C, d = d_use, coord_names = coord_names, require_all = FALSE)
  
  res <- numeric(nrow(P))
  if (exact_on_zero) {
    zero_mask <- (D == 0)
    rows_zero <- which(rowSums(zero_mask) > 0)
    if (length(rows_zero)) {
      for (i in rows_zero) res[i] <- mean(value[zero_mask[i, ]], na.rm = TRUE)
    }
    rows_rest <- setdiff(seq_len(nrow(P)), rows_zero)
  } else {
    rows_rest <- seq_len(nrow(P))
  }
  if (length(rows_rest)) {
    D_sub <- D[rows_rest, , drop = FALSE]
    D_sub[D_sub == 0] <- eps
    W <- 1 / (D_sub^p)
    denom <- rowSums(W)
    numer <- W %*% matrix(value, ncol = 1)
    res[rows_rest] <- as.vector(numer / denom)
  }
  res
}



build_poisson_formula <- function(formula,
                                  add_offset = FALSE,
                                  offset_var = "str_Klocal",
                                  data_cols = NULL) {
  # salva env originale (importantissimo per formule con costanti tipo xc_bg)
  f_env <- if (inherits(formula, "formula")) environment(formula) else parent.frame()
  if (is.null(f_env)) f_env <- parent.frame()
  
  # normalizza la formula in y_resp ~ RHS
  if (inherits(formula, "formula")) {
    if (length(formula) %in% c(2L, 3L)) {
      form_full <- update(formula, y_resp ~ .)
    } else {
      stop("Invalid 'formula': unexpected structure.")
    }
  } else if (is.character(formula)) {
    rhs_txt <- trimws(sub("^~", "", formula[1]))
    form_full <- as.formula(paste0("y_resp ~ ", rhs_txt), env = f_env)
  } else {
    stop("'formula' must be a formula or a character string.")
  }
  
  # ripristina esplicitamente l'environment
  environment(form_full) <- f_env
  
  # offset opzionale
  if (isTRUE(add_offset)) {
    form_full <- update(form_full,
                        as.formula(paste0(". ~ . + offset(", offset_var, ")"), env = f_env))
    environment(form_full) <- f_env
  }
  
  # validazione opzionale
  if (!is.null(data_cols)) {
    vars <- setdiff(all.vars(form_full), "y_resp")
    miss <- setdiff(vars, data_cols)
    if (length(miss)) {
      stop("Variables not found in data: ", paste(miss, collapse = ", "))
    }
  }
  
  form_full
}





counting.weights <- function(id, volumes) {
  id <- as.integer(id)
  fid <- factor(id, levels = seq_along(volumes))
  counts <- table(fid)
  w <- volumes[id] / counts[id]
  w <- as.vector(w)
  names(w) <- NULL
  return(w)
}

grid1.index <- function(x, xrange, nx) {
  i <- ceiling(nx * (x - xrange[1]) / diff(xrange))
  i <- pmax.int(1, i)
  i <- pmin.int(i, nx)
  i
}

rstpp.2D <- function (lambda = 500, nsim = 1, verbose = FALSE, par = NULL,
                      minX = 0, maxX = 1, minY = 0, maxY = 1) {
  if (is.numeric(lambda)) {
    par <- log(lambda)
    lambda <- function(x, y, a) exp(a[1])
  }
  if (nsim != 1) pp0 <- list(l = nsim)
  for (i in 1:nsim) {
    if (isTRUE(verbose)) progressreport(i, nsim)
    lam  <- lambda(1, 1, par)
    candn <- rpois(1, lam)
    candx <- runif(candn, minX, maxX)
    candy <- runif(candn, minY, maxY)
    d     <- runif(candn)
    lam2  <- lambda(candx, candy, par)
    lmax  <- max(lam2)
    keep  <- (d < lam2/lmax)
    lon   <- candx[keep]; lat <- candy[keep]
    if (nsim != 1) {
      pp0[[i]] <- spatstat.geom::ppp(x = lon, y = lat,
                                     window = spatstat.geom::owin(xrange = range(lon), yrange = range(lat)))
    } else {
      pp0 <- spatstat.geom::ppp(x = lon, y = lat,
                                window = spatstat.geom::owin(xrange = range(lon), yrange = range(lat)))
    }
  }
  return(pp0)
}

default.ncube.2D <- function(X){
  guess.ngrid <- floor((splancs::npts(X) / 2) ^ (1 / 3))
  max(5, guess.ngrid)
}

grid.index.2D <- function(x, y, xrange, yrange, nx, ny) {
  ix <- grid1.index(x, xrange, nx)
  iy <- grid1.index(y, yrange, ny)
  list(ix = ix, iy = iy, index = as.integer((iy - 1) * nx + ix))
}




build_covs_bg <- function(cat, 
                          covar_names,
                          process.type.bg = "s2d") {
  # cat: data.frame già passato per cat.select, con xcat.work, ycat.work, time.work, ecc.
  # covar_names: vettore di nomi delle covariate da usare (presenti in cat)
  # process.type.bg: "s2d", "st" o "s3d" (per ora usi "s2d")
  
  # controlla che le covariate esistano
  missing <- setdiff(covar_names, names(cat))
  if (length(missing) > 0) {
    stop("Le seguenti covariate non sono presenti in 'cat': ",
         paste(missing, collapse = ", "))
  }
  
  # seleziona le coordinate giuste in base al tipo di processo
  if (process.type.bg == "s2d") {
    coord_cols <- c("xcat.work", "ycat.work")
    coord_map  <- c(xcat.work = "x", ycat.work = "y")
  } else if (process.type.bg == "st") {
    coord_cols <- c("xcat.work", "ycat.work", "time.work")
    coord_map  <- c(xcat.work = "x", ycat.work = "y", time.work = "t")
  } else if (process.type.bg == "s3d") {
    # adatta questi nomi se nel tuo cat sono diversi
    coord_cols <- c("xcat.work", "ycat.work", "zcat.work")
    coord_map  <- c(xcat.work = "x", ycat.work = "y", zcat.work = "z")
  } else {
    stop("process.type.bg non riconosciuto: ", process.type.bg)
  }
  
  # controlla che le coordinate esistano
  missing_coord <- setdiff(coord_cols, names(cat))
  if (length(missing_coord) > 0) {
    stop("Le seguenti colonne di coordinate mancano in 'cat': ",
         paste(missing_coord, collapse = ", "))
  }
  
  coords <- cat[, coord_cols, drop = FALSE]
  # rinomina le colonne delle coordinate in x,y,(t/z)
  names(coords) <- unname(coord_map[names(coords)])
  
  # costruisci la lista di covariate
  covs_list <- lapply(covar_names, function(v) {
    df <- cbind(coords, cat[[v]])
    names(df)[ncol(df)] <- v
    df
  })
  names(covs_list) <- covar_names
  
  covs_list
}



##############################################
#### FUNZIONI AUSILIARIE ETAS PARAMETRICO ####
##############################################

etas.starting <- function (cat.orig, magn.threshold = 2.5, p.start = 1, gamma.start = 0.5, 
                           q.start = 2, betacov.start = 0.7, longlat.to.km = TRUE, sectoday = FALSE, 
                           onlytime = FALSE) 
{
  cat = cat.orig[cat.orig$magn1 > magn.threshold, ]
  if (sectoday) 
    cat$time = cat$time/86400
  t = cat$time
  a = cat[order(t), ]
  t = a$time
  if (longlat.to.km) {
    radius = 6371.3
    ycat.km = radius * a$lat * pi/180
    xcat.km = radius * a$long * pi/180
  }
  else {
    ycat.km = a$lat
    xcat.km = a$long
  }
  n = nrow(a)
  dt = diff(t)
  ds = diff(ycat.km)^2 + diff(xcat.km)^2
  c.start = as.numeric(quantile(dt, 0.25))
  d.start = as.numeric(quantile(ds, 0.05))
  mu.start = n * 0.5/diff(range(t))
  if (onlytime) {
    gamma.start = 0
    q.start = 0
    d.start = 0
  }
  tmax = max(t)
  it = log((c.start + tmax - t)/c.start)
  em = exp((betacov.start[1]) * (a$magn1 - magn.threshold))
  is = d.start^(1 - q.start) * exp(gamma.start * (a$magn1 - 
                                                    magn.threshold)) * pi/(q.start - 1)
  if (onlytime) 
    k0.start = n * 0.5/sum(em * it)
  else k0.start = n * 0.5/sum(em * it * is)
  return(list(mu.start = mu.start, k0.start = k0.start, c.start = c.start, 
              p.start = p.start, gamma.start = gamma.start, d.start = d.start, 
              q.start = q.start, betacov.start = betacov.start, longlat.to.km = longlat.to.km, 
              sectoday = sectoday))
}


embedding.rect.cat.epsNEW <- function (cat, cycle = FALSE, eps = 1/nrow(cat)) 
  embedding.rect.eps(cat$xcat.work, cat$ycat.work, cycle = cycle, eps = eps)


embedding.rect.eps <- function (x, y, cycle = FALSE, eps = 0) 
{
  eps.x = diff(range(x)) * eps/2
  eps.y = diff(range(y)) * eps/2
  rx = c(min(x) - eps.x, max(x) + eps.x)
  ry = c(min(y) - eps.y, max(y) + eps.y)
  rect = c(rx[2], rx[1], rx[1], rx[2], rx[2], ry[2], ry[2], 
           ry[1], ry[1], ry[2])
  rect = matrix(rect, 5, 2)
  if (!cycle) 
    rect = rect[1:4, ]
  return(rect)
}

eqcat <- function (x) 
{
  a <- c("time", "lat", "long", "z", "magn1")
  n <- 5 - length(intersect(a, names(x)))
  if (n > 0) {
    print("WRONG EARTHQUAKE CATALOG DEFINITION")
    print(c(n, " wrong variable names"))
    print("an object of class eqcat (earthquake catalog) must contains at least the five names")
    print(a)
    return(list(cat = x, ok = FALSE))
  }
  else {
    class(x) <- c("eqcat", "data.frame")
    return(list(cat = x, ok = TRUE))
  }
}


cat.select <- function (etas.obj, longlat.to.km, sectoday) 
{
  time.work = etas.obj$cat.orig$time
  if (sectoday) 
    time.work = time.work/86400
  unit = 1
  if (longlat.to.km) 
    unit = 6371.3 * pi/180
  ycat.work = etas.obj$cat.orig$lat * unit
  xcat.work = etas.obj$cat.orig$long * unit
  magnitude = etas.obj$cat.orig$magn1 - etas.obj$magn.threshold
  ind.magn1 = (etas.obj$cat.orig$magn1 >= etas.obj$magn.threshold)
  if (is.na(etas.obj$tmax)) 
    etas.obj$tmax = max(etas.obj$cat.orig$time[ind.magn1])
  if (is.na(etas.obj$long.range[1])) 
    etas.obj$long.range = range(etas.obj$cat.orig$long[ind.magn1])
  if (is.na(etas.obj$lat.range[1])) 
    etas.obj$lat.range = range(etas.obj$cat.orig$lat[ind.magn1])
  ind.time = (etas.obj$cat.orig$time <= etas.obj$tmax)
  ind.long = (etas.obj$cat.orig$long <= etas.obj$long.range[2]) & 
    (etas.obj$cat.orig$long >= etas.obj$long.range[1])
  ind.lat = (etas.obj$cat.orig$lat <= etas.obj$lat.range[2]) & 
    (etas.obj$cat.orig$lat >= etas.obj$lat.range[1])
  ind = as.logical(ind.magn1 * ind.time * ind.long * ind.lat)
  ord <- order(time.work)
  cat = data.frame(etas.obj$cat.orig, time.work, xcat.work, 
                   ycat.work, magnitude, ind, ord)
  etas.obj$cat = cat
  return(etas.obj)
}


bwd.nrd <- function (x, w = replicate(length(x), 1), d = 2) 
{
  if (length(x) < 2L) 
    stop("need at least 2 data points")
  m <- weighted.mean(x, w)
  return(sqrt(weighted.mean((x - m)^2, w)) * (length(x) * (d + 
                                                             2)/4)^(-1/(d + 4)))
}


region.P <- function (region, P, k = 100) 
{
  vert0 = cartesian2polar2D(region, P)
  vert0.theta = vert0[, 2]
  ind.theta = sort.list(vert0.theta)
  vert0 = vert0[ind.theta, ]
  n.vert0 = nrow(vert0)
  vert = matrix(0, length(vert0.theta) + 1, 2)
  vert[1:n.vert0, ] = vert0
  vert[n.vert0 + 1, ] = vert0[1, ]
  vert.theta = vert[, 2]
  vert.rho = vert[, 1]
  n.vert = nrow(vert)
  vert.seq = 1:(n.vert - 1)
  vert.gamma = angle.positive(diff(vert.theta))
  vert.alpha = triangle.abgamma.alpha(vert.rho[vert.seq + 1], 
                                      vert.rho[vert.seq], vert.gamma)
  theta = angle.positive(seq(0, 2 * pi, length.out = k + 1) + 
                           vert.theta[1] - 2 * pi)
  theta = theta[2:(k + 1)]
  ind1 = findInterval(theta, vert.theta[vert.seq])
  ind1 = ifelse(ind1 == 0, max(vert.seq), ind1)
  alpha.inf = vert.alpha[ind1]
  theta.inf = angle.positive(theta - vert.theta[ind1])
  rho.inf = vert.rho[ind1]
  rho = triangle.alphabetac.b(theta.inf, alpha.inf, rho.inf)
  return(list(vert.theta = vert.theta, vert.rho = vert.rho, 
              vert.gamma = vert.gamma, vert.alpha = vert.alpha, theta = theta, 
              vert.seq = vert.seq, theta.inf = theta.inf, rho.inf = rho.inf, 
              alpha.inf = alpha.inf, rho = rho, ind1 = ind1))
}


cartesian2polar2D <- function (x, orig = c(0, 0)) 
{
  x = x - outer(rep(1, dim(x)[1]), orig)
  y = x
  y[, 2] = angle.positive(atan2(x[, 2], x[, 1]))
  y[, 1] = sqrt(x[, 2] * x[, 2] + x[, 1] * x[, 1])
  return(y)
}

angle.positive <- function (theta) 
  ifelse(theta < 0, theta + 2 * pi, theta)


triangle.abgamma.alpha <- function (a, b, gamma) 
{
  c = triangle.abgamma.c(a, b, gamma)
  triangle.abc.gamma(b, c, a)
}

triangle.abgamma.c <- function (a, b, gamma) 
  sqrt(a * a + b * b - 2 * a * b * cos(gamma))

triangle.abc.gamma <- function (a, b, c) 
  acos((a * a + b * b - c * c)/(2 * a * b))

triangle.alphabetac.b <- function (alpha, beta, c) 
  c * sin(beta)/sin(alpha + beta)


generaloptimizationNEW <- function (etas.obj, hessian, iterlim, iprint, trace) 
{
  n = nrow(etas.obj$cat)
  params = c(etas.obj$params, etas.obj$betacov)
  etas.obj$nstep.par = etas.obj$nstep.par + 1
  if (iprint) {
    print("params in general optimization")
    print(params)
    print(etas.obj$nparams)
    print(etas.obj$nparams.etas)
  }
  etas.obj$dettrasf = replicate(n, 1)
  if (etas.obj$usenlm) {
    risult = nlm(etas.mod2NEW, params, typsize = abs(params), 
                 hessian = hessian, iterlim = iterlim, params.lim = etas.obj$params.lim, 
                 etas.obj = etas.obj)
    params.optim = risult$estimate
    l.optim = risult$minimum
  }
  else {
    risult = optim(params, etas.mod2NEW, method = etas.obj$method, 
                   params.lim = etas.obj$params.lim, hessian = hessian, 
                   control = list(trace = 2, maxit = iterlim, fnscale = n/diff(range(etas.obj$cat$time.work)), 
                                  parscale = sqrt(exp(params))), etas.obj = etas.obj)
    params.optim = risult$par
    l.optim = risult$value
  }
  return(list(params.optim = params.optim, l.optim = l.optim, 
              risult = risult, etas.obj = etas.obj))
}


etas.mod2NEW <- function (params = c(1, 1, 1, 1, 1, 1, 1), etas.obj, params.lim = c(0, 0, 0, 1, 0, 0, 1), trace = TRUE, iprint = FALSE) 
{
  params.ind = etas.obj$params.ind
  params.lim = etas.obj$params.lim
  params.fix = etas.obj$params.fix
  tmax = etas.obj$tmax
  cat = etas.obj$cat
  magn.threshold = etas.obj$magn.threshold
  back.dens = etas.obj$back.dens
  back.integral = etas.obj$back.integral
  onlytime = etas.obj$onlytime
  rho.s2 = etas.obj$rho.s2
  nparams.etas = etas.obj$nparams.etas
  nparams = etas.obj$nparams
  ntheta = etas.obj$ntheta
  params.etas = params[1:nparams.etas]
  betacov = params[(nparams.etas + 1):nparams]
  params.e = params.fix
  params.e[params.ind == 1] = exp(params.etas) + params.lim[params.ind == 
                                                              1]
  lambda = params.e[1]
  k0 = params.e[2]
  c = params.e[3]
  p = params.e[4]
  gamma = params.e[5]
  d = params.e[6]
  q = params.e[7]
  x = cat$xcat.work
  y = cat$ycat.work
  magnitudes = cat$magnitude
  predictor = as.matrix(etas.obj$cov.matrix) %*% as.vector(betacov) + 
    as.vector(etas.obj$offset)
  times = cat$time.work
  range.t = diff(range(times))
  time.init <- Sys.time()
  n = length(times)
  etas.comp = as.double(array(0, n))
  if (etas.obj$fastML) {
    ris = .Fortran("etasfull8fast", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), ind = as.integer(etas.obj$ind), 
                   nindex = as.integer(length(etas.obj$index.tot)), 
                   index = as.integer(etas.obj$index.tot), l = etas.comp)
  }
  else if (etas.obj$parallel) {
    ris = .Fortran("etasfull8newparallel", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), l = etas.comp)
  }
  else {
    ris = .Fortran("etasfull8newserial", NAOK = TRUE, tflag = as.integer(onlytime), 
                   n = as.integer(n), mu = as.double(lambda), k = as.double(k0), 
                   c = as.double(c), p = as.double(p), g = as.double(gamma), 
                   d = as.double(d), q = as.double(q), x = as.double(x), 
                   y = as.double(y), t = as.double(times), m = as.double(magnitudes), 
                   predictor = as.double(predictor), l = etas.comp)
  }
  timenow <- Sys.time()
  timeelapsed <- difftime(timenow, time.init, units = "secs")
  etas.comp = ris$l
  if (iprint) 
    cat(timeelapsed, "\n")
  ci = lambda * back.dens + etas.comp
  logL.l = sum(log(ci))
  {
    if (p == 1) {
      it = log(c + tmax - times) - log(c)
    }
    else {
      it = ((c + tmax - times)^(1 - p) - c^(1 - p))/(1 - 
                                                       p)
    }
    if (onlytime) {
      integralNEW = k0 * sum(exp(predictor) * it)
    }
    else {
      m1NEW = as.vector(exp(predictor))
      m2 = as.vector(exp(gamma * magnitudes))
      time.init <- Sys.time()
      ci = lambda * back.dens + etas.comp
      etasNEW = rowSums(m1NEW * m2 * ((rho.s2 * rho.s2 * 
                                         etas.obj$dettrasf/m2 + d)^(1 - q) - d^(1 - q)))
      spaceNEW = (pi/((1 - q) * ntheta)) * etasNEW
      integralNEW = sum(k0 * it * spaceNEW)
      if (iprint) {
        cat("integral polar", "\n")
        cat(integral, "\n")
      }
    }
    integral = integralNEW
    timenow <- Sys.time()
    timeelapsed <- difftime(timenow, time.init, units = "secs")
    integraltot = integral + lambda * range.t * back.integral
  }
  logL = -logL.l + integraltot
  if (iprint) {
    cat("whole Integral computation", "\n")
    cat(timeelapsed, "\n")
    cat(params, "\n")
    cat(exp(params), "\n")
  }
  attr(logL, "etas.vec") = etas.comp
  attr(logL, "lambda.vec") = ci
  attr(logL, "integraltot") = integraltot
  if (iprint) {
    cat("ML step: likelihoods and integral", (c(logL, logL.l, 
                                                integraltot)), "\n")
  }
  if (trace) 
    cat("*")
  return(logL)
}



###################################
#### FUNZIONE ETAS PARAMETRICO ####
###################################

etasclass.par <- function (
    cat.orig,
    time.update = FALSE,
    magn.threshold = 2.5,
    magn.threshold.back = magn.threshold + 2,
    tmax = max(cat.orig$time),
    long.range = range(cat.orig$long),
    lat.range = range(cat.orig$lat),
    mu = 1,
    k0 = 1,
    c = 0.5,
    p = 1.01,
    gamma = 0.5,
    d = 1,
    q = 1.5,
    betacov = 0.7,
    params.ind = replicate(7, TRUE),
    formula1 = "time~magnitude-1",
    offset = 0,
    hdef = c(1, 1),
    wp = replicate(nrow(cat.orig), 1),
    hvarx = replicate(nrow(cat.orig), 1),
    hvary = replicate(nrow(cat.orig), 1),
    declustering = TRUE,
    thinning = FALSE,
    flp = TRUE,
    m1 = NULL,
    ndeclust = 5,
    n.iterweight = 1,
    onlytime = FALSE,
    is.backconstant = FALSE,
    description = "",
    cat.back = NULL,
    back.smooth = 1,
    sectoday = FALSE,
    longlat.to.km = TRUE,
    usenlm = TRUE,
    method = "BFGS",
    compsqm = TRUE,
    epsmax = 1e-04,
    iterlim = 50,
    ntheta = 36,
    formula.bg = ~ s(x, k = 30) + s(y, k = 30),
    process.type.bg = "s2d",
    verbose.bg = TRUE,
    marked.bg = FALSE,
    mark.c.bg = FALSE,
    type.cov.values.bg = NULL,
    grid.bg = FALSE,
    mult.bg = 4,
    ncube.bg = NULL,
    spatial.cov.bg = FALSE,
    offset_k.bg = FALSE,
    seed.bg = NULL,
    covs.bg.user = NULL,
    min.declust.iter = 2,
    tol.rho = 1e-3,
    tol.par = 1e-3,
    stop.on.increasing.AIC = FALSE,
    verbose.convergence = TRUE
) {
  
  iprint <- FALSE
  fastML <- FALSE
  parallel <- FALSE
  fast.eps <- 0.001
  
  params.lim <- c(0, 0, 0, 1, 0, 0, 1)
  namespar <- c("mu", "k0", "c", "p", "gamma", "d", "q")
  
  this.call <- match.call()
  
  flag <- eqcat(cat.orig)
  if (!flag$ok) {
    stop("WRONG EARTHQUAKE CATALOG DEFINITION")
  }
  cat.orig <- flag$cat
  
  iter <- 0
  AIC.iter <- numeric(0)
  AIC.iter2 <- numeric(0)
  AIC.iter.raw <- numeric(0)
  AIC.iter.postbg <- numeric(0)
  
  params.iter <- numeric(0)
  sqm.iter <- numeric(0)
  rho.weights.iter <- numeric(0)
  hdef.iter <- numeric(0)
  wmat <- numeric(0)
  fl <- 0
  fl.iter <- numeric(0)
  
  AIC.decrease <- TRUE
  
  eps <- Inf
  eps.par <- Inf
  eps.rho <- Inf
  eps.bg <- Inf
  
  trace <- TRUE
  
  if (onlytime) {
    is.backconstant <- TRUE
    declustering <- FALSE
    params.ind[5:7] <- c(FALSE, FALSE, FALSE)
    gamma <- 0
    d <- 0
    q <- 0
  }
  
  if (is.backconstant) {
    declustering <- FALSE
  }
  
  if (!declustering) {
    ndeclust <- 1
    thinning <- FALSE
    flp <- FALSE
  }
  
  if (missing(tmax)) {
    is.na(tmax) <- TRUE
  }
  if (missing(long.range)) {
    is.na(long.range) <- TRUE
  }
  if (missing(lat.range)) {
    is.na(lat.range) <- TRUE
  }
  
  etas.current <- list(
    parallel = parallel,
    this.call = match.call(),
    nstep.flp = 0,
    nstep.ppp = 0,
    nstep.par = 0,
    description = description,
    time.start = Sys.time(),
    magn.threshold = magn.threshold,
    magn.threshold.back = magn.threshold.back,
    onlytime = onlytime,
    tmax = tmax,
    lat.range = lat.range,
    long.range = long.range,
    hvarx = hvarx,
    hvary = hvary,
    is.backconstant = is.backconstant,
    usenlm = usenlm,
    method = method,
    cat.orig = cat.orig,
    declustering = declustering,
    thinning = thinning,
    flp = flp,
    back.smooth = back.smooth,
    ndeclust = ndeclust,
    eps = eps,
    longlat.to.km = longlat.to.km,
    sectoday = sectoday,
    ntheta = ntheta,
    model.bg = NULL
  )
  
  class(etas.current) <- "etasclass"
  
  etas.current <- cat.select(
    etas.current,
    longlat.to.km,
    sectoday = sectoday
  )
  
  etas.current$cat <- data.frame(
    etas.current$cat,
    hvarx,
    hvary
  )
  
  cat <- etas.current$cat[etas.current$cat$ord, ][
    etas.current$cat$ind[etas.current$cat$ord],
  ]
  
  etas.current$cat <- cat
  n <- nrow(cat)
  
  missingwp <- missing(wp) || is.null(wp)
  
  if (missingwp) {
    wp_event <- rep(1, n)
  } else {
    wp_event <- as.numeric(wp)
    
    if (time.update) {
      n.old <- length(wp_event)
      if (n.old < n) {
        wp_event <- c(wp_event, rep(0.5, n - n.old))
      }
    }
    
    if (length(wp_event) != n) {
      stop(
        "Length of 'wp' must match the selected catalogue size. ",
        "length(wp) = ", length(wp_event),
        ", nrow(cat) = ", n
      )
    }
  }
  
  if (time.update) {
    n.old <- length(wp_event)
    hvarx <- c(hvarx, rep(1, n - n.old))
    hvary <- c(hvary, rep(1, n - n.old))
  } else {
    hvarx <- cat$hvarx / (prod(cat$hvarx))^(1 / n)
    hvary <- cat$hvary / (prod(cat$hvary))^(1 / n)
  }
  
  etas.current$w <- wp_event
  
  ycat.work <- cat$ycat.work
  xcat.work <- cat$xcat.work
  
  if (missing(cat.back) || is.null(cat.back)) {
    ind.back <- (cat$magn1 >= magn.threshold.back)
  } else {
    stop("cat.back argument no more allowed, computed internally")
  }
  
  xback.work <- cat$xcat.work[ind.back]
  yback.work <- cat$ycat.work[ind.back]
  
  if (!missingwp) {
    xback.work <- xcat.work
    yback.work <- ycat.work
  }
  
  if (missingwp) {
    wp_hdef <- rep(1, length(xback.work))
  } else {
    wp_hdef <- wp_event
  }
  
  starting <- etas.starting(
    cat,
    magn.threshold = magn.threshold,
    longlat.to.km = longlat.to.km,
    sectoday = sectoday,
    p.start = p,
    gamma.start = gamma,
    q.start = q,
    betacov.start = betacov[1],
    onlytime = onlytime
  )
  
  if (missing(m1) || is.null(m1)) {
    m1 <- as.integer(nrow(cat) / 2)
  }
  if (missing(mu) || is.null(mu)) {
    mu <- starting$mu.start
  }
  if (missing(k0) || is.null(k0)) {
    k0 <- starting$k0.start
  }
  if (missing(c) || is.null(c)) {
    c <- starting$c.start
  }
  if (missing(d) || is.null(d)) {
    d <- starting$d.start
  }
  if (missing(hdef) || is.null(hdef)) {
    hdef <- c(
      bwd.nrd(xback.work, wp_hdef),
      bwd.nrd(yback.work, wp_hdef)
    )
  }
  
  print("Initial ETAS params estimates after checking: ")
  vecpar.init <- c(mu, k0, c, p, gamma, d, q)
  names(vecpar.init) <- namespar
  print(round(vecpar.init, 4))
  
  etas.current$xback.work <- xback.work
  etas.current$yback.work <- yback.work
  etas.current$mu.start <- mu
  etas.current$k0.start <- k0
  etas.current$c.start <- c
  etas.current$d.start <- d
  etas.current$p.start <- p
  etas.current$q.start <- q
  etas.current$gamma.start <- gamma
  etas.current$betacov.start <- betacov
  etas.current$hdef.start <- hdef
  etas.current$formula1 <- formula1
  
  formula1 <- as.formula(formula1)
  formula1 <- update(formula1, . ~ . - 1)
  
  cov.matrix <- model.matrix(formula1, data = cat)
  offset <- model.offset(model.frame(formula1, data = cat))
  
  if (is.null(offset)) {
    offset <- rep(0, nrow(cat))
  }
  
  if (!is.element("matrix", class(cov.matrix))) {
    stop("WRONG FORMULA DEFINITION")
  }
  
  ncov <- ncol(cov.matrix)
  
  if (length(betacov) != ncov) {
    cat(
      "Wrong number of starting values for parameters linear predictor of covariates. ",
      "Correct number of zero starting values inserted",
      "\n"
    )
    betacov <- rep(0, ncov)
  }
  
  if (length(params.ind) != 7) {
    stop("Wrong number of elements of params.ind")
  }
  
  params.ind <- as.numeric(params.ind)
  
  if (sum(abs(params.ind - 0.5) == 0.5) != 7) {
    stop("WRONG params.ind DEFINITION: ONLY FALSE/TRUE ALLOWED SEE HELP")
  }
  
  params.fix <- c(mu, k0, c, p, gamma, d, q)
  nparams.etas <- sum(params.ind)
  nparams <- nparams.etas + ncov
  
  if (onlytime) {
    
    if (!(prod(params.fix[1:4] >= params.lim[1:4]))) {
      stop(
        "WRONG starting values for parameters: ",
        "mu, k0, c must be positive, p must be equal or greater than 1"
      )
    }
    
    rho.s2 <- 0
    region <- NULL
    
  } else {
    
    if (!(prod(params.fix >= params.lim))) {
      stop(
        "WRONG starting values for parameters: ",
        "mu, k0, c, gamma, d, q must be positive, ",
        "p must be equal or greater than 1"
      )
    }
    
    region <- embedding.rect.cat.epsNEW(cat)
    rho.s2 <- matrix(0, ntheta, n)
    
    for (i in seq_len(n)) {
      rho.s2[, i] <- region.P(
        region,
        c(xcat.work[i], ycat.work[i]),
        k = ntheta
      )$rho
    }
    
    rho.s2 <- t(rho.s2)
  }
  
  etas.current$cov.matrix <- cov.matrix
  etas.current$region <- region
  etas.current$rho.s2 <- rho.s2
  etas.current$offset <- offset
  
  if (!inherits(formula.bg, "formula")) {
    formula.bg <- as.formula(formula.bg)
  }
  
  if (isTRUE(spatial.cov.bg)) {
    
    vars_bg <- all.vars(formula.bg)
    
    if (is.null(type.cov.values.bg) || is.null(names(type.cov.values.bg))) {
      stop("If spatial.cov.bg = TRUE, 'type.cov.values.bg' with covariate names is needed")
    }
    
    if (!is.null(covs.bg.user)) {
      
      if (!is.list(covs.bg.user) || is.null(names(covs.bg.user))) {
        stop("covs.bg.user must be a named list when supplied.")
      }
      
      covar_bg_names <- intersect(
        vars_bg,
        intersect(names(type.cov.values.bg), names(covs.bg.user))
      )
      
      if (!length(covar_bg_names)) {
        stop("No valid user-supplied background covariates found in formula/type.cov.values.bg/covs.bg.user.")
      }
      
      covs.bg <- covs.bg.user[covar_bg_names]
      type.cov.values.bg <- type.cov.values.bg[covar_bg_names]
      
    } else {
      
      covar_bg_names <- intersect(
        vars_bg,
        intersect(names(type.cov.values.bg), names(cat))
      )
      
      if (!length(covar_bg_names)) {
        stop("No valid covariates found")
      }
      
      covs.bg <- build_covs_bg(
        cat = cat,
        covar_names = covar_bg_names,
        process.type.bg = process.type.bg
      )
      
      type.cov.values.bg <- type.cov.values.bg[covar_bg_names]
    }
    
  } else {
    covs.bg <- NULL
  }
  
  bg_fit <- NULL
  X0 <- NULL
  back_dens_obs <- rep(1, n)
  back_integral_hat <- 1
  
  if (!onlytime) {
    
    .unit.bg <- if (isTRUE(longlat.to.km)) 6371.3 * pi / 180 else 1
    .xrange.bg <- etas.current$long.range * .unit.bg
    .yrange.bg <- etas.current$lat.range * .unit.bg
    
    win.bg <- spatstat.geom::owin(
      xrange = .xrange.bg,
      yrange = .yrange.bg
    )
    
    X0 <- make_bg_ppp_with_optional_marks(
      cat = cat,
      win.bg = win.bg,
      formula.bg = formula.bg
    )
  }
  
  fit_background_parametric <- function(pb_event, force_constant = FALSE) {
    
    if (onlytime) {
      return(list(
        bg_fit = NULL,
        back_dens_obs = rep(1, n),
        back_integral_hat = 1
      ))
    }
    
    if (length(pb_event) != n) {
      stop(
        "Internal error: background weights must have length n. ",
        "length(pb_event) = ", length(pb_event),
        ", n = ", n
      )
    }
    
    formula.bg.fit <- if (isTRUE(force_constant)) {
      ~ 1
    } else {
      formula.bg
    }
    
    marked.bg.fit <- if (isTRUE(force_constant)) FALSE else marked.bg
    mark.c.bg.fit <- if (isTRUE(force_constant)) FALSE else mark.c.bg
    spatial.cov.bg.fit <- if (isTRUE(force_constant)) FALSE else spatial.cov.bg
    covs.bg.fit <- if (isTRUE(force_constant)) NULL else covs.bg
    type.cov.values.bg.fit <- if (isTRUE(force_constant)) NULL else type.cov.values.bg
    
    bg_fit_local <- WPI_bkgd.fit(
      X = X0,
      formula = formula.bg.fit,
      process.type = process.type.bg,
      verbose = verbose.bg,
      mult = mult.bg,
      ncube = ncube.bg,
      spatial.cov = spatial.cov.bg.fit,
      marked = marked.bg.fit,
      mark.c = mark.c.bg.fit,
      covs = covs.bg.fit,
      pb = pb_event,
      type.cov.values = type.cov.values.bg.fit,
      grid = grid.bg,
      mark.c.mode = c("idw", "mc"),
      seed = seed.bg
    )
    
    if (verbose.bg) {
      print(summary(bg_fit_local$mod_global))
    }
    
    bg_combined <- if (!is.null(bg_fit_local$quad.data)) {
      bg_fit_local$quad.data
    } else {
      rbind_fill_base(bg_fit_local$newdata, bg_fit_local$dummy.data)
    }
    
    lambda_all <- exp(predict(
      bg_fit_local$mod_global,
      newdata = bg_combined
    ))
    
    w_model <- bg_combined$w
    kappa_hat <- sum(w_model * lambda_all)
    
    if (!is.finite(kappa_hat) || kappa_hat <= 0) {
      stop("Background normalization failed: non-positive or non-finite kappa_hat.")
    }
    
    back_dens_all <- lambda_all / kappa_hat
    back_integral_hat <- sum(w_model * back_dens_all)
    
    lambda_obs <- bg_fit_local$l
    back_dens_obs_local <- lambda_obs / kappa_hat
    back_dens_obs_local <- pmax(back_dens_obs_local, .Machine$double.eps)
    
    list(
      bg_fit = bg_fit_local,
      back_dens_obs = back_dens_obs_local,
      back_integral_hat = back_integral_hat
    )
  }
  
  bg0 <- fit_background_parametric(
    pb_event = wp_event,
    force_constant = isTRUE(is.backconstant)
  )
  
  bg_fit <- bg0$bg_fit
  back_dens_obs <- bg0$back_dens_obs
  back_integral_hat <- bg0$back_integral_hat
  
  etas.current$model.bg <- bg_fit
  etas.current$back.integral <- back_integral_hat
  etas.current$back.dens <- back_dens_obs
  back_dens_prev_for_conv <- NULL
  
  repeat {
    
    if (iter >= ndeclust) {
      if (verbose.convergence) {
        cat("\nStopping: reached maximum number of declustering iterations: ", ndeclust, "\n")
      }
      break
    }
    
    if (
      iter >= min.declust.iter &&
      is.finite(eps.rho) &&
      is.finite(eps.par) &&
      eps.rho <= tol.rho &&
      eps.par <= tol.par
    ) {
      if (verbose.convergence) {
        cat("\nStopping: EM/declustering convergence reached.\n")
        cat("eps.rho =", signif(eps.rho, 4), "\n")
        cat("eps.par =", signif(eps.par, 4), "\n")
      }
      break
    }
    
    prev_rho <- if (iter >= 1 && length(rho.weights.iter)) {
      as.numeric(rho.weights.iter[iter, ])
    } else {
      NULL
    }
    
    prev_params <- if (iter >= 1 && length(params.iter)) {
      as.numeric(params.iter[iter, ])
    } else {
      NULL
    }
    
    prev_back_dens <- back_dens_prev_for_conv
    
    params <- log((params.fix - params.lim)[params.ind == 1])
    
    etas.current$params <- params
    etas.current$params.lim <- params.lim
    etas.current$params.ind <- as.logical(params.ind)
    etas.current$params.fix <- params.fix
    etas.current$betacov <- betacov
    etas.current$nparams <- nparams
    etas.current$nparams.etas <- nparams.etas
    
    cat("Start ML step number: ", iter + 1, "\n")
    
    etas.current$fastML <- fastML
    etas.current$fast.eps <- fast.eps
    
    if (fastML) {
      fast <- fastML.init(etas.current)
      etas.current$ind <- fast$ind
      etas.current$index.tot <- fast$index.tot
    }
    
    risult.opt <- generaloptimizationNEW(
      etas.current,
      hessian = TRUE,
      iterlim = iterlim,
      iprint = iprint,
      trace = trace
    )
    
    etas.current <- risult.opt$etas.obj
    
    params.optim <- risult.opt$params.optim[seq_len(nparams.etas)]
    betacov <- risult.opt$params.optim[(nparams.etas + 1):nparams]
    
    for (i in seq_len(n.iterweight)) {
      
      cat("weighting step n. ", i, "\n")
      
      l.optim <- risult.opt$l.optim
      risult <- risult.opt$risult
      
      l <- etas.mod2NEW(
        params = c(params.optim, betacov),
        etas.obj = etas.current,
        trace = FALSE
      )
      
      ############################################################
      ## ROBUST STANDARD ERRORS FROM HESSIAN
      ############################################################
      
      sqm <- rep(NA_real_, nparams)
      
      if (isTRUE(compsqm)) {
        
        H <- risult$hessian
        
        if (
          is.matrix(H) &&
          all(dim(H) == c(nparams, nparams)) &&
          all(is.finite(H))
        ) {
          
          invH <- try(solve(H), silent = TRUE)
          
          if (!inherits(invH, "try-error")) {
            
            diag_invH <- diag(invH)
            
            if (all(is.finite(diag_invH)) && all(diag_invH >= 0)) {
              
              scale_vec <- c(exp(params.optim), rep(1, ncov))
              sqm <- sqrt(diag_invH) * scale_vec
              
            } else {
              warning(
                "Hessian inverse has negative or non-finite diagonal entries; ",
                "standard errors set to NA."
              )
            }
            
          } else {
            warning("Hessian could not be inverted; standard errors set to NA.")
          }
          
        } else {
          warning("Invalid Hessian; standard errors set to NA.")
        }
      }
      
      sqm.etas <- sqm[seq_len(nparams.etas)]
      sqm.cov <- sqm[(nparams.etas + 1):nparams]
      
      params <- params.fix
      params[params.ind == 1] <- exp(params.optim) + params.lim[params.ind == 1]
      
      sqm.tot <- rep(NA_real_, 7)
      names(sqm.tot) <- namespar
      
      sqm.tot[params.ind == 1] <- sqm.etas
      
      ## Parametri fissati: std.err. = 0 per indicare che non sono stimati
      sqm.tot[params.ind == 0] <- 0
      
      params.MLtot <- c(params, betacov)
      sqm.MLtot <- c(sqm.tot, sqm.cov)
      
      cat("found optimum; end ML step ")
      cat(iter + 1, "\n")
      
      mu <- params[1]
      k0 <- params[2]
      c <- params[3]
      p <- params[4]
      gamma <- params[5]
      d <- params[6]
      q <- params[7]
      
      predictor <- as.matrix(etas.current$cov.matrix) %*%
        as.vector(betacov) + as.vector(etas.current$offset)
      
      lambda_vec <- attr(l, "lambda.vec")
      
      if (any(!is.finite(lambda_vec)) || any(lambda_vec <= 0)) {
        stop("Invalid lambda.vec returned by etas.mod2NEW.")
      }
      
      rho.weights.raw <- params.MLtot[1] * back_dens_obs / lambda_vec
      
      if (any(rho.weights.raw < -1e-8 | rho.weights.raw > 1 + 1e-8, na.rm = TRUE)) {
        warning(
          "Some posterior background weights are outside [0,1]. ",
          "They will be truncated for the next background fit."
        )
      }
      
      rho.weights <- pmin(pmax(rho.weights.raw, 0), 1)
      wp_event <- rho.weights
      
      xback.work <- xcat.work
      yback.work <- ycat.work
      
      if (flp) {
        
        etas.l <- attr(l, "etas.vec")
        
        ris.flp <- flp1.etas.nlmNEW(
          cat,
          h.init = hdef,
          etas.params = params.MLtot,
          etas.l = etas.l,
          wp = rho.weights,
          m1 = m1,
          m2 = as.integer(nrow(cat) - 1),
          mh = 1
        )
        
        etas.current$nstep.flp <- etas.current$nstep.flp + 1
        hdef <- ris.flp$hdef
        fl <- ris.flp$fl
      }
      
      if (thinning) {
        back.ind <- runif(n) < rho.weights
        xback.work <- xcat.work[back.ind]
        yback.work <- ycat.work[back.ind]
      }
      
      if (flp) {
        
        back.tot <- kde2dnew.fortran(
          xback.work,
          yback.work,
          xcat.work,
          ycat.work,
          eps = 1 / n,
          h = hdef,
          w = rho.weights,
          hvarx = hvarx,
          hvary = hvary
        )
        
      } else {
        
        bg_iter <- fit_background_parametric(
          pb_event = wp_event,
          force_constant = isTRUE(is.backconstant)
        )
        
        bg_fit <- bg_iter$bg_fit
        back_dens_obs <- bg_iter$back_dens_obs
        back_integral_hat <- bg_iter$back_integral_hat
        
        etas.current$nstep.ppp <- etas.current$nstep.ppp + 1
      }
      
      etas.current$model.bg <- bg_fit
      etas.current$rho.weights <- rho.weights
      etas.current$back.integral <- back_integral_hat
      etas.current$back.dens <- back_dens_obs
      
      l <- etas.mod2NEW(
        params = c(params.optim, betacov),
        etas.obj = etas.current,
        trace = FALSE
      )
    }
    
    time.res <- array(0, n)
    
    if (onlytime) {
      
      times.tot <- cat$time
      magnitudes.tot <- cat$magn1 - magn.threshold
      tmin <- min(times.tot)
      etas.t <- 0
      
      for (i in seq_len(n)) {
        
        tmax.i <- times.tot[i]
        times <- subset(times.tot, times.tot < tmax.i)
        magnitudes <- subset(magnitudes.tot, times.tot < tmax.i)
        
        if (p == 1) {
          it <- log(c + tmax.i - times) - log(c)
        } else {
          it <- ((c + tmax.i - times)^(1 - p) - c^(1 - p)) / (1 - p)
        }
        
        time.res[i] <- (tmax.i - tmin) * mu +
          k0 * sum(exp(predictor) * it)
      }
    }
    
    names(params.fix) <- namespar
    names(params.ind) <- namespar
    names(betacov) <- colnames(cov.matrix)
    names(params.MLtot) <- c(namespar, names(betacov))
    names(sqm.MLtot) <- c(namespar, names(betacov))
    
    timenow <- Sys.time()
    time.elapsed <- difftime(
      timenow,
      etas.current$time.start,
      units = "secs"
    )
    
    edf.bg <- if (!is.null(bg_fit) && !is.null(bg_fit$mod_global$edf)) {
      sum(bg_fit$mod_global$edf)
    } else {
      0
    }
    
    nparams.aic.raw <- nparams
    nparams.aic.bg <- nparams + edf.bg
    
    AIC.ML.raw <- 2 * l.optim + 2 * nparams.aic.raw
    AIC.ML.bg <- 2 * l.optim + 2 * nparams.aic.bg
    
    AIC.postbg.raw <- 2 * l + 2 * nparams.aic.raw
    AIC.postbg.bg <- 2 * l + 2 * nparams.aic.bg
    
    AIC.temp <- AIC.ML.bg
    AIC.decrease <- ((iter == 0) || (AIC.temp <= min(AIC.iter)))
    
    if (!is.null(prev_rho)) {
      eps.rho <- mean(abs(rho.weights - prev_rho), na.rm = TRUE)
    } else {
      eps.rho <- Inf
    }
    
    if (!is.null(prev_params)) {
      eps.par <- max(
        abs(params.MLtot - prev_params) /
          (abs(prev_params) + 1e-8),
        na.rm = TRUE
      )
    } else {
      eps.par <- Inf
    }
    
    if (!is.null(prev_back_dens) && length(prev_back_dens) == length(back_dens_obs)) {
      
      eps.bg <- mean(
        abs(back_dens_obs - prev_back_dens) /
          (abs(prev_back_dens) + 1e-8),
        na.rm = TRUE
      )
      
    } else {
      
      eps.bg <- Inf
      
    }
    
    back_dens_prev_for_conv <- back_dens_obs
    
    eps <- eps.rho
    
    iter <- iter + 1
    
    AIC.iter[iter] <- AIC.ML.bg
    AIC.iter2[iter] <- AIC.postbg.bg
    AIC.iter.raw[iter] <- AIC.ML.raw
    AIC.iter.postbg[iter] <- AIC.postbg.raw
    
    params.iter <- rbind(params.iter, params.MLtot)
    sqm.iter <- rbind(sqm.iter, sqm.tot)
    rho.weights.iter <- rbind(rho.weights.iter, rho.weights)
    hdef.iter <- rbind(hdef.iter, hdef)
    fl.iter <- c(fl.iter, fl)
    
    cat(
      "\n",
      "---ITERATION n.", iter,
      " ---; AIC_ML_bg = ", round(AIC.iter[iter], 3),
      "; AIC_postbg_bg = ", round(AIC.iter2[iter], 3),
      "; elapsed time: ", time.elapsed,
      " time: "
    )
    print(Sys.time())
    
    cat("\n")
    cat("Current estimates of parameters: ", "\n")
    print(round(params.MLtot, 5))
    
    if (verbose.convergence) {
      cat(
        "Convergence diagnostics:",
        "eps.rho =", signif(eps.rho, 4),
        "; eps.par =", signif(eps.par, 4),
        "; eps.bg =", signif(eps.bg, 4),
        "; mean rho =", signif(mean(rho.weights), 4),
        "; sum rho =", signif(sum(rho.weights), 4),
        "\n"
      )
    }
    
    rownames(params.iter) <- seq_len(nrow(params.iter))
    
    predictor <- as.matrix(etas.current$cov.matrix) %*%
      as.vector(betacov) + as.vector(etas.current$offset)
    
    names(params.MLtot) <- c(namespar, names(betacov))
    names(sqm.MLtot) <- names(params.MLtot)
    
    etas.current$time.elapsed <- time.elapsed
    etas.current$time.end <- timenow
    etas.current$risult <- risult
    etas.current$betacov <- betacov
    etas.current$params.MLtot <- params.MLtot
    etas.current$params <- params.MLtot[1:7]
    etas.current$predictor <- predictor
    etas.current$sqm <- sqm.MLtot
    etas.current$time.res <- time.res
    etas.current$AIC.iter <- AIC.iter
    etas.current$AIC.iter2 <- AIC.iter2
    etas.current$AIC.iter.raw <- AIC.iter.raw
    etas.current$AIC.iter.postbg <- AIC.iter.postbg
    etas.current$AIC.decrease <- AIC.decrease
    etas.current$params.iter <- params.iter
    etas.current$sqm.iter <- sqm.iter
    etas.current$rho.weights.iter <- rho.weights.iter
    etas.current$rho.weights <- rho.weights
    etas.current$n.iterweight <- n.iterweight
    etas.current$hdef <- hdef
    etas.current$hdef.iter <- hdef.iter
    etas.current$back.integral <- back_integral_hat
    etas.current$back.dens <- back_dens_obs
    etas.current$model.bg <- bg_fit
    etas.current$wmat <- wmat
    etas.current$fl.iter <- fl.iter
    etas.current$iter <- iter
    etas.current$usenlm <- usenlm
    etas.current$method <- method
    etas.current$epsmax <- epsmax
    etas.current$tol.rho <- tol.rho
    etas.current$tol.par <- tol.par
    etas.current$eps.rho <- eps.rho
    etas.current$eps.par <- eps.par
    etas.current$eps.bg <- eps.bg
    etas.current$iterlim <- iterlim
    etas.current$compsqm <- compsqm
    etas.current$ntheta <- ntheta
    etas.current$logl <- l
    etas.current$l <- attr(l, "lambda.vec")
    etas.current$integral <- attr(l, "integraltot")
    
    etas.ris <- etas.current
    
    if (!AIC.decrease && isTRUE(stop.on.increasing.AIC)) {
      print("INCREASING AIC")
      return(etas.ris)
    }
  }
  
  if (!time.update) {
    cat("\n")
    cat("\n")
    cat("--------- FINAL SUMMARY ---------", "\n")
    summary(etas.ris)
  }
  
  return(etas.ris)
}



##############################################################
#### CATALOGO ITALIA PROVA FUNZIONAMENTO ETAS PARAMETRICO ####
##############################################################

data("italycatalog")
str(catalog.withcov)

#prova_etaspar <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
#                               tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
#                               mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
#                               params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
#                               offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
#                               hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
#                               ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1, 
#                               sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
#                               formula.bg = ~ s(x, y, bs = "tp", k = 50), process.type.bg = "s2d", spatial.cov.bg = FALSE, mult.bg = 4, ncube.bg = NULL,
#                               verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)




#########################################################
#### FUNZIONE PER SIMULARE PROCESSO ETAS PARAMETRICO ####
#########################################################

############################################################
# ETAS simulation v4
# - Background shape + triggering with true parameters
# - Supports spatial covariates in background (compatible with spatial.cov.bg=TRUE)
#
# Compatibility idea with etasclass.par:
#   If cat.sim contains columns Z1, Z2, ...
#   and formula.bg includes them, you can fit with:
#     spatial.cov.bg = TRUE
#     type.cov.values.bg = list(Z1="interp", Z2="interp", ...)
#
# Notes:
# - mu controls expected number of background events per unit time (Poisson on time window)
# - background LP controls only spatial SHAPE (normalized to density)
# - background intercept is not separately identifiable from mu
############################################################

etas.par.sim_v4 <- function(
    # ---- ETAS true parameters ----
    params = c(mu = 0.5, k0 = 0.012, c = 0.01, p = 1.12, gamma = 0, d = 1.6, q = 1.85),
    m0 = 2.5,
    b  = 1.05,
    
    # ---- Time window ----
    tmin  = 0,
    t.lag = 700,
    
    # ---- Spatial window (degrees if longlat.to.km=TRUE) ----
    long.range = c(6, 20),
    lat.range  = c(36, 45),
    longlat.to.km = TRUE,
    sectoday = FALSE,
    
    # ---- Background support (quadrature-like support for simulation) ----
    n_support_obs = 1500,
    mult_support  = 4,         # dummy support = mult_support * n_support_obs
    
    # ---- Background true LP (spatial shape) ----
    # "constant" : eta_bg = 0
    # "linear"   : eta_bg = sum_j beta_j * X_j  (X_j can be x,y,x_std,y_std,r2_center and/or spatial covariates)
    # "custom"   : eta_bg = bg_lp_fun(df_features)
    bg_lp_type  = c("constant", "linear", "custom"),
    bg_lp_coefs = NULL,        # named numeric vector for linear case
    bg_lp_fun   = NULL,        # function(df_features) -> numeric vector
    
    # ---- Spatial covariates for background (NEW) ----
    # Option A (recommended): named list of deterministic functions of df_features
    #   bg_spatial_cov_funs = list(
    #     Z1 = function(df) sqrt((df$x - mean(df$x))^2 + (df$y - mean(df$y))^2),
    #     Z2 = function(df) sin(pi * df$x_std) * cos(pi * df$y_std)
    #   )
    #
    # Option B: a generator function returning a data.frame with named columns
    #   bg_spatial_cov_generator = function(df_features, window_km, seed=NULL) data.frame(...)
    #
    # If both are provided, columns are merged (generator first, then funs can overwrite same names)
    bg_spatial_cov_funs      = NULL,
    bg_spatial_cov_generator = NULL,
    
    # ---- Triggering true LP ----
    # Usually gamma=0 and trig_lp_coefs = c(m_rel = betacov_true)
    trig_lp_type  = c("linear", "custom"),
    trig_lp_coefs = c(m_rel = 0.7),   # named vector on event-level features
    trig_lp_fun   = NULL,             # function(df_event_features) -> numeric scalar/vector
    
    # ---- Safety / numerics ----
    seed = NULL,
    clip_eta_bg = 20,
    clip_eta_trig = 20,
    max_events = 50000,
    explode_guard = TRUE,
    
    # ---- Returns / plots ----
    return_bg_info = TRUE,
    return_support = TRUE,
    plot_bg = FALSE,
    plot_catalog = FALSE,
    cex_cat = 0.5
) {
  # ---------------------------
  # 0) checks and parsing
  # ---------------------------
  bg_lp_type   <- match.arg(bg_lp_type)
  trig_lp_type <- match.arg(trig_lp_type)
  
  # ETAS params
  need_names <- c("mu","k0","c","p","gamma","d","q")
  if (is.null(names(params)) || !all(need_names %in% names(params))) {
    stop("params must be a named vector containing: ", paste(need_names, collapse = ", "))
  }
  mu    <- unname(params["mu"])
  k0    <- unname(params["k0"])
  c_    <- unname(params["c"])
  p     <- unname(params["p"])
  gamma <- unname(params["gamma"])
  d     <- unname(params["d"])
  q     <- unname(params["q"])
  
  if (!is.finite(mu) || mu <= 0) stop("mu must be > 0")
  if (!is.finite(k0) || k0 <= 0) stop("k0 must be > 0")
  if (!is.finite(c_) || c_ <= 0) stop("c must be > 0")
  if (!is.finite(p)  || p <= 1) stop("p must be > 1")
  if (!is.finite(d)  || d <= 0) stop("d must be > 0")
  if (!is.finite(q)  || q <= 1) stop("q must be > 1")
  if (!is.finite(b)  || b <= 0) stop("b must be > 0")
  if (!is.finite(m0)) stop("m0 must be finite")
  
  if (!is.null(seed)) set.seed(seed)
  
  tmax <- tmin + t.lag
  if (!is.finite(tmax) || tmax <= tmin) stop("Invalid time window: tmax <= tmin")
  
  # coordinate units
  unit <- if (isTRUE(longlat.to.km)) 6371.3 * pi / 180 else 1
  
  xmin_km <- long.range[1] * unit
  xmax_km <- long.range[2] * unit
  ymin_km <- lat.range[1]  * unit
  ymax_km <- lat.range[2]  * unit
  
  if (xmin_km >= xmax_km || ymin_km >= ymax_km) stop("Invalid long.range / lat.range")
  area_km2 <- (xmax_km - xmin_km) * (ymax_km - ymin_km)
  
  # ---------------------------
  # 1) helper functions
  # ---------------------------
  clamp_vec <- function(x, lim = 20) {
    pmax(pmin(x, lim), -lim)
  }
  
  scale01_centered <- function(x, xmin, xmax) {
    # centered scale using window range (same deterministic scaling used in tests)
    (x - (xmin + xmax) / 2) / (xmax - xmin)
  }
  
  build_feature_df <- function(x, y, m_rel = NULL) {
    xc <- (xmin_km + xmax_km) / 2
    yc <- (ymin_km + ymax_km) / 2
    x_std <- scale01_centered(x, xmin_km, xmax_km)
    y_std <- scale01_centered(y, ymin_km, ymax_km)
    r2_center <- (x - xc)^2 + (y - yc)^2
    
    out <- data.frame(
      x = x,
      y = y,
      x_std = x_std,
      y_std = y_std,
      r2_center = r2_center
    )
    if (!is.null(m_rel)) out$m_rel <- m_rel
    out
  }
  
  eval_named_linear_lp <- function(df, coefs, what = "linear LP") {
    if (is.null(coefs) || !is.numeric(coefs) || is.null(names(coefs)) || any(names(coefs) == "")) {
      stop(what, ": 'coefs' must be a named numeric vector.")
    }
    miss <- setdiff(names(coefs), names(df))
    if (length(miss)) {
      stop(what, ": missing variables in feature data frame: ", paste(miss, collapse = ", "))
    }
    as.numeric(as.matrix(df[, names(coefs), drop = FALSE]) %*% as.numeric(coefs))
  }
  
  eval_custom_lp <- function(fun, df, what = "custom LP") {
    if (!is.function(fun)) stop(what, ": function required.")
    val <- fun(df)
    if (!is.numeric(val)) stop(what, ": function must return numeric.")
    if (length(val) == 1L) val <- rep(val, nrow(df))
    if (length(val) != nrow(df)) stop(what, ": function returned wrong length.")
    as.numeric(val)
  }
  
  # Generate/append spatial covariates to a feature data.frame
  append_bg_spatial_covs <- function(df_features, stage = c("support","events")) {
    stage <- match.arg(stage)
    out <- df_features
    
    # generator (returns data.frame)
    if (!is.null(bg_spatial_cov_generator)) {
      if (!is.function(bg_spatial_cov_generator)) {
        stop("bg_spatial_cov_generator must be a function(df_features, window_km, seed=NULL)")
      }
      gen_df <- bg_spatial_cov_generator(
        df_features = df_features,
        window_km = c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km),
        seed = seed
      )
      if (!is.data.frame(gen_df)) stop("bg_spatial_cov_generator must return a data.frame")
      if (nrow(gen_df) != nrow(df_features)) stop("bg_spatial_cov_generator returned wrong nrow")
      if (is.null(names(gen_df)) || any(names(gen_df) == "")) {
        stop("bg_spatial_cov_generator must return named columns")
      }
      # avoid overwriting base feature names silently
      dupn <- intersect(names(gen_df), names(out))
      if (length(dupn)) {
        warning("bg_spatial_cov_generator overwrites columns: ", paste(dupn, collapse = ", "))
        out[, dupn] <- NULL
      }
      out <- cbind(out, gen_df)
    }
    
    # named list of functions
    if (!is.null(bg_spatial_cov_funs)) {
      if (!is.list(bg_spatial_cov_funs) || is.null(names(bg_spatial_cov_funs)) ||
          any(names(bg_spatial_cov_funs) == "")) {
        stop("bg_spatial_cov_funs must be a named list of functions.")
      }
      for (nm in names(bg_spatial_cov_funs)) {
        f <- bg_spatial_cov_funs[[nm]]
        if (!is.function(f)) stop("bg_spatial_cov_funs[['", nm, "']] is not a function.")
        v <- f(out)
        if (!is.numeric(v)) stop("Spatial covariate function '", nm, "' must return numeric.")
        if (length(v) == 1L) v <- rep(v, nrow(out))
        if (length(v) != nrow(out)) stop("Spatial covariate function '", nm, "' returned wrong length.")
        out[[nm]] <- as.numeric(v)
      }
    }
    
    out
  }
  
  # background LP on support/events
  eval_bg_eta <- function(df_features_full) {
    if (bg_lp_type == "constant") {
      eta <- rep(0, nrow(df_features_full))
    } else if (bg_lp_type == "linear") {
      eta <- eval_named_linear_lp(df_features_full, bg_lp_coefs, what = "background linear LP")
    } else if (bg_lp_type == "custom") {
      eta <- eval_custom_lp(bg_lp_fun, df_features_full, what = "background custom LP")
    } else {
      stop("Unsupported bg_lp_type")
    }
    clamp_vec(eta, clip_eta_bg)
  }
  
  # triggering LP for a set of events (event-level features)
  eval_trig_eta <- function(df_event_features) {
    if (trig_lp_type == "linear") {
      eta <- eval_named_linear_lp(df_event_features, trig_lp_coefs, what = "triggering linear LP")
    } else if (trig_lp_type == "custom") {
      eta <- eval_custom_lp(trig_lp_fun, df_event_features, what = "triggering custom LP")
    } else {
      stop("Unsupported trig_lp_type")
    }
    clamp_vec(eta, clip_eta_trig)
  }
  
  # Sample offspring locations from isotropic ETAS kernel
  # Kernel ~ (d + r^2)^(-q), q>1
  sample_xy_etas <- function(n, x0, y0, d, q) {
    theta <- runif(n, 0, 2*pi)
    U <- runif(n)
    R <- sqrt(d * (U^(1/(1 - q)) - 1))
    cbind(x0 + R * cos(theta), y0 + R * sin(theta))
  }
  
  # ---------------------------
  # 2) Build support and background spatial shape
  # ---------------------------
  n_support_obs   <- as.integer(n_support_obs)
  n_support_dummy <- as.integer(mult_support * n_support_obs)
  n_support_tot   <- as.integer(n_support_obs + n_support_dummy)
  
  if (n_support_obs < 1 || n_support_dummy < 1) {
    stop("Need n_support_obs >= 1 and mult_support*n_support_obs >= 1")
  }
  
  # Uniform support points over rectangle
  x_sup_obs <- runif(n_support_obs, xmin_km, xmax_km)
  y_sup_obs <- runif(n_support_obs, ymin_km, ymax_km)
  x_sup_dum <- runif(n_support_dummy, xmin_km, xmax_km)
  y_sup_dum <- runif(n_support_dummy, ymin_km, ymax_km)
  
  x_sup <- c(x_sup_obs, x_sup_dum)
  y_sup <- c(y_sup_obs, y_sup_dum)
  is_dummy_sup <- c(rep(FALSE, n_support_obs), rep(TRUE, n_support_dummy))
  
  # simple uniform quadrature weights for simulation support
  w_sup <- rep(area_km2 / n_support_tot, n_support_tot)
  
  feat_sup_base <- build_feature_df(x_sup, y_sup, m_rel = NULL)
  feat_sup_full <- append_bg_spatial_covs(feat_sup_base, stage = "support")
  
  eta_bg_sup <- eval_bg_eta(feat_sup_full)
  lambda_shape_sup <- exp(eta_bg_sup)         # shape only
  mass_sup <- w_sup * lambda_shape_sup
  mass_tot <- sum(mass_sup)
  
  if (!is.finite(mass_tot) || mass_tot <= 0) stop("Invalid background support mass.")
  
  prob_sup <- mass_sup / mass_tot
  dens_bg_sup <- lambda_shape_sup / sum(w_sup * lambda_shape_sup)  # integrates ~1 wrt w_sup
  
  # ---------------------------
  # 3) Simulate background events (time Poisson + spatial shape)
  # ---------------------------
  muback <- mu * (tmax - tmin)
  n0 <- stats::rpois(1, muback)
  
  # output skeleton
  empty_out <- function() {
    out <- list(
      cat = NULL,
      cat.pois = data.frame(),
      cat.sim = data.frame(),
      n0 = as.integer(0),
      nson = as.integer(0),
      exploded = FALSE,
      truth = list(
        params = params,
        m0 = m0,
        b = b,
        bg = list(
          type = bg_lp_type,
          coefs = if (bg_lp_type == "linear") bg_lp_coefs else NULL,
          note = "Background intercept is not identifiable separately from mu when shape is normalized."
        ),
        trig = list(
          type = trig_lp_type,
          coefs = if (trig_lp_type == "linear") trig_lp_coefs else NULL
        ),
        window_km = c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km),
        t_window = c(tmin = tmin, tmax = tmax)
      ),
      bg_info = NULL
    )
    if (isTRUE(return_bg_info)) {
      sup_df <- cbind(
        feat_sup_full,
        data.frame(
          w = w_sup,
          eta_bg_true = eta_bg_sup,
          lambda_shape = lambda_shape_sup,
          dens_bg = dens_bg_sup,
          is_dummy = is_dummy_sup
        )
      )
      out$bg_info <- list(
        window_km = c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km),
        support = if (isTRUE(return_support)) sup_df else NULL,
        stations = NULL
      )
    }
    out
  }
  
  if (n0 == 0L) return(empty_out())
  
  sample_background_xy_continuous <- function(n_need) {
    out_x <- numeric(0)
    out_y <- numeric(0)
    eta_max <- max(eta_bg_sup, na.rm = TRUE)
    if (!is.finite(eta_max)) eta_max <- 0
    batch <- max(1000L, ceiling(1.5 * n_need))
    guard <- 0L
    while (length(out_x) < n_need) {
      guard <- guard + 1L
      if (guard > 10000L) stop("Background rejection sampler failed to terminate.")
      xcand <- stats::runif(batch, xmin_km, xmax_km)
      ycand <- stats::runif(batch, ymin_km, ymax_km)
      feat_cand <- append_bg_spatial_covs(build_feature_df(xcand, ycand), stage = "events")
      eta_cand <- eval_bg_eta(feat_cand)
      acc <- stats::runif(batch) <= exp(eta_cand - eta_max)
      if (any(acc)) {
        out_x <- c(out_x, xcand[acc])
        out_y <- c(out_y, ycand[acc])
      }
      if (length(out_x) < n_need && guard %% 20L == 0L) batch <- min(batch * 2L, 200000L)
    }
    data.frame(x = out_x[seq_len(n_need)], y = out_y[seq_len(n_need)])
  }
  xy0 <- sample_background_xy_continuous(n0)
  x0 <- xy0$x
  y0 <- xy0$y
  t0 <- runif(n0, tmin, tmax)
  
  beta_GR <- log(10) * b
  m0_bg <- m0 + rexp(n0, rate = beta_GR)
  mrel_bg <- m0_bg - m0
  
  # features and covariates at background events
  feat_bg_events_base <- build_feature_df(x0, y0, m_rel = mrel_bg)
  feat_bg_events_full <- append_bg_spatial_covs(feat_bg_events_base, stage = "events")
  
  eta_bg_bg_events <- eval_bg_eta(feat_bg_events_full)
  eta_trig_bg_events <- eval_trig_eta(feat_bg_events_full)
  
  if (isTRUE(longlat.to.km)) {
    long0 <- x0 / unit
    lat0  <- y0 / unit
  } else {
    long0 <- x0
    lat0  <- y0
  }
  
  # base background catalog
  cat.pois <- data.frame(
    event_id   = seq_len(n0),
    father_id  = 0L,
    lgen       = 0L,
    time       = t0,
    lat        = lat0,
    long       = long0,
    z          = 0,
    magn1      = m0_bg,
    x_km       = x0,
    y_km       = y0,
    m_rel      = mrel_bg,
    x_true     = x0,
    y_true     = y0,
    x_std_true = feat_bg_events_full$x_std,
    y_std_true = feat_bg_events_full$y_std,
    r2_center_true = feat_bg_events_full$r2_center,
    eta_bg_true   = eta_bg_bg_events,
    eta_trig_true = eta_trig_bg_events
  )
  
  # append explicit spatial covariates (NEW)
  bg_cov_names <- setdiff(
    names(feat_bg_events_full),
    c("x","y","x_std","y_std","r2_center","m_rel")
  )
  if (length(bg_cov_names)) {
    for (nm in bg_cov_names) cat.pois[[nm]] <- feat_bg_events_full[[nm]]
  }
  
  cat.pois <- cat.pois[order(cat.pois$time), , drop = FALSE]
  rownames(cat.pois) <- NULL
  
  # reindex event_id after sorting (and keep a map if needed)
  cat.pois$event_id <- seq_len(nrow(cat.pois))
  
  # Working catalog for branching simulation
  cat.new <- cat.pois
  
  # ---------------------------
  # 4) Simulate triggered events (branching)
  # ---------------------------
  # Mean offspring count for each parent (before time/space truncation):
  #   E[N_i] = k0 * exp(gamma*m_rel_i + eta_trig_i) * \int_t * \int_s
  # with integrals on infinite support (then we reject outside window)
  ak <- k0 * c_^(1 - p) / (p - 1)
  sk <- (pi * d^(1 - q)) / (q - 1)
  
  i <- 0L
  exploded <- FALSE
  
  while (i < nrow(cat.new)) {
    i <- i + 1L
    
    # parent features
    m_rel_i <- cat.new$m_rel[i]
    # Recompute triggering eta from current event row to guarantee consistency even for offspring
    feat_parent <- data.frame(
      x = cat.new$x_km[i],
      y = cat.new$y_km[i],
      x_std = cat.new$x_std_true[i],
      y_std = cat.new$y_std_true[i],
      r2_center = cat.new$r2_center_true[i],
      m_rel = cat.new$m_rel[i]
    )
    # add spatial covariates if present in row
    if (length(bg_cov_names)) {
      for (nm in bg_cov_names) feat_parent[[nm]] <- cat.new[[nm]][i]
    }
    
    eta_trig_i <- eval_trig_eta(feat_parent)
    # Store/overwrite eta_trig_true for parent row (for robustness)
    cat.new$eta_trig_true[i] <- eta_trig_i
    
    n_exp_i <- ak * sk * exp(gamma * m_rel_i + eta_trig_i)
    if (!is.finite(n_exp_i) || n_exp_i < 0) n_exp_i <- 0
    
    ni <- stats::rpois(1, lambda = n_exp_i)
    if (ni > 0L) {
      # offspring times
      t_child <- c_ * runif(ni)^(-1/(p - 1)) - c_ + cat.new$time[i]
      
      # offspring locations
      xy_child <- sample_xy_etas(ni, x0 = cat.new$x_km[i], y0 = cat.new$y_km[i], d = d, q = q)
      
      inside <- (t_child > tmin) & (t_child < tmax) &
        (xy_child[,1] > xmin_km) & (xy_child[,1] < xmax_km) &
        (xy_child[,2] > ymin_km) & (xy_child[,2] < ymax_km)
      
      if (any(inside)) {
        nt <- sum(inside)
        x1 <- xy_child[inside, 1]
        y1 <- xy_child[inside, 2]
        t1 <- t_child[inside]
        m1 <- m0 + rexp(nt, rate = beta_GR)
        mrel1 <- m1 - m0
        
        # event-level features/covariates
        feat_child_base <- build_feature_df(x1, y1, m_rel = mrel1)
        feat_child_full <- append_bg_spatial_covs(feat_child_base, stage = "events")
        eta_bg_child    <- eval_bg_eta(feat_child_full)
        eta_trig_child  <- eval_trig_eta(feat_child_full)
        
        if (isTRUE(longlat.to.km)) {
          long1 <- x1 / unit
          lat1  <- y1 / unit
        } else {
          long1 <- x1
          lat1  <- y1
        }
        
        child_df <- data.frame(
          event_id   = NA_integer_,     # assigned after sorting
          father_id  = as.integer(cat.new$event_id[i]),
          lgen       = as.integer(cat.new$lgen[i] + 1L),
          time       = t1,
          lat        = lat1,
          long       = long1,
          z          = 0,
          magn1      = m1,
          x_km       = x1,
          y_km       = y1,
          m_rel      = mrel1,
          x_true     = x1,
          y_true     = y1,
          x_std_true = feat_child_full$x_std,
          y_std_true = feat_child_full$y_std,
          r2_center_true = feat_child_full$r2_center,
          eta_bg_true   = eta_bg_child,
          eta_trig_true = eta_trig_child
        )
        
        if (length(bg_cov_names)) {
          for (nm in bg_cov_names) child_df[[nm]] <- feat_child_full[[nm]]
        }
        
        cat.new <- rbind_fill_base(cat.new, child_df)
        
        # safety
        if (isTRUE(explode_guard) && nrow(cat.new) > max_events) {
          warning("Maximum number of events exceeded (max_events). Stopping branching simulation.")
          exploded <- TRUE
          break
        }
        
        # Keep time-sorted and reindex event_id
        cat.new <- cat.new[order(cat.new$time), , drop = FALSE]
        rownames(cat.new) <- NULL
        old_ids <- cat.new$event_id
        cat.new$event_id <- seq_len(nrow(cat.new))
        
        # Update father_id references after reordering/reindexing:
        # old_ids contains previous IDs (including NA for just-added offspring)
        # Need map from old id -> new id for existing events.
        valid_old <- !is.na(old_ids)
        id_map <- integer(max(old_ids[valid_old], na.rm = TRUE))
        id_map[old_ids[valid_old]] <- which(valid_old)
        
        # father_id == 0 remains 0, others remapped
        idx_f <- which(cat.new$father_id > 0)
        if (length(idx_f)) {
          cat.new$father_id[idx_f] <- id_map[cat.new$father_id[idx_f]]
        }
      }
    }
    
    if (exploded) break
  }
  
  # final outputs
  nson <- nrow(cat.new) - n0
  
  # refresh cat.pois as subset father_id==0 after possible reindexing
  cat.pois_final <- cat.new[cat.new$father_id == 0, , drop = FALSE]
  cat.pois_final <- cat.pois_final[order(cat.pois_final$time), , drop = FALSE]
  rownames(cat.pois_final) <- NULL
  
  # truth object
  truth_obj <- list(
    params = params,
    m0 = m0,
    b = b,
    bg = list(
      type = bg_lp_type,
      coefs = if (bg_lp_type == "linear") bg_lp_coefs else NULL,
      note = "Background intercept is not identifiable separately from mu when shape is normalized."
    ),
    trig = list(
      type = trig_lp_type,
      coefs = if (trig_lp_type == "linear") trig_lp_coefs else NULL
    ),
    window_km = c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km),
    t_window = c(tmin = tmin, tmax = tmax)
  )
  
  bg_info <- NULL
  if (isTRUE(return_bg_info)) {
    sup_df <- cbind(
      feat_sup_full,
      data.frame(
        w = w_sup,
        eta_bg_true = eta_bg_sup,
        lambda_shape = lambda_shape_sup,
        dens_bg = dens_bg_sup,
        is_dummy = is_dummy_sup
      )
    )
    bg_info <- list(
      window_km = c(xmin = xmin_km, xmax = xmax_km, ymin = ymin_km, ymax = ymax_km),
      support = if (isTRUE(return_support)) sup_df else NULL,
      stations = NULL
    )
  }
  
  # optional plots
  if (isTRUE(plot_bg) && isTRUE(return_bg_info) && !is.null(bg_info) && !is.null(bg_info$support)) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    sup <- bg_info$support
    pal <- colorRampPalette(c("white", "yellow", "orange", "red"))(100)
    dens01 <- sup$dens_bg / max(sup$dens_bg, na.rm = TRUE)
    dens01[!is.finite(dens01)] <- 0
    colv <- pal[pmax(1, pmin(100, 1 + floor(99 * dens01)))]
    
    par(mfrow = c(1,2), mar = c(4,4,2,1))
    plot(sup$x, sup$y, pch = 16, cex = 0.35, col = colv,
         xlab = "x (km)", ylab = "y (km)",
         main = paste0("Background support (", bg_lp_type, ")"))
    hist(sup$eta_bg_true, breaks = 40, main = "eta_bg_true on support", xlab = "eta_bg_true")
  }
  
  if (isTRUE(plot_catalog)) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mfrow = c(1,1), mar = c(4,4,2,1))
    is_bg <- (cat.new$father_id == 0)
    plot(cat.new$long, cat.new$lat,
         col = ifelse(is_bg, "blue", "red"),
         pch = 16, cex = cex_cat,
         xlab = "Longitude", ylab = "Latitude",
         main = "Simulated catalog (blue=background, red=triggered)")
    legend("topright", legend = c("background", "triggered"),
           col = c("blue","red"), pch = 16, bty = "n")
  }
  
  out <- list(
    cat = NULL,
    cat.pois = cat.pois_final,
    cat.sim  = cat.new,
    n0 = as.integer(n0),
    nson = as.integer(nson),
    exploded = isTRUE(exploded),
    truth = truth_obj
  )
  if (isTRUE(return_bg_info)) out$bg_info <- bg_info
  
  out
}





##################################################

#######################################################################
# 3) PIANO DI SIMULAZIONE POSTER
#######################################################################

## ================================================================
## 2) CONFIGURAZIONE POSTER
## ================================================================

## Finestra: se catalog.withcov esiste, usa il catalogo reale; altrimenti fallback Italia.
if (!exists("catalog.withcov")) {
  try(utils::data("italycatalog", package = "etasFLP"), silent = TRUE)
}
if (exists("catalog.withcov")) {
  LONG_RANGE <- range(catalog.withcov$long, na.rm = TRUE)
  LAT_RANGE  <- range(catalog.withcov$lat,  na.rm = TRUE)
  TMIN_USE   <- min(catalog.withcov$time, na.rm = TRUE)
} else {
  LONG_RANGE <- c(6, 20)
  LAT_RANGE  <- c(36, 45)
  TMIN_USE   <- 41033
}

T_LAG_USE <- 700

TRUE_BASE <- c(
  mu    = NA_real_,
  k0    = NA_real_,
  c     = 0.010,
  p     = 1.12,
  gamma = 0.00,
  d     = 1.60,
  q     = 1.85
)

TRUE_BETACOV <- 0.70
M0_TRUE      <- 2.5
B_TRUE       <- 1.05

## Piano poster
N_REP_MAIN_POSTER    <- 30L
N_REP_PILOT_POSTER   <- 20L
N_CALIB_ITERS_POSTER <- 8L
N_CALIB_BRACKET_POSTER <- 8L
SEED_MASTER_POSTER   <- 20260504L

MAX_EVENTS <- 50000L

NDECLUST_FIT <- 4L
ITERLIM_FIT  <- 40L
MULT_BG_FIT  <- 4L
K_SPLINE_BG  <- c(12L, 12L)

TOL_BG_REL   <- 0.07
TOL_TRIG_REL <- 0.10

## Coefficienti veri background
TRUE_BETA_XSTD <- 0.90
TRUE_BETA_YSTD <- -0.70
TRUE_BETA_Z1   <- 0.80

## ================================================================
## 3) FUNZIONI DI FEATURE/COVARIATE BACKGROUND
## ================================================================

.window_km <- function(long.range = LONG_RANGE, lat.range = LAT_RANGE, longlat.to.km = TRUE) {
  unit <- if (isTRUE(longlat.to.km)) 6371.3 * pi / 180 else 1
  c(
    xmin = long.range[1] * unit,
    xmax = long.range[2] * unit,
    ymin = lat.range[1]  * unit,
    ymax = lat.range[2]  * unit
  )
}

background_features_from_xy <- function(x, y, window_km = .window_km()) {
  xc <- (window_km["xmin"] + window_km["xmax"]) / 2
  yc <- (window_km["ymin"] + window_km["ymax"]) / 2
  data.frame(
    x = x,
    y = y,
    x_std = (x - xc) / (window_km["xmax"] - window_km["xmin"]),
    y_std = (y - yc) / (window_km["ymax"] - window_km["ymin"]),
    r2_center = (x - xc)^2 + (y - yc)^2
  )
}

bg_cov_Z1_truth <- function(df_feat) {
  sin(pi * df_feat$x_std) * cos(pi * df_feat$y_std)
}

bg_fun_smooth_truth <- function(df_feat) {
  1.1 * sin(pi * df_feat$x_std) -
    0.9 * cos(pi * df_feat$y_std) +
    0.6 * (df_feat$x_std * df_feat$y_std)
}

exact_covariate_functions <- function(bg_case, window_km = .window_km()) {
  f_xstd <- function(df) background_features_from_xy(df$x, df$y, window_km)$x_std
  f_ystd <- function(df) background_features_from_xy(df$x, df$y, window_km)$y_std
  f_Z1 <- function(df) {
    ff <- background_features_from_xy(df$x, df$y, window_km)
    bg_cov_Z1_truth(ff)
  }

  if (bg_case == "linear_xy") {
    list(x_std = f_xstd, y_std = f_ystd)
  } else if (bg_case == "cov1") {
    list(Z1 = f_Z1)
  } else {
    NULL
  }
}

BACKGROUND_SCENARIOS_POSTER <- list(
  linear_xy = list(
    bg_case = "linear_xy",
    label = "Linear background: x + y",
    sim_spec = list(
      bg_lp_type  = "linear",
      bg_lp_coefs = c(x_std = TRUE_BETA_XSTD, y_std = TRUE_BETA_YSTD)
    ),
    has_bg_coef = TRUE,
    true_bg_coef = c(x_std = TRUE_BETA_XSTD, y_std = TRUE_BETA_YSTD)
  ),

  cov1 = list(
    bg_case = "cov1",
    label = "Covariate-driven background: Z1",
    sim_spec = list(
      bg_lp_type = "linear",
      bg_lp_coefs = c(Z1 = TRUE_BETA_Z1),
      bg_spatial_cov_funs = list(Z1 = bg_cov_Z1_truth)
    ),
    has_bg_coef = TRUE,
    true_bg_coef = c(Z1 = TRUE_BETA_Z1)
  ),

  smooth = list(
    bg_case = "smooth",
    label = "Smooth background: f(x,y)",
    sim_spec = list(
      bg_lp_type = "custom",
      bg_lp_fun = bg_fun_smooth_truth
    ),
    has_bg_coef = FALSE,
    true_bg_coef = numeric(0)
  )
)

make_count_regimes_poster <- function(T_lag = T_LAG_USE,
                                      c_true = TRUE_BASE["c"],
                                      p_true = TRUE_BASE["p"],
                                      d_true = TRUE_BASE["d"],
                                      q_true = TRUE_BASE["q"],
                                      beta_mag_true = TRUE_BETACOV,
                                      b_true = B_TRUE) {
  grid <- data.frame(
    count_case = c("N1000_bal", "N1000_trigdom"),
    N_tot_target = c(1000, 1000),
    r_trig_bg = c(1, 3),
    stringsAsFactors = FALSE
  )

  grid$N_bg_target <- with(grid, N_tot_target / (1 + r_trig_bg))
  grid$N_trig_target <- with(grid, N_tot_target - N_bg_target)
  grid$mu_init <- grid$N_bg_target / T_lag

  beta_GR <- log(10) * b_true
  if (beta_mag_true >= beta_GR) {
    stop("TRUE_BETACOV deve essere < log(10)*B_TRUE per avere E[exp(beta*M_rel)] finito.")
  }
  EexpM <- beta_GR / (beta_GR - beta_mag_true)

  C_k0 <- (c_true^(1 - p_true) / (p_true - 1)) *
    (pi * d_true^(1 - q_true) / (q_true - 1)) *
    EexpM

  grid$nbar_target <- with(grid, r_trig_bg / (1 + r_trig_bg))
  grid$k0_init <- grid$nbar_target / C_k0
  grid$T_lag <- T_lag
  grid
}

COUNT_REGIMES_POSTER <- make_count_regimes_poster()

## ================================================================
## 4) SIMULAZIONE, CALIBRAZIONE, FIT
## ================================================================

get_sim_core <- function(obj) {
  
  ## Se per caso viene passato un safe_run object completo
  if (is.list(obj) && !is.null(obj$ok) && !is.null(obj$value)) {
    obj <- obj$value
  }
  
  ## Se la simulazione è dentro $sim
  if (is.list(obj) && !is.null(obj$sim) && is.list(obj$sim)) {
    obj <- obj$sim
  }
  
  obj
}

get_cat_sim <- function(sim_obj) {
  
  s <- get_sim_core(sim_obj)
  
  if (!is.null(s$cat.sim) && is.data.frame(s$cat.sim)) {
    return(s$cat.sim)
  }
  
  if (!is.null(s$catalogue) && is.data.frame(s$catalogue)) {
    return(s$catalogue)
  }
  
  if (!is.null(s$cat) && is.data.frame(s$cat)) {
    return(s$cat)
  }
  
  stop(
    "Catalogo simulato non trovato. Nomi disponibili dopo get_sim_core(): ",
    paste(names(s), collapse = ", ")
  )
}

get_bg_info <- function(sim_obj) {
  get_sim_core(sim_obj)$bg_info
}

get_counts <- function(sim_obj) {
  s <- get_sim_core(sim_obj)
  c(n0 = s$n0, nson = s$nson, exploded = isTRUE(s$exploded))
}

make_eqcat_from_sim <- function(sim_obj) {
  cat_df <- get_cat_sim(sim_obj)
  need <- c("time", "lat", "long", "z", "magn1")
  miss <- setdiff(need, names(cat_df))
  if (length(miss)) stop("Mancano colonne per eqcat: ", paste(miss, collapse = ", "))
  cat_df <- cat_df[order(cat_df$time), , drop = FALSE]
  rownames(cat_df) <- NULL
  cat_df
}

simulate_one_from_cell_poster <- function(bg_scenario, count_regime, seed,
                                          mu_value, k0_value,
                                          return_bg_info = TRUE,
                                          return_support = TRUE,
                                          plot_bg = FALSE,
                                          plot_catalog = FALSE) {
  params_now <- TRUE_BASE
  params_now["mu"] <- mu_value
  params_now["k0"] <- k0_value

  sim_args <- list(
    params = params_now,
    m0 = M0_TRUE,
    b = B_TRUE,
    tmin = TMIN_USE,
    t.lag = count_regime$T_lag,
    long.range = LONG_RANGE,
    lat.range = LAT_RANGE,
    longlat.to.km = TRUE,
    sectoday = FALSE,
    trig_lp_type = "linear",
    trig_lp_coefs = c(m_rel = TRUE_BETACOV),
    n_support_obs = 1500,
    mult_support = 4,
    max_events = MAX_EVENTS,
    return_bg_info = return_bg_info,
    return_support = return_support,
    plot_bg = plot_bg,
    plot_catalog = plot_catalog,
    seed = seed
  )

  sim_args <- c(sim_args, bg_scenario$sim_spec)
  do.call(etas.par.sim_v4, sim_args)
}

calibrate_cell_poster <- function(bg_scenario, count_regime,
                                  nrep_pilot = N_REP_PILOT_POSTER,
                                  n_iters = N_CALIB_ITERS_POSTER,
                                  n_bracket = N_CALIB_BRACKET_POSTER,
                                  seed_base = 10000L,
                                  verbose = TRUE,
                                  max_exploded_rate = 0.20) {
  ## Calibrazione robusta per il piano poster.
  ##
  ## Nota metodologica:
  ## - mu controlla direttamente il numero medio di eventi background,
  ##   quindi lo teniamo vicino a target_bg / T_lag;
  ## - k0 controlla la branching/triggering cascade. Il numero totale di
  ##   discendenti e' molto non lineare in k0, soprattutto nel caso
  ##   triggered-dominant. Per questo NON usiamo piu' aggiornamenti
  ##   moltiplicativi iterativi, ma una ricerca a bracket + bisezione
  ##   sul numero medio simulato di triggered events.

  mu_cur <- count_regime$mu_init
  target_bg <- count_regime$N_bg_target
  target_trig <- count_regime$N_trig_target

  eval_id <- 0L
  logs <- list()

  eval_k0 <- function(k0_value, phase, iter) {
    eval_id <<- eval_id + 1L
    n0_vec <- numeric(nrep_pilot)
    nson_vec <- numeric(nrep_pilot)
    exploded_vec <- logical(nrep_pilot)

    for (r in seq_len(nrep_pilot)) {
      seed_r <- seed_base + 100000L * eval_id + r
      sim_r <- try(
        simulate_one_from_cell_poster(
          bg_scenario = bg_scenario,
          count_regime = count_regime,
          seed = seed_r,
          mu_value = mu_cur,
          k0_value = k0_value,
          return_bg_info = FALSE,
          return_support = FALSE
        ),
        silent = TRUE
      )

      if (inherits(sim_r, "try-error")) {
        n0_vec[r] <- NA_real_
        nson_vec[r] <- NA_real_
        exploded_vec[r] <- TRUE
      } else {
        cc <- get_counts(sim_r)
        n0_vec[r] <- cc["n0"]
        nson_vec[r] <- cc["nson"]
        exploded_vec[r] <- isTRUE(cc["exploded"])
      }
    }

    mean_n0 <- mean(n0_vec, na.rm = TRUE)
    mean_nson <- mean(nson_vec, na.rm = TRUE)
    exploded_rate <- mean(exploded_vec, na.rm = TRUE)
    rel_err_bg <- (mean_n0 - target_bg) / target_bg
    rel_err_trig <- (mean_nson - target_trig) / target_trig

    row <- data.frame(
      eval_id = eval_id,
      phase = phase,
      iter = iter,
      bg_case = bg_scenario$bg_case,
      count_case = count_regime$count_case,
      mu = mu_cur,
      k0 = k0_value,
      mean_n0 = mean_n0,
      mean_nson = mean_nson,
      target_bg = target_bg,
      target_trig = target_trig,
      rel_err_bg = rel_err_bg,
      rel_err_trig = rel_err_trig,
      exploded_rate = exploded_rate
    )
    logs[[length(logs) + 1L]] <<- row

    if (isTRUE(verbose)) {
      cat(sprintf(
        "[CALIB:%s] %-12s | %-14s | iter %02d | mu=%.5f k0=%.6f | n0 %.1f/%.1f | nson %.1f/%.1f | rel %.2f | expl %.2f\n",
        phase, bg_scenario$bg_case, count_regime$count_case, iter,
        mu_cur, k0_value, mean_n0, target_bg, mean_nson, target_trig,
        rel_err_trig, exploded_rate
      ))
    }
    row
  }

  ## Piccolo controllo su mu: il numero di background e' Poisson(mu*T_lag).
  ## Usiamo il valore analitico target_bg/T_lag; questo evita che la
  ## calibrazione di k0 venga contaminata da rumore su mu.
  mu_cur <- max(target_bg / count_regime$T_lag, 1e-8)

  k0_guess <- max(count_regime$k0_init, 1e-8)
  low_k0 <- 0
  high_k0 <- k0_guess

  high_eval <- eval_k0(high_k0, phase = "bracket", iter = 1L)

  ## Allarga il bracket finche' la media triggered supera il target oppure
  ## finche' compare instabilita'/explosion. In entrambi i casi high_k0 e'
  ## considerato un upper bracket.
  b_iter <- 1L
  while (is.finite(high_eval$mean_nson) &&
         high_eval$mean_nson < target_trig &&
         high_eval$exploded_rate <= max_exploded_rate &&
         b_iter < n_bracket) {
    low_k0 <- high_k0
    high_k0 <- high_k0 * 1.45
    b_iter <- b_iter + 1L
    high_eval <- eval_k0(high_k0, phase = "bracket", iter = b_iter)
  }

  ## Se anche dopo il bracket non si supera il target, scegliamo il migliore
  ## tra i valori valutati e segnaliamo la cosa nel log.
  if (!is.finite(high_eval$mean_nson) || high_eval$exploded_rate > max_exploded_rate) {
    ## high_k0 rimane upper perché instabile/exploded.
  } else if (high_eval$mean_nson < target_trig) {
    if (isTRUE(verbose)) {
      cat("[CALIB WARNING] Non sono riuscito a bracketare il target triggered; seleziono il migliore tra i valori provati.\n")
    }
  }

  ## Bisezione. Se il midpoint produce meno triggered del target, alza low;
  ## se ne produce troppi o esplode, abbassa high.
  for (it in seq_len(n_iters)) {
    mid_k0 <- (low_k0 + high_k0) / 2
    mid_eval <- eval_k0(mid_k0, phase = "bisect", iter = it)

    if (!is.finite(mid_eval$mean_nson) || mid_eval$exploded_rate > max_exploded_rate) {
      high_k0 <- mid_k0
    } else if (mid_eval$mean_nson < target_trig) {
      low_k0 <- mid_k0
    } else {
      high_k0 <- mid_k0
    }

    if (is.finite(mid_eval$rel_err_trig) && abs(mid_eval$rel_err_trig) < TOL_TRIG_REL &&
        is.finite(mid_eval$rel_err_bg) && abs(mid_eval$rel_err_bg) < TOL_BG_REL) {
      break
    }
  }

  pilot_log <- do.call(rbind_fill_base, logs)

  ## Scegli il valore migliore valutato, penalizzando explosion.
  cand <- pilot_log[is.finite(pilot_log$mean_nson) & pilot_log$exploded_rate <= max_exploded_rate, , drop = FALSE]
  if (!nrow(cand)) {
    cand <- pilot_log[is.finite(pilot_log$mean_nson), , drop = FALSE]
  }
  if (!nrow(cand)) {
    stop("Calibrazione fallita: nessun valore finito di mean_nson.")
  }
  cand$abs_rel_err_trig <- abs((cand$mean_nson - target_trig) / target_trig)
  cand$abs_rel_err_bg <- abs((cand$mean_n0 - target_bg) / target_bg)
  cand$score <- cand$abs_rel_err_trig + 0.25 * cand$abs_rel_err_bg + 2 * cand$exploded_rate
  best <- cand[which.min(cand$score), , drop = FALSE]

  if (isTRUE(verbose)) {
    cat(sprintf(
      "[CALIB:SELECT] %-12s | %-14s | mu=%.5f k0=%.6f | n0 %.1f/%.1f | nson %.1f/%.1f | rel %.2f | expl %.2f\n",
      bg_scenario$bg_case, count_regime$count_case,
      best$mu, best$k0, best$mean_n0, target_bg, best$mean_nson, target_trig,
      best$rel_err_trig, best$exploded_rate
    ))
  }

  list(mu_cal = as.numeric(best$mu), k0_cal = as.numeric(best$k0), pilot_log = pilot_log)
}

fit_etas_classic_wrapper_poster <- function(cat_sim_df,
                                            m0_true = M0_TRUE,
                                            starts = NULL,
                                            ndeclust = NDECLUST_FIT,
                                            iterlim = ITERLIM_FIT) {
  if (is.null(starts)) {
    starts <- list(mu = 0.4, k0 = 0.01, c = 0.02, p = 1.2, gamma = 0, d = 1.2, q = 1.7, betacov = 0.5)
  }

  etasclass(
    cat.orig = cat_sim_df,
    magn.threshold = m0_true,
    magn.threshold.back = m0_true,
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
    flp = FALSE,
    ndeclust = ndeclust,
    onlytime = FALSE,
    is.backconstant = FALSE,
    sectoday = FALSE,
    longlat.to.km = TRUE,
    usenlm = TRUE,
    compsqm = TRUE,
    epsmax = 1e-4,
    iterlim = iterlim,
    ntheta = 36
  )
}

fit_etas_parametric_wrapper_poster <- function(cat_sim_df,
                                               m0_true = M0_TRUE,
                                               bg_case,
                                               starts = NULL,
                                               ndeclust = NDECLUST_FIT,
                                               iterlim = ITERLIM_FIT,
                                               mult.bg = MULT_BG_FIT,
                                               verbose.bg = FALSE) {
  if (is.null(starts)) {
    starts <- list(mu = 0.4, k0 = 0.01, c = 0.02, p = 1.2, gamma = 0, d = 1.2, q = 1.7, betacov = 0.5)
  }

  covs.bg.user <- exact_covariate_functions(bg_case, .window_km())

  if (bg_case == "linear_xy") {
    formula.bg <- ~ x_std + y_std
    spatial.cov.bg <- TRUE
    type.cov.values.bg <- c(x_std = "exact", y_std = "exact")
  } else if (bg_case == "cov1") {
    formula.bg <- ~ Z1
    spatial.cov.bg <- TRUE
    type.cov.values.bg <- c(Z1 = "exact")
  } else if (bg_case == "smooth") {
    formula.bg <- as.formula(sprintf("~ te(x, y, k = c(%d,%d))", K_SPLINE_BG[1], K_SPLINE_BG[2]))
    spatial.cov.bg <- FALSE
    type.cov.values.bg <- NULL
    covs.bg.user <- NULL
  } else {
    stop("bg_case non riconosciuto: ", bg_case)
  }

  etasclass.par(
    cat.orig = cat_sim_df,
    time.update = FALSE,
    magn.threshold = m0_true,
    magn.threshold.back = m0_true,
    tmax = max(cat_sim_df$time),
    long.range = LONG_RANGE,
    lat.range = LAT_RANGE,
    mu = starts$mu,
    k0 = starts$k0,
    c = starts$c,
    p = starts$p,
    gamma = 0.0,
    d = starts$d,
    q = max(starts$q, 1.01),
    betacov = starts$betacov,
    params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    formula1 = "time ~ magnitude - 1",
    offset = 0,
    hdef = c(1, 1),
    wp = rep(1, nrow(cat_sim_df)),
    hvarx = rep(1, nrow(cat_sim_df)),
    hvary = rep(1, nrow(cat_sim_df)),
    declustering = TRUE,
    thinning = FALSE,
    flp = FALSE,
    m1 = NULL,
    ndeclust = ndeclust,
    n.iterweight = 1,
    onlytime = FALSE,
    is.backconstant = FALSE,
    description = "",
    cat.back = NULL,
    back.smooth = 1,
    sectoday = FALSE,
    longlat.to.km = TRUE,
    usenlm = TRUE,
    method = "BFGS",
    compsqm = TRUE,
    epsmax = 1e-4,
    iterlim = iterlim,
    ntheta = 36,
    formula.bg = formula.bg,
    process.type.bg = "s2d",
    verbose.bg = verbose.bg,
    marked.bg = FALSE,
    mark.c.bg = FALSE,
    type.cov.values.bg = type.cov.values.bg,
    grid.bg = FALSE,
    mult.bg = mult.bg,
    ncube.bg = NULL,
    spatial.cov.bg = spatial.cov.bg,
    offset_k.bg = FALSE,
    seed.bg = NULL,
    covs.bg.user = covs.bg.user
  )
}

## ================================================================
## 5) METRICHE
## ================================================================

rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))

safe_cor <- function(x, y) {
  out <- suppressWarnings(stats::cor(x, y, use = "complete.obs"))
  if (is.na(out)) NA_real_ else out
}

compute_true_event_intensity <- function(sim_obj) {
  s <- get_sim_core(sim_obj)
  catsim <- s$cat.sim[order(s$cat.sim$time), , drop = FALSE]
  rownames(catsim) <- NULL
  truth <- s$truth
  bginfo <- s$bg_info

  req_ev <- c("time", "x_km", "y_km", "m_rel", "eta_bg_true", "eta_trig_true")
  miss_ev <- setdiff(req_ev, names(catsim))
  if (length(miss_ev)) stop("cat.sim manca: ", paste(miss_ev, collapse = ", "))

  par <- truth$params
  mu_true <- unname(par["mu"])
  k0_true <- unname(par["k0"])
  c_true <- unname(par["c"])
  p_true <- unname(par["p"])
  gamma_true <- unname(par["gamma"])
  d_true <- unname(par["d"])
  q_true <- unname(par["q"])

  sup <- bginfo$support
  if (is.null(sup)) stop("support background mancante")
  kappa_shape <- if ("lambda_shape" %in% names(sup)) {
    sum(sup$w * sup$lambda_shape, na.rm = TRUE)
  } else {
    sum(sup$w * exp(sup$eta_bg_true), na.rm = TRUE)
  }

  lambda_bg_true <- mu_true * exp(catsim$eta_bg_true) / kappa_shape

  n <- nrow(catsim)
  lambda_trig_true <- numeric(n)
  tvec <- catsim$time
  xvec <- catsim$x_km
  yvec <- catsim$y_km
  mrel <- catsim$m_rel
  eta_prod <- catsim$eta_trig_true

  for (i in seq_len(n)) {
    if (i == 1L) next
    dt <- tvec[i] - tvec[seq_len(i - 1L)]
    j <- which(dt > 0)
    if (!length(j)) next
    dx <- xvec[i] - xvec[j]
    dy <- yvec[i] - yvec[j]
    r2 <- dx^2 + dy^2
    contrib <- k0_true * exp(gamma_true * mrel[j] + eta_prod[j]) *
      (c_true + dt[j])^(-p_true) * (d_true + r2)^(-q_true)
    lambda_trig_true[i] <- sum(contrib, na.rm = TRUE)
  }

  data.frame(
    lambda_bg_true = lambda_bg_true,
    eta_trig_true = lambda_trig_true,
    lambda_tot_true = lambda_bg_true + lambda_trig_true
  )
}

extract_fit_trigger_params <- function(fit_obj) {
  est <- fit_obj$params.MLtot
  idx_mag <- grep("magnitude", names(est), fixed = TRUE)
  betacov_hat <- if (length(idx_mag)) unname(est[idx_mag[1]]) else NA_real_
  c(
    mu = unname(est["mu"]),
    k0 = unname(est["k0"]),
    c = unname(est["c"]),
    p = unname(est["p"]),
    gamma = unname(est["gamma"]),
    d = unname(est["d"]),
    q = unname(est["q"]),
    betacov = betacov_hat
  )
}

extract_fit_event_intensity <- function(fit_obj) {
  if (is.null(fit_obj$l)) stop("fit_obj$l mancante")
  if (is.null(fit_obj$back.dens)) stop("fit_obj$back.dens mancante")
  if (is.null(fit_obj$params.MLtot)) stop("fit_obj$params.MLtot mancante")
  mu_hat <- unname(fit_obj$params.MLtot["mu"])
  lambda_tot_hat <- as.numeric(fit_obj$l)
  lambda_bg_hat <- mu_hat * as.numeric(fit_obj$back.dens)
  data.frame(
    lambda_bg_hat = lambda_bg_hat,
    eta_trig_hat = lambda_tot_hat - lambda_bg_hat,
    lambda_tot_hat = lambda_tot_hat
  )
}

compute_intensity_metrics <- function(true_int, fit_int) {
  stopifnot(nrow(true_int) == nrow(fit_int))
  eps <- 1e-12
  c(
    rmse_lambda_tot = rmse(fit_int$lambda_tot_hat - true_int$lambda_tot_true),
    rmse_lambda_bg = rmse(fit_int$lambda_bg_hat - true_int$lambda_bg_true),
    rmse_eta_trig = rmse(fit_int$eta_trig_hat - true_int$eta_trig_true),
    rmse_log_lambda_tot = rmse(log(fit_int$lambda_tot_hat + eps) - log(true_int$lambda_tot_true + eps)),
    rmse_log_lambda_bg = rmse(log(pmax(fit_int$lambda_bg_hat, 0) + eps) - log(true_int$lambda_bg_true + eps)),
    rmse_log_eta_trig = rmse(log(pmax(fit_int$eta_trig_hat, 0) + eps) - log(true_int$eta_trig_true + eps)),
    cor_log_lambda_tot = safe_cor(log(fit_int$lambda_tot_hat + eps), log(true_int$lambda_tot_true + eps)),
    cor_log_lambda_bg = safe_cor(log(pmax(fit_int$lambda_bg_hat, 0) + eps), log(true_int$lambda_bg_true + eps)),
    cor_log_eta_trig = safe_cor(log(pmax(fit_int$eta_trig_hat, 0) + eps), log(true_int$eta_trig_true + eps))
  )
}

compare_background_shape_param <- function(sim_obj, fit_param_obj) {
  s <- get_sim_core(sim_obj)
  sup <- s$bg_info$support
  if (!all(c("x", "y", "w", "dens_bg") %in% names(sup))) {
    stop("support simulazione non contiene x,y,w,dens_bg")
  }
  lambda_hat <- as.numeric(stats::predict(fit_param_obj$model.bg$mod_global, newdata = sup, type = "response"))
  kappa_hat <- sum(sup$w * lambda_hat, na.rm = TRUE)
  dens_hat <- lambda_hat / kappa_hat
  dens_true <- sup$dens_bg
  w <- sup$w
  eps <- 1e-20
  log_true <- log(dens_true + eps)
  log_hat <- log(dens_hat + eps)
  log_true_c <- log_true - stats::weighted.mean(log_true, w = w)
  log_hat_c <- log_hat - stats::weighted.mean(log_hat, w = w)
  c(
    rmse_dens = sqrt(stats::weighted.mean((dens_hat - dens_true)^2, w = w)),
    mae_dens = stats::weighted.mean(abs(dens_hat - dens_true), w = w),
    corr_log_shape = safe_cor(log_true_c, log_hat_c),
    rmse_log_shape = sqrt(stats::weighted.mean((log_hat_c - log_true_c)^2, w = w))
  )
}

extract_param_bg_coefs <- function(fit_param_obj, coef_names) {
  cf <- stats::coef(fit_param_obj$model.bg$mod_global)
  out <- setNames(rep(NA_real_, length(coef_names)), coef_names)
  for (nm in coef_names) {
    ix <- which(names(cf) == nm)
    if (!length(ix)) ix <- grep(nm, names(cf), fixed = TRUE)
    if (length(ix)) out[nm] <- unname(cf[ix[1]])
  }
  out
}

trigger_param_names7 <- c("mu", "k0", "c", "p", "d", "q", "betacov")

INT_METRIC_NAMES <- c(
  "rmse_lambda_tot", "rmse_lambda_bg", "rmse_eta_trig",
  "rmse_log_lambda_tot", "rmse_log_lambda_bg", "rmse_log_eta_trig",
  "cor_log_lambda_tot", "cor_log_lambda_bg", "cor_log_eta_trig"
)

empty_int_metrics <- function() {
  x <- rep(NA_real_, length(INT_METRIC_NAMES))
  names(x) <- INT_METRIC_NAMES
  x
}

compute_trigger_errors_one <- function(est_vec, true_vec_named, true_betacov) {
  truth <- c(true_vec_named[c("mu", "k0", "c", "p", "d", "q")], betacov = true_betacov)
  est <- c(est_vec[c("mu", "k0", "c", "p", "d", "q")], betacov = est_vec["betacov"])
  est - truth
}

rmse_trigger7_one <- function(err_named) {
  sqrt(mean(as.numeric(err_named[trigger_param_names7])^2, na.rm = TRUE))
}

## ================================================================
## 6) FIT CON CHECKPOINT DOPO OGNI FIT
## ================================================================

run_fit_stage <- function(model_name, expr_fun, rep_dir, metadata,
                          progress_env = NULL, rerun_failed = FALSE) {
  fit_file <- file.path(rep_dir, paste0("fit_", model_name, ".rds"))
  err_file <- file.path(rep_dir, paste0("fit_", model_name, "_ERROR.rds"))
  log_file <- file.path(rep_dir, paste0("fit_", model_name, ".log"))

  if (file.exists(fit_file)) {
    return(readRDS(fit_file))
  }
  if (file.exists(err_file) && !isTRUE(rerun_failed)) {
    return(readRDS(err_file))
  }

  t0 <- Sys.time()
  res <- NULL
  out <- capture.output({
    res <- try(expr_fun(), silent = TRUE)
  }, type = "output")
  writeLines(out, con = log_file)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(res, "try-error")) {
    obj <- c(metadata, list(ok = FALSE, model = model_name, elapsed_sec = elapsed, error = as.character(res)))
    saveRDS_atomic(obj, err_file)
  } else {
    obj <- c(metadata, list(ok = TRUE, model = model_name, elapsed_sec = elapsed, fit = res))
    saveRDS_atomic(obj, fit_file)
  }

  if (!is.null(progress_env)) {
    progress_env$current <- progress_env$current + 1L
    if (!is.null(progress_env$pb)) utils::setTxtProgressBar(progress_env$pb, progress_env$current)
  }

  obj
}

run_one_rep_poster <- function(bg_scenario, count_regime, rep_id,
                               mu_cal, k0_cal, seed_rep,
                               cell_dir,
                               progress_env = NULL,
                               verbose = TRUE,
                               rerun_failed = FALSE) {
  rep_dir <- file.path(cell_dir, sprintf("rep_%03d", rep_id))
  .safe_dir_create(rep_dir)
  progress_file <- file.path(dirname(cell_dir), "fit_progress.csv")

  if (isTRUE(verbose)) {
    cat(sprintf("\n[REP] %s | %s | rep %03d | seed=%d\n",
                bg_scenario$bg_case, count_regime$count_case, rep_id, seed_rep))
  }

  sim_file <- file.path(rep_dir, "simulation.rds")
  sim_err_file <- file.path(rep_dir, "simulation_ERROR.rds")

  if (file.exists(sim_file)) {
    sim_obj <- readRDS(sim_file)
  } else if (file.exists(sim_err_file) && !isTRUE(rerun_failed)) {
    return(list(ok = FALSE, stage = "simulation", error = readRDS(sim_err_file)$error))
  } else {
    sim_obj <- try(
      simulate_one_from_cell_poster(
        bg_scenario = bg_scenario,
        count_regime = count_regime,
        seed = seed_rep,
        mu_value = mu_cal,
        k0_value = k0_cal,
        return_bg_info = TRUE,
        return_support = TRUE
      ),
      silent = TRUE
    )
    if (inherits(sim_obj, "try-error")) {
      saveRDS_atomic(list(ok = FALSE, error = as.character(sim_obj)), sim_err_file)
      return(list(ok = FALSE, stage = "simulation", error = as.character(sim_obj)))
    }
    saveRDS_atomic(sim_obj, sim_file)
  }

  cnt <- get_counts(sim_obj)
  if (isTRUE(cnt["exploded"])) {
    return(list(ok = FALSE, stage = "simulation", error = "exploded catalog"))
  }

  cat_fit <- make_eqcat_from_sim(sim_obj)
  starts <- list(
    mu = mu_cal * 0.85,
    k0 = k0_cal * 1.15,
    c = TRUE_BASE["c"] * 1.5,
    p = 1.20,
    gamma = 0,
    d = TRUE_BASE["d"] * 0.8,
    q = 1.70,
    betacov = TRUE_BETACOV * 0.7
  )

  metadata <- list(
    bg_case = bg_scenario$bg_case,
    bg_label = bg_scenario$label,
    count_case = count_regime$count_case,
    rep = rep_id,
    seed = seed_rep,
    n0 = unname(cnt["n0"]),
    nson = unname(cnt["nson"]),
    nobs = nrow(cat_fit),
    mu_cal = mu_cal,
    k0_cal = k0_cal
  )

  if (isTRUE(verbose)) cat("  - fit classic ETAS...\n")
  classic_obj <- run_fit_stage(
    model_name = "classic",
    expr_fun = function() fit_etas_classic_wrapper_poster(cat_fit, starts = starts),
    rep_dir = rep_dir,
    metadata = metadata,
    progress_env = progress_env,
    rerun_failed = rerun_failed
  )
  append_csv_row(c(metadata, list(model = "classic", ok = classic_obj$ok, elapsed_sec = classic_obj$elapsed_sec %||% NA_real_)), progress_file)

  if (isTRUE(verbose)) cat("  - fit parametric-background ETAS...\n")
  param_obj <- run_fit_stage(
    model_name = "parametric",
    expr_fun = function() fit_etas_parametric_wrapper_poster(cat_fit, bg_case = bg_scenario$bg_case, starts = starts),
    rep_dir = rep_dir,
    metadata = metadata,
    progress_env = progress_env,
    rerun_failed = rerun_failed
  )
  append_csv_row(c(metadata, list(model = "parametric", ok = param_obj$ok, elapsed_sec = param_obj$elapsed_sec %||% NA_real_)), progress_file)

  metrics_file <- file.path(rep_dir, "metrics.rds")
  if (file.exists(metrics_file) && !isTRUE(rerun_failed)) {
    return(readRDS(metrics_file))
  }

  true_params_full <- c(
    mu = unname(mu_cal),
    k0 = unname(k0_cal),
    c = unname(TRUE_BASE["c"]),
    p = unname(TRUE_BASE["p"]),
    gamma = 0,
    d = unname(TRUE_BASE["d"]),
    q = unname(TRUE_BASE["q"])
  )

  true_int <- try(compute_true_event_intensity(sim_obj), silent = TRUE)

  model_metrics <- list()
  trig_rows <- list()
  int_rows <- list()
  bg_shape_rows <- list()
  bg_coef_rows <- list()

  for (obj in list(classic_obj, param_obj)) {
    model_name <- if (identical(obj$model, "classic")) "classic" else "parametric"
    if (!isTRUE(obj$ok)) {
      model_metrics[[model_name]] <- data.frame(
        c(metadata, list(
          model = model_name,
          ok_fit = FALSE,
          elapsed_sec = obj$elapsed_sec %||% NA_real_,
          rmse_trigger7 = NA_real_,
          error = obj$error %||% NA_character_
        )),
        as.list(empty_int_metrics())
      )
      next
    }

    fit_obj <- obj$fit
    est_par <- try(extract_fit_trigger_params(fit_obj), silent = TRUE)
    if (!inherits(est_par, "try-error")) {
      err_trig <- compute_trigger_errors_one(est_par, true_params_full, TRUE_BETACOV)
      rmse7 <- rmse_trigger7_one(err_trig)
      trig_rows[[model_name]] <- data.frame(
        c(metadata, list(model = model_name)),
        parameter = names(err_trig),
        estimate = as.numeric(c(est_par[c("mu", "k0", "c", "p", "d", "q")], betacov = est_par["betacov"])),
        truth = as.numeric(c(true_params_full[c("mu", "k0", "c", "p", "d", "q")], betacov = TRUE_BETACOV)),
        error = as.numeric(err_trig),
        sq_error = as.numeric(err_trig)^2
      )
    } else {
      rmse7 <- NA_real_
    }

    int_metrics <- empty_int_metrics()
    if (!inherits(true_int, "try-error")) {
      fit_int <- try(extract_fit_event_intensity(fit_obj), silent = TRUE)
      if (!inherits(fit_int, "try-error")) {
        m <- min(nrow(true_int), nrow(fit_int))
        int_metrics <- compute_intensity_metrics(true_int[seq_len(m), , drop = FALSE], fit_int[seq_len(m), , drop = FALSE])
        int_rows[[model_name]] <- data.frame(c(metadata, list(model = model_name)), as.list(int_metrics))
      }
    }

    if (model_name == "parametric") {
      bg_shape <- try(compare_background_shape_param(sim_obj, fit_obj), silent = TRUE)
      if (!inherits(bg_shape, "try-error")) {
        bg_shape_rows[[model_name]] <- data.frame(c(metadata, list(model = model_name)), as.list(bg_shape))
      }

      if (isTRUE(bg_scenario$has_bg_coef)) {
        true_cf <- bg_scenario$true_bg_coef
        est_cf <- extract_param_bg_coefs(fit_obj, names(true_cf))
        bg_coef_rows[[model_name]] <- data.frame(
          c(metadata, list(model = model_name)),
          coef = names(true_cf),
          estimate = as.numeric(est_cf),
          truth = as.numeric(true_cf),
          error = as.numeric(est_cf - true_cf),
          sq_error = as.numeric((est_cf - true_cf)^2)
        )
      }
    }

    model_metrics[[model_name]] <- data.frame(
      c(metadata, list(
        model = model_name,
        ok_fit = TRUE,
        elapsed_sec = obj$elapsed_sec,
        rmse_trigger7 = rmse7,
        error = NA_character_
      )),
      as.list(int_metrics)
    )
  }

  out <- list(
    ok = TRUE,
    metadata = metadata,
    main = if (length(model_metrics)) do.call(rbind_fill_base, model_metrics) else data.frame(),
    trig_err = if (length(trig_rows)) do.call(rbind_fill_base, trig_rows) else data.frame(),
    intensity = if (length(int_rows)) do.call(rbind_fill_base, int_rows) else data.frame(),
    bg_shape = if (length(bg_shape_rows)) do.call(rbind_fill_base, bg_shape_rows) else data.frame(),
    bg_coef = if (length(bg_coef_rows)) do.call(rbind_fill_base, bg_coef_rows) else data.frame()
  )

  saveRDS_atomic(out, metrics_file)
  out
}

## ================================================================
## 7) ESECUZIONE DEL PIANO COMPLETO POSTER
## ================================================================

.count_completed_fit_attempts <- function(checkpoint_dir) {
  if (!dir.exists(checkpoint_dir)) return(0L)
  length(list.files(checkpoint_dir, pattern = "fit_(classic|parametric)(|_ERROR)\\.rds$", recursive = TRUE, full.names = TRUE))
}

run_poster_plan <- function(background_scenarios = BACKGROUND_SCENARIOS_POSTER,
                            count_regimes = COUNT_REGIMES_POSTER,
                            nrep_main = N_REP_MAIN_POSTER,
                            seed_master = SEED_MASTER_POSTER,
                            checkpoint_dir = "chk_poster_etas",
                            resume = TRUE,
                            rerun_failed = FALSE,
                            rerun_calibration = FALSE,
                            verbose = TRUE) {
  .safe_dir_create(checkpoint_dir)
  global_file <- file.path(checkpoint_dir, "poster_plan_checkpoint.rds")

  cells <- expand.grid(
    bg_case = names(background_scenarios),
    count_case = count_regimes$count_case,
    stringsAsFactors = FALSE
  )

  total_fits <- nrow(cells) * nrep_main * 2L
  completed0 <- if (isTRUE(resume)) .count_completed_fit_attempts(checkpoint_dir) else 0L
  completed0 <- min(completed0, total_fits)
  pb <- utils::txtProgressBar(min = 0, max = total_fits, initial = completed0, style = 3)
  progress_env <- new.env(parent = emptyenv())
  progress_env$current <- completed0
  progress_env$pb <- pb
  on.exit(close(pb), add = TRUE)

  cell_results <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    bg_case_i <- cells$bg_case[i]
    count_case_i <- cells$count_case[i]
    bg_scen <- background_scenarios[[bg_case_i]]
    cnt_reg <- count_regimes[count_regimes$count_case == count_case_i, , drop = FALSE][1, ]

    cell_name <- paste(bg_case_i, count_case_i, sep = "__")
    cell_dir <- file.path(checkpoint_dir, "cells", cell_name)
    .safe_dir_create(cell_dir)
    calib_file <- file.path(cell_dir, "calibration.rds")
    cell_file <- file.path(cell_dir, "cell_result.rds")

    if (isTRUE(verbose)) {
      cat("\n\n============================================================\n")
      cat(sprintf("[CELL %d/%d] %s | %s\n", i, nrow(cells), bg_scen$label, count_case_i))
      cat("============================================================\n")
    }

    if (isTRUE(resume) && file.exists(calib_file) && !isTRUE(rerun_calibration)) {
      calib <- readRDS(calib_file)
      if (isTRUE(verbose)) cat(sprintf("[RESUME CALIB] mu=%.5f k0=%.6f\n", calib$mu_cal, calib$k0_cal))
    } else {
      calib <- calibrate_cell_poster(
        bg_scenario = bg_scen,
        count_regime = cnt_reg,
        seed_base = seed_master + i * 100000L,
        verbose = verbose
      )
      saveRDS_atomic(calib, calib_file)
    }

    reps <- vector("list", nrep_main)
    for (r in seq_len(nrep_main)) {
      seed_rep <- seed_master + i * 100000L + 10000L + r
      reps[[r]] <- run_one_rep_poster(
        bg_scenario = bg_scen,
        count_regime = cnt_reg,
        rep_id = r,
        mu_cal = calib$mu_cal,
        k0_cal = calib$k0_cal,
        seed_rep = seed_rep,
        cell_dir = cell_dir,
        progress_env = progress_env,
        verbose = verbose,
        rerun_failed = rerun_failed
      )
    }

    cell_res <- list(
      bg_case = bg_case_i,
      bg_label = bg_scen$label,
      count_case = count_case_i,
      count_regime = cnt_reg,
      mu_cal = calib$mu_cal,
      k0_cal = calib$k0_cal,
      pilot_log = calib$pilot_log,
      reps = reps,
      finished_at = Sys.time()
    )
    saveRDS_atomic(cell_res, cell_file)
    cell_results[[i]] <- cell_res

    saveRDS_atomic(
      list(cells = cells, cell_results = cell_results, checkpoint_dir = checkpoint_dir, updated_at = Sys.time()),
      global_file
    )
  }

  res <- collect_poster_results(checkpoint_dir = checkpoint_dir)
  saveRDS_atomic(res, file.path(checkpoint_dir, "poster_plan_final_result.rds"))
  res
}

collect_poster_results <- function(checkpoint_dir = "chk_poster_etas") {
  metric_files <- list.files(checkpoint_dir, pattern = "metrics\\.rds$", recursive = TRUE, full.names = TRUE)
  metrics <- lapply(metric_files, readRDS)

  bind_piece <- function(name) {
    pieces <- lapply(metrics, `[[`, name)
    pieces <- pieces[vapply(pieces, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
    if (!length(pieces)) data.frame() else do.call(rbind_fill_base, pieces)
  }

  calib_files <- list.files(checkpoint_dir, pattern = "calibration\\.rds$", recursive = TRUE, full.names = TRUE)
  pilot_logs <- lapply(calib_files, function(f) readRDS(f)$pilot_log)
  pilot_logs <- pilot_logs[vapply(pilot_logs, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]

  list(
    tables = list(
      main = bind_piece("main"),
      trig_err = bind_piece("trig_err"),
      intensity = bind_piece("intensity"),
      bg_shape = bind_piece("bg_shape"),
      bg_coef = bind_piece("bg_coef"),
      pilot_log = if (length(pilot_logs)) do.call(rbind_fill_base, pilot_logs) else data.frame()
    ),
    metadata = list(
      checkpoint_dir = checkpoint_dir,
      collected_at = Sys.time(),
      n_metric_files = length(metric_files),
      background_scenarios = names(BACKGROUND_SCENARIOS_POSTER),
      count_regimes = COUNT_REGIMES_POSTER$count_case,
      nrep_main = N_REP_MAIN_POSTER
    )
  )
}

summarise_poster_results <- function(plan_res) {
  tb_main <- plan_res$tables$main %||% data.frame()
  tb_trig <- plan_res$tables$trig_err %||% data.frame()
  tb_int <- plan_res$tables$intensity %||% data.frame()
  tb_bgshape <- plan_res$tables$bg_shape %||% data.frame()
  tb_bgcoef <- plan_res$tables$bg_coef %||% data.frame()

  agg_mean <- function(formula, data) {
    if (!is.data.frame(data) || !nrow(data)) return(data.frame())
    stats::aggregate(formula, data = data, FUN = function(z) mean(z, na.rm = TRUE))
  }

  count_summary <- if (nrow(tb_main)) {
    tmp_counts <- tb_main[!duplicated(tb_main[c("bg_case", "bg_label", "count_case", "rep")]),
                          c("bg_case", "bg_label", "count_case", "rep", "n0", "nson", "nobs"), drop = FALSE]
    out_counts <- stats::aggregate(cbind(n0, nson, nobs) ~ bg_label + count_case,
                                   data = tmp_counts,
                                   FUN = function(z) mean(z, na.rm = TRUE))
    out_counts <- merge(out_counts,
                        COUNT_REGIMES_POSTER[, c("count_case", "N_bg_target", "N_trig_target", "N_tot_target")],
                        by = "count_case", all.x = TRUE, sort = FALSE)
    out_counts$err_n0 <- out_counts$n0 - out_counts$N_bg_target
    out_counts$err_nson <- out_counts$nson - out_counts$N_trig_target
    out_counts$trig_fraction <- out_counts$nson / pmax(out_counts$nobs, 1)
    out_counts$target_trig_fraction <- out_counts$N_trig_target / out_counts$N_tot_target
    out_counts
  } else data.frame()

  list(
    count_summary = count_summary,
    fit_success = if (nrow(tb_main)) stats::aggregate(ok_fit ~ bg_label + count_case + model, data = tb_main, FUN = mean) else data.frame(),
    rmse_trigger7 = if (nrow(tb_main)) agg_mean(rmse_trigger7 ~ bg_label + count_case + model, tb_main) else data.frame(),
    trig_param_rmse = if (nrow(tb_trig)) {
      tmp <- tb_trig
      tmp$rmse <- sqrt(tmp$sq_error)
      agg_mean(rmse ~ bg_label + count_case + model + parameter, tmp)
    } else data.frame(),
    intensity = if (nrow(tb_int)) {
      cols <- intersect(c("rmse_log_lambda_tot", "rmse_log_lambda_bg", "cor_log_lambda_bg"), names(tb_int))
      pieces <- lapply(cols, function(cc) {
        out <- stats::aggregate(tb_int[[cc]], by = tb_int[c("bg_label", "count_case", "model")], FUN = function(z) mean(z, na.rm = TRUE))
        names(out)[ncol(out)] <- "mean"
        out$metric <- cc
        out
      })
      if (length(pieces)) do.call(rbind_fill_base, pieces) else data.frame()
    } else data.frame(),
    bg_shape = if (nrow(tb_bgshape)) agg_mean(cbind(rmse_log_shape, corr_log_shape, rmse_dens) ~ bg_label + count_case + model, tb_bgshape) else data.frame(),
    bg_coef = if (nrow(tb_bgcoef)) {
      tmp <- tb_bgcoef
      tmp$rmse <- sqrt(tmp$sq_error)
      stats::aggregate(cbind(estimate, error, rmse) ~ bg_label + count_case + coef, data = tmp, FUN = function(z) mean(z, na.rm = TRUE))
    } else data.frame()
  )
}


## ================================================================
## TEST RAPIDO DELLA NUOVA CALIBRAZIONE + SIMULAZIONE PATTERN
## ================================================================

test_calibration_patterns <- function(
    background_scenarios = BACKGROUND_SCENARIOS_POSTER,
    count_regimes = COUNT_REGIMES_POSTER,
    nrep_pilot = 8,          # per test rapido; per il piano finale puoi usare 20
    n_iters = 5,             # per test rapido; per il piano finale puoi usare 8
    n_bracket = 6,
    n_check = 10,            # numero di pattern simulati dopo calibrazione
    checkpoint_dir = "chk_calibration_only_TEST",
    seed_base = 987000L,
    plot_first_pattern = TRUE,
    verbose = TRUE
) {
  
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  calib_list <- list()
  calib_log_list <- list()
  count_rows <- list()
  
  n_cells <- length(background_scenarios) * nrow(count_regimes)
  total_steps <- n_cells * (1 + n_check)
  pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0L
  
  if (isTRUE(plot_first_pattern)) {
    pdf_file <- file.path(checkpoint_dir, "first_pattern_per_cell.pdf")
    grDevices::pdf(pdf_file, width = 12, height = 8)
    op <- par(no.readonly = TRUE)
    on.exit({
      par(op)
      grDevices::dev.off()
    }, add = TRUE)
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  }
  
  cell_id <- 0L
  
  for (bg_name in names(background_scenarios)) {
    
    bg_scenario <- background_scenarios[[bg_name]]
    
    for (j in seq_len(nrow(count_regimes))) {
      
      count_regime <- count_regimes[j, , drop = FALSE]
      cell_id <- cell_id + 1L
      
      if (isTRUE(verbose)) {
        cat("\n\n============================================================\n")
        cat(sprintf(
          "[CALIBRATION TEST CELL %d/%d] %s | %s\n",
          cell_id, n_cells, bg_scenario$label, count_regime$count_case
        ))
        cat("============================================================\n")
      }
      
      ## ------------------------------------------------------------
      ## 1) Calibrazione solo mu/k0
      ## ------------------------------------------------------------
      
      calib <- calibrate_cell_poster(
        bg_scenario = bg_scenario,
        count_regime = count_regime,
        nrep_pilot = nrep_pilot,
        n_iters = n_iters,
        n_bracket = n_bracket,
        seed_base = seed_base + 10000L * cell_id,
        verbose = verbose
      )
      
      calib_key <- paste(bg_name, count_regime$count_case, sep = "__")
      calib_list[[calib_key]] <- calib
      calib_log_list[[calib_key]] <- calib$pilot_log
      
      saveRDS(
        calib,
        file = file.path(checkpoint_dir, paste0("calibration_", calib_key, ".rds"))
      )
      
      step <- step + 1L
      utils::setTxtProgressBar(pb, step)
      
      ## ------------------------------------------------------------
      ## 2) Simula n_check pattern con i parametri calibrati
      ## ------------------------------------------------------------
      
      for (r in seq_len(n_check)) {
        
        seed_r <- seed_base + 100000L * cell_id + r
        
        sim_r <- simulate_one_from_cell_poster(
          bg_scenario = bg_scenario,
          count_regime = count_regime,
          seed = seed_r,
          mu_value = calib$mu_cal,
          k0_value = calib$k0_cal,
          return_bg_info = TRUE,
          return_support = r == 1L,
          plot_bg = FALSE,
          plot_catalog = FALSE
        )
        
        cc <- get_counts(sim_r)
        cat_r <- get_cat_sim(sim_r)
        
        n0_r <- as.numeric(cc["n0"])
        nson_r <- as.numeric(cc["nson"])
        nobs_r <- n0_r + nson_r
        
        count_rows[[length(count_rows) + 1L]] <- data.frame(
          bg_case = bg_scenario$bg_case,
          bg_label = bg_scenario$label,
          count_case = count_regime$count_case,
          rep = r,
          seed = seed_r,
          mu_cal = calib$mu_cal,
          k0_cal = calib$k0_cal,
          n0 = n0_r,
          nson = nson_r,
          nobs = nobs_r,
          frac_triggered = ifelse(nobs_r > 0, nson_r / nobs_r, NA_real_),
          target_bg = count_regime$N_bg_target,
          target_trig = count_regime$N_trig_target,
          target_total = count_regime$N_tot_target,
          target_frac_triggered = count_regime$N_trig_target / count_regime$N_tot_target,
          exploded = isTRUE(cc["exploded"]),
          stringsAsFactors = FALSE
        )
        
        ## Plot del primo pattern della cella
        if (isTRUE(plot_first_pattern) && r == 1L) {
          is_bg <- cat_r$father_id == 0
          
          plot(
            cat_r$long,
            cat_r$lat,
            pch = 16,
            cex = 0.45,
            col = ifelse(is_bg, "blue", "red"),
            xlab = "Longitude",
            ylab = "Latitude",
            main = paste0(bg_scenario$bg_case, " | ", count_regime$count_case)
          )
          
          legend(
            "topright",
            legend = c(
              paste0("background: ", sum(is_bg)),
              paste0("triggered: ", sum(!is_bg))
            ),
            col = c("blue", "red"),
            pch = 16,
            bty = "n",
            cex = 0.85
          )
        }
        
        step <- step + 1L
        utils::setTxtProgressBar(pb, step)
      }
    }
  }
  
  close(pb)
  
  counts_df <- do.call(rbind, count_rows)
  calib_log_df <- do.call(rbind, calib_log_list)
  
  ## ------------------------------------------------------------
  ## 3) Riassunto conteggi simulati
  ## ------------------------------------------------------------
  
  split_counts <- split(
    counts_df,
    interaction(counts_df$bg_case, counts_df$count_case, drop = TRUE)
  )
  
  summary_counts <- do.call(
    rbind,
    lapply(split_counts, function(d) {
      data.frame(
        bg_case = d$bg_case[1],
        bg_label = d$bg_label[1],
        count_case = d$count_case[1],
        n_patterns = nrow(d),
        
        target_bg = d$target_bg[1],
        target_trig = d$target_trig[1],
        target_total = d$target_total[1],
        target_frac_triggered = d$target_frac_triggered[1],
        
        mean_n0 = mean(d$n0, na.rm = TRUE),
        sd_n0 = stats::sd(d$n0, na.rm = TRUE),
        rel_error_n0 = (mean(d$n0, na.rm = TRUE) - d$target_bg[1]) / d$target_bg[1],
        
        mean_nson = mean(d$nson, na.rm = TRUE),
        sd_nson = stats::sd(d$nson, na.rm = TRUE),
        rel_error_nson = (mean(d$nson, na.rm = TRUE) - d$target_trig[1]) / d$target_trig[1],
        
        mean_nobs = mean(d$nobs, na.rm = TRUE),
        sd_nobs = stats::sd(d$nobs, na.rm = TRUE),
        rel_error_nobs = (mean(d$nobs, na.rm = TRUE) - d$target_total[1]) / d$target_total[1],
        
        mean_frac_triggered = mean(d$frac_triggered, na.rm = TRUE),
        sd_frac_triggered = stats::sd(d$frac_triggered, na.rm = TRUE),
        
        exploded_rate = mean(d$exploded, na.rm = TRUE),
        mu_cal = d$mu_cal[1],
        k0_cal = d$k0_cal[1],
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(summary_counts) <- NULL
  
  out <- list(
    calibrations = calib_list,
    calibration_log = calib_log_df,
    counts = counts_df,
    summary_counts = summary_counts
  )
  
  saveRDS(
    out,
    file = file.path(checkpoint_dir, "calibration_pattern_test_results.rds")
  )
  
  utils::write.csv(
    counts_df,
    file = file.path(checkpoint_dir, "calibration_pattern_test_counts.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    summary_counts,
    file = file.path(checkpoint_dir, "calibration_pattern_test_summary.csv"),
    row.names = FALSE
  )
  
  if (isTRUE(verbose)) {
    cat("\n\n==================== SUMMARY COUNTS ====================\n")
    print(summary_counts)
    cat("\nRisultati salvati in:\n")
    cat(checkpoint_dir, "\n")
    if (isTRUE(plot_first_pattern)) {
      cat("PDF pattern:", file.path(checkpoint_dir, "first_pattern_per_cell.pdf"), "\n")
    }
  }
  
  invisible(out)
}



#######################################################################
# GLOBAL SETTINGS FOR THE FINAL SIMULATION PLAN
#######################################################################

BASE_DIR <- "/home/nicolettadangelo/sim_paper1_marco/etas parametrico"
OUT_ROOT <- file.path(BASE_DIR, "etas_parametric_final_no_smooth_outputs_v1")

SAVE_FULL_FIT <- TRUE
GRID_N <- 100L

## Replications
NREP_MAIN      <- 50L
NREP_SAMPLE    <- 30L
NREP_MISSPEC   <- 30L
NREP_EXTREME   <- 30L

## Calibration settings
NREP_CALIB      <- 20L
N_CALIB_ITERS   <- 8L
N_CALIB_BRACKET <- 8L

## Fit settings inherited from the poster plan, but made explicit here.
NDECLUST_FIT_NEW <- 4L
ITERLIM_FIT_NEW  <- 40L
MULT_BG_FIT_NEW  <- 4L

SEED_MASTER_FULL <- 20260514L

## Default autorun mode.
RUN_MODE <- Sys.getenv("ETAS_RUN_MODE", unset = "none")

## Required helper for marked background PPP construction inside etasclass.par.
## The etasclass.par definition above has been patched to call this function.
make_bg_ppp_with_optional_marks <- function(cat, win.bg, formula.bg) {
  vars <- all.vars(formula.bg)
  non_mark_vars <- c(
    "x", "y", "z", "t",
    "x_std", "y_std", "Z1", "r2_center",
    "time", "lat", "long", "magn1", "magnitude", "time.work",
    "xcat.work", "ycat.work", "ind", "ord", "hvarx", "hvary"
  )
  mark_vars <- intersect(setdiff(vars, non_mark_vars), names(cat))

  if (length(mark_vars)) {
    mk <- as.data.frame(cat[, mark_vars, drop = FALSE])
    for (j in seq_along(mk)) mk[[j]] <- as.factor(mk[[j]])
    spatstat.geom::ppp(
      x = cat$xcat.work,
      y = cat$ycat.work,
      window = win.bg,
      marks = mk,
      check = FALSE
    )
  } else {
    spatstat.geom::ppp(
      x = cat$xcat.work,
      y = cat$ycat.work,
      window = win.bg,
      check = FALSE
    )
  }
}

## Optional plotting package used only by downstream plots, not by core fitting.
suppressPackageStartupMessages({
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})


#######################################################################
# OPTIONAL POSTER PLOTTING HELPER COPIED FROM PREVIOUS FILE
#######################################################################

make_bg_coef_plot <- function(data,
                              coef_name,
                              y_lab,
                              legend_lab,
                              title_lab = NULL,
                              color_limits = NULL) {
  
  df <- data %>%
    filter(coef == coef_name) %>%
    mutate(
      error_coef = estimate - truth
    ) %>%
    group_by(regime_plot) %>%
    mutate(
      median_bias = median(error_coef, na.rm = TRUE),
      abs_median_bias = abs(median_bias)
    ) %>%
    ungroup()
  
  truth_df <- df %>%
    distinct(truth)
  
  if (is.null(color_limits)) {
    color_limits <- c(0, max(df$abs_median_bias, na.rm = TRUE))
  }
  
  y_range <- range(c(df$estimate, df$truth), na.rm = TRUE)
  y_pad <- diff(y_range) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1
  
  ggplot(
    df,
    aes(
      x = regime_plot,
      y = estimate,
      fill = abs_median_bias
    )
  ) +
    geom_hline(
      data = truth_df,
      aes(yintercept = truth),
      linetype = "dashed",
      linewidth = 0.9,
      color = "grey20",
      inherit.aes = FALSE
    ) +
    geom_boxplot(
      width = 0.48,
      linewidth = 0.85,
      outlier.size = 1.5,
      outlier.alpha = 0.75,
      color = "grey15"
    ) +
    scale_x_discrete(
      labels = function(x) parse(text = x)
    ) +
    scale_fill_gradient(
      low = "#c6dbef",
      high = "#b2182b",
      limits = color_limits,
      name = legend_lab
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(5.2, "cm"),
        barwidth  = grid::unit(1.25, "cm"),
        ticks = TRUE
      )
    ) +
    coord_cartesian(
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)
    ) +
    labs(
      x = "",
      y = y_lab
    ) +
    theme_bw(base_size = 26) +
    theme(
      plot.title = element_text(
        size = 32,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      
      axis.title.y = element_text(
        size = 38,
        face = "bold",
        margin = margin(r = 14)
      ),
      axis.text.y = element_text(
        size = 31,
        color = "grey20"
      ),
      axis.text.x = element_text(
        size = 30,
        angle = 0,
        color = "grey20"
      ),
      
      legend.position = "right",
      legend.title = element_text(
        size = 28,
        face = "bold"
      ),
      legend.text = element_text(
        size = 24
      ),
      legend.key.height = grid::unit(1.3, "cm"),
      legend.key.width  = grid::unit(1.2, "cm"),
      
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "grey25",
        linewidth = 0.85
      ),
      
      plot.margin = margin(12, 24, 12, 18)
    )
}


#######################################################################
# FINAL FULL SIMULATION CONTROLLER
#######################################################################

## 2. SMALL UTILITIES
## ================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

saveRDS_atomic <- function(object, file, compress = "gzip") {
  safe_dir_create(dirname(file))
  tmp <- paste0(file, ".tmp_", Sys.getpid())
  saveRDS(object, tmp, compress = compress)
  ok <- file.rename(tmp, file)
  if (!ok) {
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(file)
}

append_csv_row <- function(row, file) {
  safe_dir_create(dirname(file))
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

log_msg <- function(..., .level = "INFO") {
  cat(sprintf("[%s] [%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), .level),
      sprintf(...), "\n", sep = "")
  flush.console()
}

safe_run <- function(expr_fun) {
  warnings <- character(0)
  t0 <- Sys.time()
  res <- withCallingHandlers(
    tryCatch(
      expr_fun(),
      error = function(e) structure(list(error = conditionMessage(e)), class = "safe_error")
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  t1 <- Sys.time()
  if (inherits(res, "safe_error")) {
    list(ok = FALSE, value = NULL, error_message = res$error, warnings = warnings,
         elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")))
  } else {
    list(ok = TRUE, value = res, error_message = NA_character_, warnings = warnings,
         elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")))
  }
}

make_id <- function(...) {
  x <- paste(..., sep = "__")
  x <- gsub("[^A-Za-z0-9_=-]+", "_", x)
  x
}

## ================================================================
## 3. TRUE BACKGROUND DEFINITIONS
## ================================================================

## These constants are inherited from the base file unless overwritten here.
TRUE_BASE_FULL <- TRUE_BASE
TRUE_BETACOV_FULL <- TRUE_BETACOV
M0_FULL <- M0_TRUE
B_FULL  <- B_TRUE

## True background coefficients.
TRUE_BETA_XSTD_FULL <- TRUE_BETA_XSTD
TRUE_BETA_YSTD_FULL <- TRUE_BETA_YSTD
TRUE_BETA_Z1_FULL   <- TRUE_BETA_Z1
TRUE_BETA_MARK_FULL <- 0.80

## Smooth and hotspot definitions.
bg_fun_smooth_truth_full <- function(df_feat) {
  1.1 * sin(pi * df_feat$x_std) -
    0.9 * cos(pi * df_feat$y_std) +
    0.6 * (df_feat$x_std * df_feat$y_std)
}

bg_fun_hotspot_truth_full <- function(df_feat) {
  ## Two localized peaks in standardized coordinates.
  a1 <- 1.5; x1 <- -0.6; y1 <-  0.5; sig1 <- 0.35
  a2 <- 1.2; x2 <-  0.5; y2 <- -0.4; sig2 <- 0.25
  a1 * exp(-((df_feat$x_std - x1)^2 + (df_feat$y_std - y1)^2) / (2 * sig1^2)) +
    a2 * exp(-((df_feat$x_std - x2)^2 + (df_feat$y_std - y2)^2) / (2 * sig2^2))
}

bg_cov_Z1_truth_full <- function(df_feat) {
  sin(pi * df_feat$x_std) * cos(pi * df_feat$y_std)
}

exact_covariate_functions_full <- function(names_needed = c("x_std", "y_std", "Z1"), window_km = .window_km()) {
  f_xstd <- function(df) background_features_from_xy(df$x, df$y, window_km)$x_std
  f_ystd <- function(df) background_features_from_xy(df$x, df$y, window_km)$y_std
  f_Z1 <- function(df) {
    ff <- background_features_from_xy(df$x, df$y, window_km)
    bg_cov_Z1_truth_full(ff)
  }
  out <- list(x_std = f_xstd, y_std = f_ystd, Z1 = f_Z1)
  out[intersect(names_needed, names(out))]
}

BACKGROUND_SCENARIOS_FULL <- list(
  constant = list(
    bg_case = "constant",
    label = "Constant background",
    sim_spec = list(bg_lp_type = "constant"),
    has_bg_coef = FALSE,
    true_bg_coef = numeric(0),
    correct_formula_key = "constant"
  ),
  linear_xy = list(
    bg_case = "linear_xy",
    label = "Linear background: x + y",
    sim_spec = list(
      bg_lp_type  = "linear",
      bg_lp_coefs = c(x_std = TRUE_BETA_XSTD_FULL, y_std = TRUE_BETA_YSTD_FULL)
    ),
    has_bg_coef = TRUE,
    true_bg_coef = c(x_std = TRUE_BETA_XSTD_FULL, y_std = TRUE_BETA_YSTD_FULL),
    correct_formula_key = "linear_xy"
  ),
  cov1 = list(
    bg_case = "cov1",
    label = "Covariate-driven background: Z1",
    sim_spec = list(
      bg_lp_type = "linear",
      bg_lp_coefs = c(Z1 = TRUE_BETA_Z1_FULL),
      bg_spatial_cov_funs = list(Z1 = bg_cov_Z1_truth_full)
    ),
    has_bg_coef = TRUE,
    true_bg_coef = c(Z1 = TRUE_BETA_Z1_FULL),
    correct_formula_key = "cov1"
  ),
  mark_cat = list(
    bg_case = "mark_cat",
    label = "Categorical mark background: two levels",
    ## Spatial simulation is uniform; the mark distribution of background events
    ## is altered after simulation according to beta_mark. This represents a
    ## background intensity on W x {A,B} with a categorical mark effect.
    sim_spec = list(bg_lp_type = "constant"),
    has_bg_coef = TRUE,
    true_bg_coef = c(mark_catB = TRUE_BETA_MARK_FULL),
    correct_formula_key = "mark_cat",
    mark_beta = TRUE_BETA_MARK_FULL,
    mark_levels = c("A", "B")
  ),
  smooth = list(
    bg_case = "smooth",
    label = "Smooth background: f(x,y)",
    sim_spec = list(
      bg_lp_type = "custom",
      bg_lp_fun = bg_fun_smooth_truth_full
    ),
    has_bg_coef = FALSE,
    true_bg_coef = numeric(0),
    correct_formula_key = "smooth"
  ),
  hotspot = list(
    bg_case = "hotspot",
    label = "Hotspot background",
    sim_spec = list(
      bg_lp_type = "custom",
      bg_lp_fun = bg_fun_hotspot_truth_full
    ),
    has_bg_coef = FALSE,
    true_bg_coef = numeric(0),
    correct_formula_key = "hotspot"
  )
)

## Final scenarios used in the paper-level simulation plan.
## The smooth background case is intentionally excluded because pilot/debug runs
## showed numerical degeneracy in the ETAS fit when an unconstrained smooth
## background is estimated jointly with all triggering parameters.
FINAL_BG_CASES <- c("constant", "linear_xy", "cov1", "mark_cat", "hotspot")
BACKGROUND_SCENARIOS_FULL <- BACKGROUND_SCENARIOS_FULL[FINAL_BG_CASES]

## ================================================================
## 4. COUNT REGIMES AND FINAL PLANS
## ================================================================

make_count_regimes_full <- function(N_total = 1000L, T_lag = T_LAG_USE,
                                    regimes = c("bgdom", "balanced", "trigdom"),
                                    c_true = TRUE_BASE_FULL["c"],
                                    p_true = TRUE_BASE_FULL["p"],
                                    d_true = TRUE_BASE_FULL["d"],
                                    q_true = TRUE_BASE_FULL["q"],
                                    beta_mag_true = TRUE_BETACOV_FULL,
                                    b_true = B_FULL) {
  defs <- list(
    bgdom    = list(count_case = paste0("N", N_total, "_bgdom_75_25"),    N_tot_target = N_total, r_trig_bg = 1/3),
    balanced = list(count_case = paste0("N", N_total, "_balanced_50_50"), N_tot_target = N_total, r_trig_bg = 1),
    trigdom  = list(count_case = paste0("N", N_total, "_trigdom_25_75"),  N_tot_target = N_total, r_trig_bg = 3),
    extreme  = list(count_case = paste0("N", N_total, "_extreme_10_90"),  N_tot_target = N_total, r_trig_bg = 9)
  )
  x <- defs[regimes]
  grid <- do.call(rbind, lapply(x, as.data.frame, stringsAsFactors = FALSE))
  rownames(grid) <- NULL

  grid$N_bg_target <- with(grid, N_tot_target / (1 + r_trig_bg))
  grid$N_trig_target <- with(grid, N_tot_target - N_bg_target)
  grid$mu_init <- grid$N_bg_target / T_lag

  beta_GR <- log(10) * b_true
  if (beta_mag_true >= beta_GR) {
    stop("TRUE_BETACOV must be < log(10)*B_TRUE for E[exp(beta*M_rel)] to be finite.")
  }
  EexpM <- beta_GR / (beta_GR - beta_mag_true)

  C_k0 <- (c_true^(1 - p_true) / (p_true - 1)) *
    (pi * d_true^(1 - q_true) / (q_true - 1)) * EexpM

  grid$nbar_target <- with(grid, r_trig_bg / (1 + r_trig_bg))
  grid$k0_init <- grid$nbar_target / C_k0
  grid$T_lag <- T_lag
  grid
}

COUNT_MAIN_1000 <- make_count_regimes_full(1000L, regimes = c("bgdom", "balanced", "trigdom"))
COUNT_EXTREME_1000 <- make_count_regimes_full(1000L, regimes = c("extreme"))
COUNT_SAMPLE_500 <- make_count_regimes_full(500L, regimes = c("balanced", "trigdom"))
COUNT_SAMPLE_2000 <- make_count_regimes_full(2000L, regimes = c("balanced", "trigdom"))

MAIN_PLAN <- expand.grid(
  bg_case = FINAL_BG_CASES,
  count_case = COUNT_MAIN_1000$count_case,
  rep = seq_len(NREP_MAIN),
  stringsAsFactors = FALSE
)

SAMPLE_BG_CASES <- c("constant", "cov1", "mark_cat", "hotspot")
SAMPLE_COUNTS <- rbind(COUNT_SAMPLE_500, COUNT_SAMPLE_2000)
SAMPLE_PLAN <- expand.grid(
  bg_case = SAMPLE_BG_CASES,
  count_case = SAMPLE_COUNTS$count_case,
  rep = seq_len(NREP_SAMPLE),
  stringsAsFactors = FALSE
)

MISSPEC_BG_CASES <- c("linear_xy", "cov1", "mark_cat", "hotspot")
MISSPEC_COUNTS <- COUNT_MAIN_1000[COUNT_MAIN_1000$count_case %in% c("N1000_balanced_50_50", "N1000_trigdom_25_75"), , drop = FALSE]
MISSPEC_PLAN <- expand.grid(
  bg_case = MISSPEC_BG_CASES,
  count_case = MISSPEC_COUNTS$count_case,
  rep = seq_len(NREP_MISSPEC),
  stringsAsFactors = FALSE
)

EXTREME_BG_CASES <- c("constant", "cov1", "mark_cat", "hotspot")
EXTREME_PLAN <- expand.grid(
  bg_case = EXTREME_BG_CASES,
  count_case = COUNT_EXTREME_1000$count_case,
  rep = seq_len(NREP_EXTREME),
  stringsAsFactors = FALSE
)

## Pilot plan removed from the final script: all single-case checks were performed during debugging.
PILOT_PLAN <- data.frame(
  bg_case = character(0),
  count_case = character(0),
  rep = integer(0),
  stringsAsFactors = FALSE
)

count_lookup <- function(count_case) {
  all_counts <- rbind(COUNT_MAIN_1000, COUNT_EXTREME_1000, COUNT_SAMPLE_500, COUNT_SAMPLE_2000)
  out <- all_counts[all_counts$count_case == count_case, , drop = FALSE]
  if (!nrow(out)) stop("Unknown count_case: ", count_case)
  out[1, , drop = FALSE]
}

## ================================================================
## 5. SIMULATION WRAPPERS
## ================================================================

augment_mark_cat_sim <- function(sim_obj, beta_mark = TRUE_BETA_MARK_FULL,
                                 p_trigger_B = 0.50,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  s <- get_sim_core(sim_obj)
  catsim <- s$cat.sim
  if (!nrow(catsim)) return(sim_obj)

  is_bg <- catsim$father_id == 0
  p_bg_B <- exp(beta_mark) / (1 + exp(beta_mark))
  mark_B <- logical(nrow(catsim))
  mark_B[is_bg]  <- stats::rbinom(sum(is_bg), 1, p_bg_B) == 1
  mark_B[!is_bg] <- stats::rbinom(sum(!is_bg), 1, p_trigger_B) == 1

  catsim$mark_cat <- factor(ifelse(mark_B, "B", "A"), levels = c("A", "B"))
  ## Event-level true background LP on W x {A,B}. The normalizer is stored in truth.
  catsim$eta_bg_true <- ifelse(catsim$mark_cat == "B", beta_mark, 0)

  s$cat.sim <- catsim
  s$cat.pois <- catsim[catsim$father_id == 0, , drop = FALSE]
  s$truth$bg$mark <- list(
    variable = "mark_cat",
    levels = c("A", "B"),
    beta_B = beta_mark,
    p_bg_B = p_bg_B,
    p_trigger_B = p_trigger_B,
    normalizer_mark = 1 + exp(beta_mark),
    note = "Categorical mark effect is defined on the immigrant/background component. Triggered-event marks are generated independently with p_trigger_B."
  )
  sim_obj <- s
  sim_obj
}

simulate_one_full <- function(bg_scenario, count_regime, seed,
                              mu_value, k0_value,
                              return_bg_info = TRUE,
                              return_support = TRUE) {
  params_now <- TRUE_BASE_FULL
  params_now["mu"] <- mu_value
  params_now["k0"] <- k0_value

  sim_args <- list(
    params = params_now,
    m0 = M0_FULL,
    b = B_FULL,
    tmin = TMIN_USE,
    t.lag = count_regime$T_lag,
    long.range = LONG_RANGE,
    lat.range = LAT_RANGE,
    longlat.to.km = TRUE,
    sectoday = FALSE,
    trig_lp_type = "linear",
    trig_lp_coefs = c(m_rel = TRUE_BETACOV_FULL),
    n_support_obs = 1500,
    mult_support = 4,
    max_events = MAX_EVENTS,
    return_bg_info = return_bg_info,
    return_support = return_support,
    plot_bg = FALSE,
    plot_catalog = FALSE,
    seed = seed
  )
  sim_args <- c(sim_args, bg_scenario$sim_spec)
  sim_obj <- do.call(etas.par.sim_v4, sim_args)

  if (identical(bg_scenario$bg_case, "mark_cat")) {
    sim_obj <- augment_mark_cat_sim(
      sim_obj = sim_obj,
      beta_mark = bg_scenario$mark_beta,
      seed = seed + 999L
    )
  }
  sim_obj
}

calibrate_cell_full <- function(bg_scenario, count_regime, seed_base,
                                nrep_pilot = NREP_CALIB,
                                n_iters = N_CALIB_ITERS,
                                n_bracket = N_CALIB_BRACKET,
                                verbose = TRUE) {
  calibrate_cell_poster(
    bg_scenario = bg_scenario,
    count_regime = count_regime,
    nrep_pilot = nrep_pilot,
    n_iters = n_iters,
    n_bracket = n_bracket,
    seed_base = seed_base,
    verbose = verbose,
    max_exploded_rate = 0.20
  )
}

## ================================================================
## 6. FIT WRAPPERS
## ================================================================

background_fit_spec <- function(bg_case, formula_key = c("correct", "misspecified")) {
  formula_key <- match.arg(formula_key)

  if (formula_key == "correct") {
    if (bg_case == "constant") {
      return(list(formula.bg = ~ 1, spatial.cov.bg = FALSE, marked.bg = FALSE,
                  type.cov.values.bg = NULL, covs.bg.user = NULL))
    }
    if (bg_case == "linear_xy") {
      return(list(formula.bg = ~ x_std + y_std, spatial.cov.bg = TRUE, marked.bg = FALSE,
                  type.cov.values.bg = c(x_std = "exact", y_std = "exact"),
                  covs.bg.user = exact_covariate_functions_full(c("x_std", "y_std"))))
    }
    if (bg_case == "cov1") {
      return(list(formula.bg = ~ Z1, spatial.cov.bg = TRUE, marked.bg = FALSE,
                  type.cov.values.bg = c(Z1 = "exact"),
                  covs.bg.user = exact_covariate_functions_full(c("Z1"))))
    }
    if (bg_case == "mark_cat") {
      return(list(formula.bg = ~ s(mark_cat, bs = "re"), spatial.cov.bg = FALSE, marked.bg = TRUE,
                  type.cov.values.bg = NULL, covs.bg.user = NULL))
    }
    if (bg_case == "hotspot") {
      return(list(formula.bg = ~ s(x, y, k = 30), spatial.cov.bg = FALSE, marked.bg = FALSE,
                  type.cov.values.bg = NULL, covs.bg.user = NULL))
    }
  }

  if (formula_key == "misspecified") {
    if (bg_case == "linear_xy") {
      return(list(formula.bg = ~ 1, spatial.cov.bg = FALSE, marked.bg = FALSE,
                  type.cov.values.bg = NULL, covs.bg.user = NULL))
    }
    if (bg_case == "cov1") {
      return(list(formula.bg = ~ x_std + y_std, spatial.cov.bg = TRUE, marked.bg = FALSE,
                  type.cov.values.bg = c(x_std = "exact", y_std = "exact"),
                  covs.bg.user = exact_covariate_functions_full(c("x_std", "y_std"))))
    }
    if (bg_case == "mark_cat") {
      return(list(formula.bg = ~ 1, spatial.cov.bg = FALSE, marked.bg = FALSE,
                  type.cov.values.bg = NULL, covs.bg.user = NULL))
    }
    if (bg_case == "hotspot") {
      return(list(formula.bg = ~ x_std + y_std + I(x_std * y_std), spatial.cov.bg = TRUE, marked.bg = FALSE,
                  type.cov.values.bg = c(x_std = "exact", y_std = "exact"),
                  covs.bg.user = exact_covariate_functions_full(c("x_std", "y_std"))))
    }
  }

  stop("No background fit specification for bg_case=", bg_case, ", formula_key=", formula_key)
}

fit_etas_classic_full <- function(cat_sim_df, starts = NULL) {
  if (is.null(starts)) {
    starts <- list(mu = 0.4, k0 = 0.01, c = TRUE_BASE_FULL["c"], p = TRUE_BASE_FULL["p"],
                   gamma = 0, d = TRUE_BASE_FULL["d"], q = TRUE_BASE_FULL["q"],
                   betacov = TRUE_BETACOV_FULL)
  }
  fit_etas_classic_wrapper_poster(
    cat_sim_df = cat_sim_df,
    m0_true = M0_FULL,
    starts = starts,
    ndeclust = NDECLUST_FIT_NEW,
    iterlim = ITERLIM_FIT_NEW
  )
}

fit_etas_parametric_full <- function(cat_sim_df, bg_case,
                                     formula_key = c("correct", "misspecified"),
                                     starts = NULL,
                                     verbose.bg = FALSE,
                                     seed.bg = NULL) {
  formula_key <- match.arg(formula_key)
  if (is.null(starts)) {
    starts <- list(mu = 0.4, k0 = 0.01, c = TRUE_BASE_FULL["c"], p = TRUE_BASE_FULL["p"],
                   gamma = 0, d = TRUE_BASE_FULL["d"], q = TRUE_BASE_FULL["q"],
                   betacov = TRUE_BETACOV_FULL)
  }
  spec <- background_fit_spec(bg_case, formula_key)

  if (is.null(seed.bg)) {
    seed.bg <- if (identical(bg_case, "hotspot")) SEED_MASTER_FULL else NULL
  }

  etasclass.par(
    cat.orig = cat_sim_df,
    time.update = FALSE,
    magn.threshold = M0_FULL,
    magn.threshold.back = M0_FULL,
    tmax = max(cat_sim_df$time),
    long.range = LONG_RANGE,
    lat.range = LAT_RANGE,
    mu = starts$mu,
    k0 = starts$k0,
    c = starts$c,
    p = starts$p,
    gamma = 0.0,
    d = starts$d,
    q = max(starts$q, 1.01),
    betacov = starts$betacov,
    params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    formula1 = "time ~ magnitude - 1",
    offset = 0,
    hdef = c(1, 1),
    wp = rep(1, nrow(cat_sim_df)),
    hvarx = rep(1, nrow(cat_sim_df)),
    hvary = rep(1, nrow(cat_sim_df)),
    declustering = TRUE,
    thinning = FALSE,
    flp = FALSE,
    m1 = NULL,
    ndeclust = NDECLUST_FIT_NEW,
    n.iterweight = 1,
    onlytime = FALSE,
    is.backconstant = FALSE,
    description = paste("bg_case", bg_case, "formula_key", formula_key),
    cat.back = NULL,
    back.smooth = 1,
    sectoday = FALSE,
    longlat.to.km = TRUE,
    usenlm = TRUE,
    method = "BFGS",
    compsqm = TRUE,
    epsmax = 1e-4,
    iterlim = ITERLIM_FIT_NEW,
    ntheta = 36,
    formula.bg = spec$formula.bg,
    process.type.bg = "s2d",
    verbose.bg = verbose.bg,
    marked.bg = spec$marked.bg,
    mark.c.bg = FALSE,
    type.cov.values.bg = spec$type.cov.values.bg,
    grid.bg = FALSE,
    mult.bg = MULT_BG_FIT_NEW,
    ncube.bg = NULL,
    spatial.cov.bg = spec$spatial.cov.bg,
    offset_k.bg = FALSE,
    seed.bg = seed.bg,
    covs.bg.user = spec$covs.bg.user
  )
}

## ================================================================
## 7. DIAGNOSTIC PAYLOAD AND METRIC-READY EXTRACTION
## ================================================================

extract_basic_fit_summary <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  list(
    params = fit_obj$params.MLtot %||% fit_obj$params,
    params_iter = fit_obj$params.iter %||% NULL,
    sqm = fit_obj$sqm %||% NULL,
    AIC_iter = fit_obj$AIC.iter %||% NULL,
    logl = fit_obj$logl %||% NA,
    integral = fit_obj$integral %||% NA,
    rho_weights = fit_obj$rho.weights %||% NULL,
    back_dens_events = fit_obj$back.dens %||% NULL,
    lambda_events = fit_obj$l %||% NULL,
    etas_events = if (!is.null(fit_obj$logl)) attr(fit_obj$logl, "etas.vec") else NULL,
    model_bg_summary = if (!is.null(fit_obj$model.bg) && !is.null(fit_obj$model.bg$mod_global)) {
      list(
        coef = stats::coef(fit_obj$model.bg$mod_global),
        edf = tryCatch(summary(fit_obj$model.bg$mod_global)$s.table, error = function(e) NULL),
        family = tryCatch(fit_obj$model.bg$mod_global$family$family, error = function(e) NA_character_),
        formula = tryCatch(as.character(formula(fit_obj$model.bg$mod_global)), error = function(e) NA_character_)
      )
    } else NULL
  )
}

make_grid_eval <- function(bg_case, grid_n = GRID_N) {
  win <- .window_km()
  xg <- seq(win["xmin"], win["xmax"], length.out = grid_n)
  yg <- seq(win["ymin"], win["ymax"], length.out = grid_n)
  grid <- expand.grid(x = xg, y = yg)
  feat <- background_features_from_xy(grid$x, grid$y, win)
  grid$x_std <- feat$x_std
  grid$y_std <- feat$y_std
  grid$Z1 <- bg_cov_Z1_truth_full(feat)

  if (bg_case == "constant") {
    eta <- rep(0, nrow(grid))
  } else if (bg_case == "linear_xy") {
    eta <- TRUE_BETA_XSTD_FULL * grid$x_std + TRUE_BETA_YSTD_FULL * grid$y_std
  } else if (bg_case == "cov1") {
    eta <- TRUE_BETA_Z1_FULL * grid$Z1
  } else if (bg_case == "smooth") {
    eta <- bg_fun_smooth_truth_full(feat)
  } else if (bg_case == "hotspot") {
    eta <- bg_fun_hotspot_truth_full(feat)
  } else if (bg_case == "mark_cat") {
    eta <- rep(0, nrow(grid))
  } else {
    eta <- rep(NA_real_, nrow(grid))
  }

  area <- diff(range(xg)) * diff(range(yg))
  ## More stable grid-cell weight based on true window area.
  w_cell <- ((win["xmax"] - win["xmin"]) * (win["ymax"] - win["ymin"])) / nrow(grid)
  grid$w <- w_cell
  grid$eta_bg_true_spatial <- eta
  grid$lambda_shape_true <- exp(eta)
  grid$dens_bg_true_spatial <- grid$lambda_shape_true / sum(grid$w * grid$lambda_shape_true, na.rm = TRUE)
  grid
}


auc_base <- function(y, score) {
  y <- as.numeric(y)
  score <- as.numeric(score)
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]
  score <- score[ok]
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

fit_classification_metrics <- function(fit_obj, event_truth, threshold = 0.5) {
  if (is.null(fit_obj) || is.null(fit_obj$rho.weights)) return(NULL)
  rho <- as.numeric(fit_obj$rho.weights)
  if (is.null(fit_obj$cat)) return(NULL)
  df_fit <- as.data.frame(fit_obj$cat)
  if (length(rho) != nrow(df_fit)) return(NULL)

  if ("event_id" %in% names(df_fit) && "event_id" %in% names(event_truth)) {
    m <- match(df_fit$event_id, event_truth$event_id)
    if (any(is.na(m))) return(NULL)
    truth <- event_truth[m, , drop = FALSE]
  } else {
    if (nrow(event_truth) != nrow(df_fit)) return(NULL)
    truth <- event_truth
  }

  y_bg <- as.numeric(truth$is_background_true)
  p_bg <- pmin(pmax(rho, 1e-12), 1 - 1e-12)
  pred_bg <- p_bg >= threshold
  true_bg <- y_bg == 1

  TP <- sum(pred_bg & true_bg, na.rm = TRUE)
  TN <- sum(!pred_bg & !true_bg, na.rm = TRUE)
  FP <- sum(pred_bg & !true_bg, na.rm = TRUE)
  FN <- sum(!pred_bg & true_bg, na.rm = TRUE)
  n <- length(y_bg)

  safe_div <- function(a, b) ifelse(b == 0, NA_real_, a / b)
  accuracy <- safe_div(TP + TN, n)
  sensitivity_bg <- safe_div(TP, TP + FN)
  specificity_bg <- safe_div(TN, TN + FP)
  precision_bg <- safe_div(TP, TP + FP)
  f1_bg <- safe_div(2 * precision_bg * sensitivity_bg, precision_bg + sensitivity_bg)

  n_true_bg <- sum(true_bg, na.rm = TRUE)
  n_true_trig <- sum(!true_bg, na.rm = TRUE)
  n_hat_bg_expected <- sum(p_bg, na.rm = TRUE)
  n_hat_trig_expected <- sum(1 - p_bg, na.rm = TRUE)

  ord_bg <- order(p_bg, decreasing = TRUE)
  n_hat_bg_round <- max(0L, min(length(p_bg), round(n_hat_bg_expected)))
  pred_bg_top <- rep(FALSE, length(p_bg))
  if (n_hat_bg_round > 0) pred_bg_top[ord_bg[seq_len(n_hat_bg_round)]] <- TRUE
  accuracy_top <- mean(pred_bg_top == true_bg, na.rm = TRUE)

  confusion_05 <- matrix(
    c(TP, FN, FP, TN),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(True = c("background", "triggered"), Predicted = c("background", "triggered"))
  )

  list(
    count_summary = data.frame(
      class = c("background", "triggered"),
      true_count = c(n_true_bg, n_true_trig),
      expected_count_hat = c(n_hat_bg_expected, n_hat_trig_expected),
      difference = c(n_hat_bg_expected - n_true_bg, n_hat_trig_expected - n_true_trig),
      stringsAsFactors = FALSE
    ),
    confusion_05 = confusion_05,
    metrics = data.frame(
      metric = c(
        "accuracy_05", "sensitivity_bg_05", "specificity_bg_05", "precision_bg_05", "f1_bg_05",
        "brier_bg", "logloss_bg", "auc_bg", "auc_triggered", "accuracy_top_expected"
      ),
      value = c(
        accuracy, sensitivity_bg, specificity_bg, precision_bg, f1_bg,
        mean((p_bg - y_bg)^2, na.rm = TRUE),
        -mean(y_bg * log(p_bg) + (1 - y_bg) * log(1 - p_bg), na.rm = TRUE),
        auc_base(y_bg, p_bg),
        auc_base(1 - y_bg, 1 - p_bg),
        accuracy_top
      ),
      stringsAsFactors = FALSE
    )
  )
}

fit_parameter_comparison <- function(fit_obj, true_params_vec, cat_sim_df = NULL) {
  if (is.null(fit_obj) || is.null(fit_obj$params.MLtot)) return(NULL)
  true_params <- c(
    mu = as.numeric(true_params_vec[["mu"]]),
    k0 = as.numeric(true_params_vec[["k0"]]),
    c = as.numeric(TRUE_BASE_FULL[["c"]]),
    p = as.numeric(TRUE_BASE_FULL[["p"]]),
    gamma = 0,
    d = as.numeric(TRUE_BASE_FULL[["d"]]),
    q = as.numeric(TRUE_BASE_FULL[["q"]]),
    magnitude = as.numeric(TRUE_BETACOV_FULL)
  )
  est <- fit_obj$params.MLtot[names(true_params)]
  se <- if (!is.null(fit_obj$sqm)) fit_obj$sqm[names(true_params)] else rep(NA_real_, length(true_params))
  out <- data.frame(
    parameter = names(true_params),
    true = as.numeric(true_params),
    estimate = as.numeric(est),
    std_error = as.numeric(se),
    error = as.numeric(est - true_params),
    rel_error = as.numeric((est - true_params) / true_params),
    stringsAsFactors = FALSE
  )
  out$rel_error[!is.finite(out$rel_error)] <- NA_real_
  out
}

fit_trigger_components <- function(fit_obj, true_params_vec, event_truth) {
  if (is.null(fit_obj) || is.null(fit_obj$params.MLtot)) return(NULL)
  par_hat <- fit_obj$params.MLtot
  required <- c("k0", "c", "p", "d", "q", "magnitude")
  if (!all(required %in% names(par_hat))) return(NULL)
  if (!("m_rel" %in% names(event_truth))) return(NULL)

  k0_true <- as.numeric(true_params_vec[["k0"]])
  c_true  <- as.numeric(TRUE_BASE_FULL[["c"]])
  p_true  <- as.numeric(TRUE_BASE_FULL[["p"]])
  d_true  <- as.numeric(TRUE_BASE_FULL[["d"]])
  q_true  <- as.numeric(TRUE_BASE_FULL[["q"]])
  b_true  <- as.numeric(TRUE_BETACOV_FULL)

  k0_hat <- as.numeric(par_hat[["k0"]])
  c_hat  <- as.numeric(par_hat[["c"]])
  p_hat  <- as.numeric(par_hat[["p"]])
  d_hat  <- as.numeric(par_hat[["d"]])
  q_hat  <- as.numeric(par_hat[["q"]])
  b_hat  <- as.numeric(par_hat[["magnitude"]])

  m_rel <- event_truth$m_rel
  temporal_integral <- function(c, p) if (!is.finite(c) || !is.finite(p) || c <= 0 || p <= 1) NA_real_ else c^(1 - p) / (p - 1)
  spatial_integral <- function(d, q) if (!is.finite(d) || !is.finite(q) || d <= 0 || q <= 1) NA_real_ else pi * d^(1 - q) / (q - 1)

  mean_prod_true <- mean(k0_true * exp(b_true * m_rel), na.rm = TRUE)
  mean_prod_hat  <- mean(k0_hat  * exp(b_hat  * m_rel), na.rm = TRUE)
  It_true <- temporal_integral(c_true, p_true)
  It_hat  <- temporal_integral(c_hat, p_hat)
  Is_true <- spatial_integral(d_true, q_true)
  Is_hat  <- spatial_integral(d_hat, q_hat)

  out <- data.frame(
    component = c("mean_productivity_k0_exp_beta_m", "temporal_integral", "spatial_integral", "approx_branching_component", "mean_prod_x_spatial_integral"),
    true = c(mean_prod_true, It_true, Is_true, mean_prod_true * It_true * Is_true, mean_prod_true * Is_true),
    estimated = c(mean_prod_hat, It_hat, Is_hat, mean_prod_hat * It_hat * Is_hat, mean_prod_hat * Is_hat),
    stringsAsFactors = FALSE
  )
  out$error <- out$estimated - out$true
  out$rel_error <- out$error / out$true
  out$rel_error[!is.finite(out$rel_error)] <- NA_real_
  out
}

fit_background_metrics <- function(fit_obj, bg_case, grid_truth = NULL) {
  if (is.null(fit_obj) || is.null(fit_obj$model.bg) || is.null(fit_obj$model.bg$mod_global)) return(NULL)
  mod <- fit_obj$model.bg$mod_global
  coef_bg <- tryCatch(stats::coef(mod), error = function(e) NULL)
  out <- list(coef = coef_bg)

  if (bg_case == "linear_xy") {
    out$coef_comparison <- data.frame(
      coefficient = c("x_std", "y_std"),
      true = c(TRUE_BETA_XSTD_FULL, TRUE_BETA_YSTD_FULL),
      estimate = as.numeric(coef_bg[c("x_std", "y_std")]),
      stringsAsFactors = FALSE
    )
    out$coef_comparison$error <- out$coef_comparison$estimate - out$coef_comparison$true
  }

  if (bg_case == "cov1") {
    out$coef_comparison <- data.frame(
      coefficient = "Z1",
      true = TRUE_BETA_Z1_FULL,
      estimate = as.numeric(coef_bg["Z1"]),
      stringsAsFactors = FALSE
    )
    out$coef_comparison$error <- out$coef_comparison$estimate - out$coef_comparison$true
  }

  if (bg_case == "mark_cat") {
    re_names <- grep("^s\\(mark_cat\\)", names(coef_bg), value = TRUE)
    if (length(re_names) >= 2) {
      contrast_hat <- as.numeric(coef_bg[re_names[2]] - coef_bg[re_names[1]])
    } else {
      contrast_hat <- NA_real_
    }
    out$mark_contrast <- data.frame(
      contrast = "B_minus_A",
      true = TRUE_BETA_MARK_FULL,
      estimate = contrast_hat,
      error = contrast_hat - TRUE_BETA_MARK_FULL,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(grid_truth) && bg_case != "mark_cat") {
    pred_eta <- tryCatch(as.numeric(stats::predict(mod, newdata = grid_truth, type = "link")), error = function(e) NULL)
    if (!is.null(pred_eta) && length(pred_eta) == nrow(grid_truth) && all(is.finite(pred_eta))) {
      w <- grid_truth$w
      dens_hat <- exp(pred_eta)
      dens_hat <- dens_hat / sum(w * dens_hat, na.rm = TRUE)
      dens_true <- grid_truth$dens_bg_true_spatial
      eta_hat_c <- pred_eta - stats::weighted.mean(pred_eta, w = w, na.rm = TRUE)
      eta_true_c <- grid_truth$eta_bg_true_spatial - stats::weighted.mean(grid_truth$eta_bg_true_spatial, w = w, na.rm = TRUE)
      out$grid_shape <- data.frame(
        rmse_dens = sqrt(stats::weighted.mean((dens_hat - dens_true)^2, w = w, na.rm = TRUE)),
        mae_dens = stats::weighted.mean(abs(dens_hat - dens_true), w = w, na.rm = TRUE),
        cor_dens = suppressWarnings(stats::cor(dens_hat, dens_true, use = "complete.obs")),
        rmse_log_shape = sqrt(stats::weighted.mean((eta_hat_c - eta_true_c)^2, w = w, na.rm = TRUE)),
        cor_log_shape = suppressWarnings(stats::cor(eta_hat_c, eta_true_c, use = "complete.obs")),
        stringsAsFactors = FALSE
      )
    }
  }

  out
}

make_metric_payload <- function(sim_obj, fit_result_list = list(), bg_case = NULL, calibration = NULL) {
  s <- get_sim_core(sim_obj)
  catsim <- s$cat.sim[order(s$cat.sim$time), , drop = FALSE]
  rownames(catsim) <- NULL
  event_truth <- data.frame(
    event_id = catsim$event_id %||% seq_len(nrow(catsim)),
    time = catsim$time,
    long = catsim$long,
    lat = catsim$lat,
    x_km = catsim$x_km %||% NA_real_,
    y_km = catsim$y_km %||% NA_real_,
    magn1 = catsim$magn1,
    m_rel = catsim$m_rel %||% (catsim$magn1 - M0_FULL),
    father_id = catsim$father_id %||% NA_integer_,
    is_background_true = if ("father_id" %in% names(catsim)) catsim$father_id == 0 else NA,
    eta_bg_true = catsim$eta_bg_true %||% NA_real_,
    eta_trig_true = catsim$eta_trig_true %||% NA_real_,
    mark_cat = if ("mark_cat" %in% names(catsim)) as.character(catsim$mark_cat) else NA_character_,
    Z1 = if ("Z1" %in% names(catsim)) catsim$Z1 else NA_real_,
    stringsAsFactors = FALSE
  )

  true_params <- s$truth$params
  grid_truth <- if (!is.null(bg_case)) make_grid_eval(bg_case, grid_n = GRID_N) else NULL

  fit_summaries <- lapply(fit_result_list, function(x) {
    if (is.null(x) || !isTRUE(x$ok)) return(NULL)
    extract_basic_fit_summary(x$value)
  })

  fit_metrics <- lapply(fit_result_list, function(x) {
    if (is.null(x) || !isTRUE(x$ok)) return(NULL)
    fit <- x$value
    list(
      classification = fit_classification_metrics(fit, event_truth = event_truth, threshold = 0.5),
      parameter_comparison = fit_parameter_comparison(fit, true_params_vec = true_params, cat_sim_df = event_truth),
      trigger_components = fit_trigger_components(fit, true_params_vec = true_params, event_truth = event_truth),
      background = fit_background_metrics(fit, bg_case = bg_case, grid_truth = grid_truth)
    )
  })

  list(
    event_truth = event_truth,
    fit_summaries = fit_summaries,
    fit_metrics = fit_metrics,
    true_params = true_params,
    calibration = calibration,
    truth = s$truth,
    counts = c(n0 = s$n0, nson = s$nson, nobs = nrow(catsim), exploded = isTRUE(s$exploded))
  )
}

## ================================================================
## 8. CHECKPOINT LOGIC
## ================================================================

calibration_file <- function(block_dir, bg_case, count_case) {
  file.path(block_dir, "calibration", paste0(make_id("calib", bg_case, count_case), ".rds"))
}

result_file <- function(block_dir, bg_case, count_case, rep, suffix = NULL) {
  id <- make_id(bg_case, count_case, sprintf("rep%03d", as.integer(rep)), suffix %||% "")
  file.path(block_dir, "fits", paste0(id, ".rds"))
}

get_or_run_calibration <- function(block_dir, bg_case, count_case, seed_base,
                                   rerun = FALSE, verbose = TRUE) {
  bg_scenario <- BACKGROUND_SCENARIOS_FULL[[bg_case]]
  count_regime <- count_lookup(count_case)
  f <- calibration_file(block_dir, bg_case, count_case)

  if (file.exists(f) && !isTRUE(rerun)) {
    if (isTRUE(verbose)) log_msg("Calibration exists: %s | %s", bg_case, count_case)
    return(readRDS(f))
  }

  if (isTRUE(verbose)) log_msg("Calibrating: %s | %s", bg_case, count_case)
  cal <- calibrate_cell_full(
    bg_scenario = bg_scenario,
    count_regime = count_regime,
    seed_base = seed_base,
    nrep_pilot = NREP_CALIB,
    n_iters = N_CALIB_ITERS,
    n_bracket = N_CALIB_BRACKET,
    verbose = verbose
  )
  saveRDS_atomic(cal, f)
  cal
}

run_one_catalog_and_fits <- function(block_dir, bg_case, count_case, rep,
                                     fit_classic = TRUE,
                                     fit_param_correct = TRUE,
                                     fit_param_misspec = FALSE,
                                     force = FALSE,
                                     verbose = TRUE) {
  f <- result_file(block_dir, bg_case, count_case, rep)
  if (file.exists(f) && !isTRUE(force)) {
    obj <- readRDS(f)
    missing_classic <- fit_classic && is.null(obj$fits$classic)
    missing_param   <- fit_param_correct && is.null(obj$fits$param_correct)
    missing_mis     <- fit_param_misspec && is.null(obj$fits$param_misspec)
    if (!missing_classic && !missing_param && !missing_mis) {
      if (isTRUE(verbose)) log_msg("Skipping completed: %s | %s | rep %03d", bg_case, count_case, rep)
      return(obj)
    }
  } else {
    obj <- NULL
  }

  bg_scenario <- BACKGROUND_SCENARIOS_FULL[[bg_case]]
  count_regime <- count_lookup(count_case)

  seed_base <- SEED_MASTER_FULL + 1000000L * match(bg_case, names(BACKGROUND_SCENARIOS_FULL)) +
    10000L * match(count_case, unique(c(COUNT_MAIN_1000$count_case, COUNT_EXTREME_1000$count_case, SAMPLE_COUNTS$count_case)))

  cal <- get_or_run_calibration(block_dir, bg_case, count_case, seed_base = seed_base, verbose = verbose)

  if (is.null(obj)) {
    seed_rep <- seed_base + as.integer(rep)
    if (isTRUE(verbose)) log_msg("Simulating: %s | %s | rep %03d", bg_case, count_case, rep)
    sim_res <- safe_run(function() {
      simulate_one_full(
        bg_scenario = bg_scenario,
        count_regime = count_regime,
        seed = seed_rep,
        mu_value = cal$mu_cal,
        k0_value = cal$k0_cal,
        return_bg_info = TRUE,
        return_support = TRUE
      )
    })

    obj <- list(
      metadata = list(
        block_dir = block_dir,
        bg_case = bg_case,
        bg_label = bg_scenario$label,
        count_case = count_case,
        rep = as.integer(rep),
        seed_rep = seed_rep,
        created_at = Sys.time(),
        grid_n = GRID_N
      ),
      calibration = cal,
      simulation = sim_res,
      fits = list(),
      payload = NULL
    )
    saveRDS_atomic(obj, f)
  }

  if (!isTRUE(obj$simulation$ok)) {
    log_msg("Simulation failed: %s | %s | rep %03d | %s", bg_case, count_case, rep, obj$simulation$error_message, .level = "ERROR")
    return(obj)
  }

  cat_sim_df <- make_eqcat_from_sim(obj$simulation$value)
  if ("mark_cat" %in% names(get_cat_sim(obj$simulation$value))) {
    cat_sim_df$mark_cat <- factor(get_cat_sim(obj$simulation$value)$mark_cat, levels = c("A", "B"))
  }

  starts <- list(
    mu = obj$calibration$mu_cal,
    k0 = obj$calibration$k0_cal,
    c = TRUE_BASE_FULL["c"],
    p = TRUE_BASE_FULL["p"],
    gamma = 0,
    d = TRUE_BASE_FULL["d"],
    q = TRUE_BASE_FULL["q"],
    betacov = TRUE_BETACOV_FULL
  )

  if (isTRUE(fit_classic) && is.null(obj$fits$classic)) {
    if (isTRUE(verbose)) log_msg("Fitting CLASSIC: %s | %s | rep %03d", bg_case, count_case, rep)
    obj$fits$classic <- safe_run(function() fit_etas_classic_full(cat_sim_df, starts = starts))
    saveRDS_atomic(obj, f)
    log_msg("CLASSIC done: ok=%s | elapsed=%.2f sec", obj$fits$classic$ok, obj$fits$classic$elapsed_sec)
  }

  if (isTRUE(fit_param_correct) && is.null(obj$fits$param_correct)) {
    if (isTRUE(verbose)) log_msg("Fitting PARAM correct: %s | %s | rep %03d", bg_case, count_case, rep)
    obj$fits$param_correct <- safe_run(function() {
      fit_etas_parametric_full(cat_sim_df, bg_case = bg_case, formula_key = "correct", starts = starts, verbose.bg = FALSE, seed.bg = obj$metadata$seed_rep + 500000L)
    })
    saveRDS_atomic(obj, f)
    log_msg("PARAM correct done: ok=%s | elapsed=%.2f sec", obj$fits$param_correct$ok, obj$fits$param_correct$elapsed_sec)
  }

  if (isTRUE(fit_param_misspec) && is.null(obj$fits$param_misspec)) {
    if (isTRUE(verbose)) log_msg("Fitting PARAM misspecified: %s | %s | rep %03d", bg_case, count_case, rep)
    obj$fits$param_misspec <- safe_run(function() {
      fit_etas_parametric_full(cat_sim_df, bg_case = bg_case, formula_key = "misspecified", starts = starts, verbose.bg = FALSE, seed.bg = obj$metadata$seed_rep + 700000L)
    })
    saveRDS_atomic(obj, f)
    log_msg("PARAM misspecified done: ok=%s | elapsed=%.2f sec", obj$fits$param_misspec$ok, obj$fits$param_misspec$elapsed_sec)
  }

  ## Store metric-ready payload after the requested fits.
  obj$payload <- make_metric_payload(obj$simulation$value, obj$fits, bg_case = bg_case, calibration = obj$calibration)
  obj$grid_truth <- make_grid_eval(bg_case, grid_n = GRID_N)
  saveRDS_atomic(obj, f)

  log_file <- file.path(block_dir, "progress_log.csv")
  sim_counts <- obj$payload$counts
  append_csv_row(data.frame(
    time = as.character(Sys.time()),
    bg_case = bg_case,
    count_case = count_case,
    rep = as.integer(rep),
    n0 = as.numeric(sim_counts["n0"]),
    nson = as.numeric(sim_counts["nson"]),
    nobs = as.numeric(sim_counts["nobs"]),
    classic_ok = if (!is.null(obj$fits$classic)) obj$fits$classic$ok else NA,
    classic_elapsed = if (!is.null(obj$fits$classic)) obj$fits$classic$elapsed_sec else NA,
    param_correct_ok = if (!is.null(obj$fits$param_correct)) obj$fits$param_correct$ok else NA,
    param_correct_elapsed = if (!is.null(obj$fits$param_correct)) obj$fits$param_correct$elapsed_sec else NA,
    param_misspec_ok = if (!is.null(obj$fits$param_misspec)) obj$fits$param_misspec$ok else NA,
    param_misspec_elapsed = if (!is.null(obj$fits$param_misspec)) obj$fits$param_misspec$elapsed_sec else NA,
    stringsAsFactors = FALSE
  ), log_file)

  obj
}

run_plan_dataframe <- function(plan_df, block_name,
                               fit_classic = TRUE,
                               fit_param_correct = TRUE,
                               fit_param_misspec = FALSE,
                               force = FALSE,
                               verbose = TRUE) {
  block_dir <- file.path(OUT_ROOT, block_name)
  safe_dir_create(file.path(block_dir, "fits"))
  safe_dir_create(file.path(block_dir, "calibration"))

  log_msg("Starting block '%s' with %d rows", block_name, nrow(plan_df))

  for (ii in seq_len(nrow(plan_df))) {
    row <- plan_df[ii, ]
    log_msg("Block %s | row %d/%d | %s | %s | rep %03d",
            block_name, ii, nrow(plan_df), row$bg_case, row$count_case, as.integer(row$rep))
    run_one_catalog_and_fits(
      block_dir = block_dir,
      bg_case = row$bg_case,
      count_case = row$count_case,
      rep = row$rep,
      fit_classic = fit_classic,
      fit_param_correct = fit_param_correct,
      fit_param_misspec = fit_param_misspec,
      force = force,
      verbose = verbose
    )
    gc(verbose = FALSE)
  }

  log_msg("Completed block '%s'", block_name)
  invisible(block_dir)
}

## ================================================================
## 9. PUBLIC RUNNERS
## ================================================================

run_pilot_checks_full <- function(force = FALSE, verbose = TRUE) {
  stop("Pilot checks were removed from the final no-smooth simulation script. Use RUN_MODE='main', 'sample', 'misspec', 'extreme', or 'all'.")
}

run_main_plan_full <- function(force = FALSE, verbose = TRUE) {
  run_plan_dataframe(
    plan_df = MAIN_PLAN,
    block_name = "01_main_plan",
    fit_classic = TRUE,
    fit_param_correct = TRUE,
    fit_param_misspec = FALSE,
    force = force,
    verbose = verbose
  )
}

run_secondary_sample_size_full <- function(force = FALSE, verbose = TRUE) {
  run_plan_dataframe(
    plan_df = SAMPLE_PLAN,
    block_name = "02_secondary_sample_size",
    fit_classic = TRUE,
    fit_param_correct = TRUE,
    fit_param_misspec = FALSE,
    force = force,
    verbose = verbose
  )
}

run_one_misspec_from_main <- function(block_dir, bg_case, count_case, rep,
                                      main_block_dir = file.path(OUT_ROOT, "01_main_plan"),
                                      force = FALSE,
                                      verbose = TRUE) {
  out_file <- result_file(block_dir, bg_case, count_case, rep, suffix = "misspec")
  if (file.exists(out_file) && !isTRUE(force)) {
    obj <- readRDS(out_file)
    if (!is.null(obj$fits$param_misspec)) {
      if (isTRUE(verbose)) log_msg("Skipping completed MISSPEC: %s | %s | rep %03d", bg_case, count_case, rep)
      return(obj)
    }
  }

  main_file <- result_file(main_block_dir, bg_case, count_case, rep)
  if (!file.exists(main_file)) {
    stop(
      "Main-plan result not found for misspecification block: ", main_file,
      "\nRun run_main_plan_full() first, because misspecification is fitted on the same catalogues as the main plan."
    )
  }

  main_obj <- readRDS(main_file)
  if (is.null(main_obj$simulation) || !isTRUE(main_obj$simulation$ok)) {
    stop("Main-plan simulation is missing or failed: ", main_file)
  }

  obj <- list(
    metadata = c(main_obj$metadata, list(
      block_dir = block_dir,
      source_main_file = main_file,
      misspecification = TRUE,
      created_at_misspec = Sys.time()
    )),
    calibration = main_obj$calibration,
    simulation = main_obj$simulation,
    fits = list(),
    payload = NULL,
    grid_truth = main_obj$grid_truth %||% make_grid_eval(bg_case, grid_n = GRID_N)
  )
  saveRDS_atomic(obj, out_file)

  cat_sim_df <- make_eqcat_from_sim(obj$simulation$value)
  if ("mark_cat" %in% names(get_cat_sim(obj$simulation$value))) {
    cat_sim_df$mark_cat <- factor(get_cat_sim(obj$simulation$value)$mark_cat, levels = c("A", "B"))
  }

  starts <- list(
    mu = obj$calibration$mu_cal,
    k0 = obj$calibration$k0_cal,
    c = TRUE_BASE_FULL["c"],
    p = TRUE_BASE_FULL["p"],
    gamma = 0,
    d = TRUE_BASE_FULL["d"],
    q = TRUE_BASE_FULL["q"],
    betacov = TRUE_BETACOV_FULL
  )

  if (isTRUE(verbose)) log_msg("Fitting PARAM misspecified on MAIN catalogue: %s | %s | rep %03d", bg_case, count_case, rep)
  obj$fits$param_misspec <- safe_run(function() {
    fit_etas_parametric_full(cat_sim_df, bg_case = bg_case, formula_key = "misspecified", starts = starts, verbose.bg = FALSE, seed.bg = obj$metadata$seed_rep + 700000L)
  })
  obj$payload <- make_metric_payload(obj$simulation$value, obj$fits, bg_case = bg_case, calibration = obj$calibration)
  saveRDS_atomic(obj, out_file)

  log_file <- file.path(block_dir, "progress_log.csv")
  sim_counts <- obj$payload$counts
  append_csv_row(data.frame(
    time = as.character(Sys.time()),
    bg_case = bg_case,
    count_case = count_case,
    rep = as.integer(rep),
    n0 = as.numeric(sim_counts["n0"]),
    nson = as.numeric(sim_counts["nson"]),
    nobs = as.numeric(sim_counts["nobs"]),
    classic_ok = NA,
    classic_elapsed = NA,
    param_correct_ok = NA,
    param_correct_elapsed = NA,
    param_misspec_ok = obj$fits$param_misspec$ok,
    param_misspec_elapsed = obj$fits$param_misspec$elapsed_sec,
    source_main_file = main_file,
    stringsAsFactors = FALSE
  ), log_file)

  log_msg("PARAM misspecified done: ok=%s | elapsed=%.2f sec", obj$fits$param_misspec$ok, obj$fits$param_misspec$elapsed_sec)
  obj
}

run_misspec_plan_dataframe <- function(plan_df, block_name = "03_secondary_misspecification",
                                       force = FALSE, verbose = TRUE) {
  block_dir <- file.path(OUT_ROOT, block_name)
  safe_dir_create(file.path(block_dir, "fits"))
  log_msg("Starting block '%s' with %d rows", block_name, nrow(plan_df))

  for (ii in seq_len(nrow(plan_df))) {
    row <- plan_df[ii, ]
    log_msg("Block %s | row %d/%d | %s | %s | rep %03d",
            block_name, ii, nrow(plan_df), row$bg_case, row$count_case, as.integer(row$rep))
    run_one_misspec_from_main(
      block_dir = block_dir,
      bg_case = row$bg_case,
      count_case = row$count_case,
      rep = row$rep,
      force = force,
      verbose = verbose
    )
    gc(verbose = FALSE)
  }

  log_msg("Completed block '%s'", block_name)
  invisible(block_dir)
}

run_secondary_misspecification_full <- function(force = FALSE, verbose = TRUE) {
  ## This block intentionally fits only the misspecified parametric model
  ## on the SAME catalogues generated in 01_main_plan. The classic and
  ## correctly specified parametric benchmarks are therefore reused from
  ## the main-plan RDS files and are not refitted here.
  run_misspec_plan_dataframe(
    plan_df = MISSPEC_PLAN,
    block_name = "03_secondary_misspecification",
    force = force,
    verbose = verbose
  )
}

run_secondary_extreme_triggering_full <- function(force = FALSE, verbose = TRUE) {
  run_plan_dataframe(
    plan_df = EXTREME_PLAN,
    block_name = "04_secondary_extreme_triggering",
    fit_classic = TRUE,
    fit_param_correct = TRUE,
    fit_param_misspec = FALSE,
    force = force,
    verbose = verbose
  )
}

## ================================================================
## 10. QUICK SUMMARY HELPERS
## ================================================================

collect_progress_logs <- function(out_root = OUT_ROOT) {
  files <- list.files(out_root, pattern = "progress_log\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(data.frame())
  out <- lapply(files, function(f) {
    x <- utils::read.csv(f, stringsAsFactors = FALSE)
    x$block <- basename(dirname(f))
    x
  })
  do.call(rbind_fill_base, out)
}

quick_counts_summary <- function(out_root = OUT_ROOT) {
  x <- collect_progress_logs(out_root)
  if (!nrow(x)) return(x)
  stats::aggregate(
    cbind(n0, nson, nobs, classic_elapsed, param_correct_elapsed, param_misspec_elapsed) ~ block + bg_case + count_case,
    data = x,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
}

quick_success_summary <- function(out_root = OUT_ROOT) {
  x <- collect_progress_logs(out_root)
  if (!nrow(x)) return(x)
  stats::aggregate(
    cbind(classic_ok, param_correct_ok, param_misspec_ok) ~ block + bg_case + count_case,
    data = x,
    FUN = function(z) mean(as.numeric(z), na.rm = TRUE)
  )
}

print_plan_sizes <- function() {
  cat("\n==================== PLAN SIZES ====================\n")
  cat("Scenarios:       ", paste(FINAL_BG_CASES, collapse = ", "), "\n")
  cat("Pilot rows:      removed from final plan\n")
  cat("Main rows:       ", nrow(MAIN_PLAN), "\n")
  cat("Sample rows:     ", nrow(SAMPLE_PLAN), "\n")
  cat("Misspec rows:    ", nrow(MISSPEC_PLAN), "\n")
  cat("Extreme rows:    ", nrow(EXTREME_PLAN), "\n")
  cat("Output root:     ", OUT_ROOT, "\n")
  cat("====================================================\n\n")
}

## ================================================================
## 11. AUTORUN SWITCH
## ================================================================

print_plan_sizes()

RUN_MODE = "main" # FATTO
RUN_MODE = "sample" # DA ESPLORARE


OUT_ROOT <- "/home/nicolettadangelo/sim_paper1_marco/etas parametrico/etas_parametric_final_no_smooth_outputs_v1"

RUN_MODE = "misspec"

if (identical(RUN_MODE, "none")) {
  log_msg("RUN_MODE='none': functions loaded, no simulations launched.")
} else if (identical(RUN_MODE, "main")) {
  log_msg("RUN_MODE='main': running main plan only.")
  run_main_plan_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "sample")) {
  log_msg("RUN_MODE='sample': running secondary sample-size plan only.")
  run_secondary_sample_size_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "misspec")) {
  log_msg("RUN_MODE='misspec': running secondary misspecification plan only.")
  run_secondary_misspecification_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "extreme")) {
  log_msg("RUN_MODE='extreme': running secondary extreme-triggering plan only.")
  run_secondary_extreme_triggering_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "all")) {
  log_msg("RUN_MODE='all': running all final blocks sequentially: main, sample-size, misspecification, extreme-triggering.")
  run_main_plan_full(force = FALSE, verbose = TRUE)
  run_secondary_sample_size_full(force = FALSE, verbose = TRUE)
  run_secondary_misspecification_full(force = FALSE, verbose = TRUE)
  run_secondary_extreme_triggering_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "pilot")) {
  stop("RUN_MODE='pilot' is not available in the final no-smooth script: the pilot block was removed.")
} else {
  stop("Unknown ETAS_RUN_MODE: ", RUN_MODE)
}


RUN_MODE = "extreme"

if (identical(RUN_MODE, "none")) {
  log_msg("RUN_MODE='none': functions loaded, no simulations launched.")
} else if (identical(RUN_MODE, "main")) {
  log_msg("RUN_MODE='main': running main plan only.")
  run_main_plan_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "sample")) {
  log_msg("RUN_MODE='sample': running secondary sample-size plan only.")
  run_secondary_sample_size_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "misspec")) {
  log_msg("RUN_MODE='misspec': running secondary misspecification plan only.")
  run_secondary_misspecification_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "extreme")) {
  log_msg("RUN_MODE='extreme': running secondary extreme-triggering plan only.")
  run_secondary_extreme_triggering_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "all")) {
  log_msg("RUN_MODE='all': running all final blocks sequentially: main, sample-size, misspecification, extreme-triggering.")
  run_main_plan_full(force = FALSE, verbose = TRUE)
  run_secondary_sample_size_full(force = FALSE, verbose = TRUE)
  run_secondary_misspecification_full(force = FALSE, verbose = TRUE)
  run_secondary_extreme_triggering_full(force = FALSE, verbose = TRUE)
} else if (identical(RUN_MODE, "pilot")) {
  stop("RUN_MODE='pilot' is not available in the final no-smooth script: the pilot block was removed.")
} else {
  stop("Unknown ETAS_RUN_MODE: ", RUN_MODE)
}



############################################################
## MAIN PLAN RESULTS EXPLORATION
############################################################

rbind_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && is.data.frame(z) && nrow(z) > 0, x)
  if (!length(x)) return(data.frame())
  
  all_names <- unique(unlist(lapply(x, names)))
  
  x <- lapply(x, function(df) {
    miss <- setdiff(all_names, names(df))
    for (m in miss) df[[m]] <- NA
    df <- df[, all_names, drop = FALSE]
    rownames(df) <- NULL
    df
  })
  
  do.call(rbind, x)
}

safe_readRDS <- function(file) {
  tryCatch(
    readRDS(file),
    error = function(e) {
      warning("Could not read: ", file, "\n", conditionMessage(e))
      NULL
    }
  )
}

extract_main_one <- function(obj, file) {
  
  meta <- obj$metadata %||% list()
  payload <- obj$payload %||% list()
  fits <- obj$fits %||% list()
  
  bg_case <- as.character(meta$bg_case %||% NA)
  count_case <- as.character(meta$count_case %||% NA)
  rep_id <- as.integer(meta$rep %||% NA)
  
  counts <- payload$counts %||% c(n0 = NA, nson = NA, nobs = NA, exploded = NA)
  
  index <- data.frame(
    file = file,
    bg_case = bg_case,
    count_case = count_case,
    rep = rep_id,
    seed_rep = as.integer(meta$seed_rep %||% NA),
    n0 = as.numeric(counts["n0"] %||% NA),
    nson = as.numeric(counts["nson"] %||% NA),
    nobs = as.numeric(counts["nobs"] %||% NA),
    frac_bg = as.numeric(counts["n0"] %||% NA) / as.numeric(counts["nobs"] %||% NA),
    frac_trig = as.numeric(counts["nson"] %||% NA) / as.numeric(counts["nobs"] %||% NA),
    exploded = as.logical(counts["exploded"] %||% NA),
    stringsAsFactors = FALSE
  )
  
  fit_status <- rbind_fill(lapply(names(fits), function(model_name) {
    fit_i <- fits[[model_name]]
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      ok = isTRUE(fit_i$ok),
      elapsed_sec = as.numeric(fit_i$elapsed_sec %||% NA),
      error = as.character(fit_i$error %||% NA),
      stringsAsFactors = FALSE
    )
  }))
  
  fit_metrics <- payload$fit_metrics %||% list()
  
  parameter_comparison <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    pc <- fit_metrics[[model_name]]$parameter_comparison
    if (is.null(pc) || !is.data.frame(pc)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      pc,
      stringsAsFactors = FALSE
    )
  }))
  
  classification_counts <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    cs <- fit_metrics[[model_name]]$classification$count_summary
    if (is.null(cs) || !is.data.frame(cs)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      cs,
      stringsAsFactors = FALSE
    )
  }))
  
  classification_metrics <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    cm <- fit_metrics[[model_name]]$classification$metrics %||%
      fit_metrics[[model_name]]$classification$summary_metrics
    if (is.null(cm) || !is.data.frame(cm)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      cm,
      stringsAsFactors = FALSE
    )
  }))
  
  trigger_components <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    tc <- fit_metrics[[model_name]]$trigger_components
    if (is.null(tc) || !is.data.frame(tc)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      tc,
      stringsAsFactors = FALSE
    )
  }))
  
  background_coef_comparison <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    bg <- fit_metrics[[model_name]]$background
    cc <- bg$coef_comparison
    if (is.null(cc) || !is.data.frame(cc)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      cc,
      stringsAsFactors = FALSE
    )
  }))
  
  mark_contrast <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    bg <- fit_metrics[[model_name]]$background
    mc <- bg$mark_contrast
    if (is.null(mc) || !is.data.frame(mc)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      mc,
      stringsAsFactors = FALSE
    )
  }))
  
  background_grid_shape <- rbind_fill(lapply(names(fit_metrics), function(model_name) {
    bg <- fit_metrics[[model_name]]$background
    gs <- bg$grid_shape
    if (is.null(gs) || !is.data.frame(gs)) return(NULL)
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      gs,
      stringsAsFactors = FALSE
    )
  }))
  
  fit_summaries <- payload$fit_summaries %||% list()
  
  iterations <- rbind_fill(lapply(names(fit_summaries), function(model_name) {
    fs <- fit_summaries[[model_name]]
    pi <- fs$params_iter
    if (is.null(pi)) return(NULL)
    
    pi <- as.data.frame(pi, check.names = FALSE)
    pi$iteration <- seq_len(nrow(pi))
    
    aic <- fs$AIC_iter %||% rep(NA_real_, nrow(pi))
    if (length(aic) < nrow(pi)) {
      aic <- c(aic, rep(NA_real_, nrow(pi) - length(aic)))
    }
    
    pi$AIC <- as.numeric(aic[seq_len(nrow(pi))])
    
    data.frame(
      file = file,
      bg_case = bg_case,
      count_case = count_case,
      rep = rep_id,
      model = model_name,
      pi,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  
  list(
    index = index,
    fit_status = fit_status,
    parameter_comparison = parameter_comparison,
    classification_counts = classification_counts,
    classification_metrics = classification_metrics,
    trigger_components = trigger_components,
    background_coef_comparison = background_coef_comparison,
    mark_contrast = mark_contrast,
    background_grid_shape = background_grid_shape,
    iterations = iterations
  )
}

collect_main_results <- function(out_root = OUT_ROOT) {
  
  main_dir <- file.path(out_root, "01_main_plan", "fits")
  
  files <- list.files(
    main_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  
  if (!length(files)) {
    stop("No .rds files found in: ", main_dir)
  }
  
  message("Reading ", length(files), " main-plan files...")
  
  objs <- lapply(files, safe_readRDS)
  
  extracted <- Map(
    function(obj, file) {
      if (is.null(obj)) return(NULL)
      extract_main_one(obj, file)
    },
    objs,
    files
  )
  
  extracted <- Filter(Negate(is.null), extracted)
  
  tabs <- list(
    index = rbind_fill(lapply(extracted, `[[`, "index")),
    fit_status = rbind_fill(lapply(extracted, `[[`, "fit_status")),
    parameter_comparison = rbind_fill(lapply(extracted, `[[`, "parameter_comparison")),
    classification_counts = rbind_fill(lapply(extracted, `[[`, "classification_counts")),
    classification_metrics = rbind_fill(lapply(extracted, `[[`, "classification_metrics")),
    trigger_components = rbind_fill(lapply(extracted, `[[`, "trigger_components")),
    background_coef_comparison = rbind_fill(lapply(extracted, `[[`, "background_coef_comparison")),
    mark_contrast = rbind_fill(lapply(extracted, `[[`, "mark_contrast")),
    background_grid_shape = rbind_fill(lapply(extracted, `[[`, "background_grid_shape")),
    iterations = rbind_fill(lapply(extracted, `[[`, "iterations"))
  )
  
  tabs
}

main_tabs <- collect_main_results(OUT_ROOT)

index_tbl <- main_tabs$index
fit_status_tbl <- main_tabs$fit_status
param_tbl <- main_tabs$parameter_comparison
class_counts_tbl <- main_tabs$classification_counts
class_metrics_tbl <- main_tabs$classification_metrics
trigger_tbl <- main_tabs$trigger_components
bg_coef_tbl <- main_tabs$background_coef_comparison
mark_contrast_tbl <- main_tabs$mark_contrast
bg_shape_tbl <- main_tabs$background_grid_shape
iter_tbl <- main_tabs$iterations


############################################################
## 1. COMPLETION AND SUCCESS TABLES
############################################################

expected_main <- expand.grid(
  bg_case = c("constant", "linear_xy", "cov1", "mark_cat", "hotspot"),
  count_case = c(
    "N1000_bgdom_75_25",
    "N1000_balanced_50_50",
    "N1000_trigdom_25_75"
  ),
  rep = 1:50,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

completion_tbl <- expected_main %>%
  left_join(
    index_tbl %>% mutate(done = TRUE) %>% select(bg_case, count_case, rep, done),
    by = c("bg_case", "count_case", "rep")
  ) %>%
  mutate(done = ifelse(is.na(done), FALSE, done))

completion_summary <- completion_tbl %>%
  group_by(bg_case, count_case) %>%
  summarise(
    expected_reps = n(),
    completed_reps = sum(done),
    missing_reps = sum(!done),
    completion_rate = mean(done),
    .groups = "drop"
  )

completion_summary

fit_success_summary <- fit_status_tbl %>%
  group_by(bg_case, count_case, model) %>%
  summarise(
    n_fit = n(),
    n_ok = sum(ok, na.rm = TRUE),
    success_rate = mean(ok, na.rm = TRUE),
    mean_elapsed_sec = mean(elapsed_sec, na.rm = TRUE),
    median_elapsed_sec = median(elapsed_sec, na.rm = TRUE),
    .groups = "drop"
  )

fit_success_summary

failed_fits <- fit_status_tbl %>%
  filter(!ok | is.na(ok)) %>%
  arrange(bg_case, count_case, rep, model)

failed_fits



############################################################
## 2. SIMULATED COUNTS
############################################################

sim_counts_summary <- index_tbl %>%
  group_by(bg_case, count_case) %>%
  summarise(
    n_rep = n(),
    mean_n0 = mean(n0, na.rm = TRUE),
    sd_n0 = sd(n0, na.rm = TRUE),
    mean_nson = mean(nson, na.rm = TRUE),
    sd_nson = sd(nson, na.rm = TRUE),
    mean_nobs = mean(nobs, na.rm = TRUE),
    sd_nobs = sd(nobs, na.rm = TRUE),
    mean_frac_bg = mean(frac_bg, na.rm = TRUE),
    sd_frac_bg = sd(frac_bg, na.rm = TRUE),
    mean_frac_trig = mean(frac_trig, na.rm = TRUE),
    sd_frac_trig = sd(frac_trig, na.rm = TRUE),
    n_exploded = sum(exploded, na.rm = TRUE),
    .groups = "drop"
  )

sim_counts_summary


############################################################
## 3. ETAS PARAMETER PERFORMANCE
############################################################

param_summary <- param_tbl %>%
  group_by(bg_case, count_case, model, parameter) %>%
  summarise(
    n = sum(!is.na(estimate)),
    true = mean(true, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    abs_bias = mean(abs(error), na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    median_rel_error = median(rel_error, na.rm = TRUE),
    .groups = "drop"
  )

param_summary

param_summary %>%
  filter(model == "param_correct") %>%
  arrange(bg_case, count_case, parameter)


############################################################
## 4. CLASSIFICATION PERFORMANCE
############################################################

class_metric_summary <- class_metrics_tbl %>%
  group_by(bg_case, count_case, model, metric) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

class_metric_summary

class_metric_summary %>%
  filter(metric %in% c("accuracy_05", "auc_bg", "brier_bg", "logloss_bg")) %>%
  arrange(bg_case, count_case, model, metric)

count_error_summary <- class_counts_tbl %>%
  group_by(bg_case, count_case, model, class) %>%
  summarise(
    n = n(),
    mean_true_count = mean(true_count, na.rm = TRUE),
    mean_expected_count_hat = mean(expected_count_hat, na.rm = TRUE),
    mean_difference = mean(difference, na.rm = TRUE),
    sd_difference = sd(difference, na.rm = TRUE),
    rmse_difference = sqrt(mean(difference^2, na.rm = TRUE)),
    .groups = "drop"
  )

count_error_summary


############################################################
## 5. TRIGGERING COMPONENTS
############################################################

trigger_summary <- trigger_tbl %>%
  group_by(bg_case, count_case, model, component) %>%
  summarise(
    n = n(),
    true = mean(true, na.rm = TRUE),
    mean_estimated = mean(estimated, na.rm = TRUE),
    sd_estimated = sd(estimated, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    median_rel_error = median(rel_error, na.rm = TRUE),
    .groups = "drop"
  )

trigger_summary

trigger_summary %>%
  filter(component == "approx_branching_component") %>%
  arrange(bg_case, count_case, model)


############################################################
## 6. BACKGROUND PERFORMANCE
############################################################

bg_coef_summary <- bg_coef_tbl %>%
  group_by(bg_case, count_case, model, coefficient) %>%
  summarise(
    n = n(),
    true = mean(true, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    .groups = "drop"
  )

bg_coef_summary

mark_contrast_summary <- mark_contrast_tbl %>%
  group_by(bg_case, count_case, model, contrast) %>%
  summarise(
    n = n(),
    true = mean(true, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    .groups = "drop"
  )

mark_contrast_summary

bg_shape_summary <- bg_shape_tbl %>%
  group_by(bg_case, count_case, model) %>%
  summarise(
    n = n(),
    mean_rmse_dens = mean(rmse_dens, na.rm = TRUE),
    sd_rmse_dens = sd(rmse_dens, na.rm = TRUE),
    mean_mae_dens = mean(mae_dens, na.rm = TRUE),
    mean_cor_dens = mean(cor_dens, na.rm = TRUE),
    mean_rmse_log_shape = mean(rmse_log_shape, na.rm = TRUE),
    mean_cor_log_shape = mean(cor_log_shape, na.rm = TRUE),
    .groups = "drop"
  )

bg_shape_summary


############################################################
## PLOT SETTINGS
############################################################

count_case_labs <- c(
  N1000_bgdom_75_25 = "BG-dominant 75/25",
  N1000_balanced_50_50 = "Balanced 50/50",
  N1000_trigdom_25_75 = "TR-dominant 25/75"
)

bg_case_labs <- c(
  constant = "Constant",
  linear_xy = "Linear x/y",
  cov1 = "Covariate Z1",
  mark_cat = "Categorical mark",
  hotspot = "Hotspot"
)

plot_dir <- file.path(OUT_ROOT, "01_main_plan", "diagnostic_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts_long <- index_tbl %>%
  select(bg_case, count_case, rep, n0, nson, nobs) %>%
  pivot_longer(
    cols = c(n0, nson, nobs),
    names_to = "count_type",
    values_to = "count"
  )

p_counts <- ggplot(
  counts_long,
  aes(x = count_case, y = count)
) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_grid(count_type ~ bg_case, scales = "free_y", labeller = labeller(
    bg_case = bg_case_labs
  )) +
  scale_x_discrete(labels = count_case_labs) +
  labs(
    title = "Simulated event counts by scenario",
    x = NULL,
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_counts


class_plot_data <- class_metrics_tbl %>%
  filter(metric %in% c("accuracy_05", "auc_bg", "brier_bg", "logloss_bg")) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("accuracy_05", "auc_bg", "brier_bg", "logloss_bg"),
      labels = c("Accuracy", "AUC", "Brier score", "Log loss")
    )
  )

p_class <- ggplot(
  class_plot_data,
  aes(x = count_case, y = value, fill = model)
) +
  geom_boxplot(outlier.alpha = 0.25, position = position_dodge(width = 0.8)) +
  facet_grid(metric ~ bg_case, scales = "free_y", labeller = labeller(
    bg_case = bg_case_labs
  )) +
  scale_x_discrete(labels = count_case_labs) +
  labs(
    title = "Background/triggered classification performance",
    x = NULL,
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_class


param_plot_data <- param_tbl %>%
  filter(parameter %in% c("mu", "k0", "c", "p", "d", "q", "magnitude")) %>%
  mutate(
    rel_error_clip = pmax(pmin(rel_error, 2), -2)
  )

p_param_rel <- ggplot(
  param_plot_data,
  aes(x = parameter, y = rel_error_clip, fill = model)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.alpha = 0.25, position = position_dodge(width = 0.8)) +
  facet_grid(count_case ~ bg_case, labeller = labeller(
    bg_case = bg_case_labs,
    count_case = count_case_labs
  )) +
  labs(
    title = "Relative error of ETAS parameter estimates",
    subtitle = "Relative errors clipped to [-2, 2] for readability",
    x = NULL,
    y = "Relative error",
    fill = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_param_rel


p_trigger_branch <- trigger_tbl %>%
  filter(component == "approx_branching_component") %>%
  ggplot(aes(x = count_case, y = rel_error, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.alpha = 0.25, position = position_dodge(width = 0.8)) +
  facet_wrap(~ bg_case, labeller = labeller(bg_case = bg_case_labs)) +
  scale_x_discrete(labels = count_case_labs) +
  scale_y_continuous(limits = c(-1,1)) +
  labs(
    title = "Relative error in integrated triggering component",
    x = NULL,
    y = "Relative error",
    fill = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_trigger_branch



############################################################
## BACKGROUND SHAPE RECOVERY: SEPARATE RMSE PLOTS
############################################################

if (nrow(bg_shape_tbl) > 0) {
  
  ############################################################
  ## 1. RMSE density
  ############################################################
  
  bg_rmse_dens <- bg_shape_tbl %>%
    select(bg_case, count_case, rep, model, rmse_dens) %>%
    filter(is.finite(rmse_dens))
  
  y_rmse_dens_max <- quantile(bg_rmse_dens$rmse_dens, 0.95, na.rm = TRUE) * 1.15
  
  if (!is.finite(y_rmse_dens_max) || y_rmse_dens_max <= 0) {
    y_rmse_dens_max <- max(bg_rmse_dens$rmse_dens, na.rm = TRUE)
  }
  
  p_rmse_dens <- ggplot(
    bg_rmse_dens,
    aes(x = count_case, y = rmse_dens, fill = model)
  ) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_boxplot(
      outlier.alpha = 0.25,
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap(
      ~ bg_case,
      scales = "free_y",
      labeller = labeller(bg_case = bg_case_labs)
    ) +
    scale_x_discrete(labels = count_case_labs) +
    coord_cartesian(ylim = c(0, y_rmse_dens_max)) +
    labs(
      title = "Background density recovery: RMSE",
      subtitle = "Values closer to 0 indicate better recovery",
      x = NULL,
      y = "RMSE of density",
      fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    )
  
  p_rmse_dens
  
  
  ############################################################
  ## 2. RMSE log-shape
  ############################################################
  
  bg_rmse_log <- bg_shape_tbl %>%
    select(bg_case, count_case, rep, model, rmse_log_shape) %>%
    filter(is.finite(rmse_log_shape))
  
  y_rmse_log_max <- quantile(bg_rmse_log$rmse_log_shape, 0.95, na.rm = TRUE) * 1.15
  
  if (!is.finite(y_rmse_log_max) || y_rmse_log_max <= 0) {
    y_rmse_log_max <- max(bg_rmse_log$rmse_log_shape, na.rm = TRUE)
  }
  
  p_rmse_log <- ggplot(
    bg_rmse_log,
    aes(x = count_case, y = rmse_log_shape, fill = model)
  ) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_boxplot(
      outlier.alpha = 0.25,
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap(
      ~ bg_case,
      scales = "free_y",
      labeller = labeller(bg_case = bg_case_labs)
    ) +
    scale_x_discrete(labels = count_case_labs) +
    coord_cartesian(ylim = c(0, y_rmse_log_max)) +
    labs(
      title = "Background log-shape recovery: RMSE",
      subtitle = "Values closer to 0 indicate better recovery",
      x = NULL,
      y = "RMSE of log-shape",
      fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    )
  
  p_rmse_log
}


if (nrow(bg_coef_tbl) > 0) {
  
  p_bg_coef <- bg_coef_tbl %>%
    ggplot(aes(x = coefficient, y = estimate, fill = model)) +
    geom_hline(aes(yintercept = true), linetype = 2) +
    geom_boxplot(outlier.alpha = 0.25, position = position_dodge(width = 0.8)) +
    scale_y_continuous(limits = c(-1,1)) +
    facet_grid(count_case ~ bg_case, scales = "free_x", labeller = labeller(
      bg_case = bg_case_labs,
      count_case = count_case_labs
    )) +
    labs(
      title = "Background coefficient recovery",
      subtitle = "Dashed lines are true values",
      x = NULL,
      y = "Estimate",
      fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  p_bg_coef
  
}


if (nrow(mark_contrast_tbl) > 0) {
  
  p_mark_contrast <- ggplot(
    mark_contrast_tbl,
    aes(x = count_case, y = estimate, fill = model)
  ) +
    geom_hline(aes(yintercept = true), linetype = 2) +
    geom_boxplot(outlier.alpha = 0.25, position = position_dodge(width = 0.8)) +
    scale_y_continuous(limits = c(0,1)) +
    facet_wrap(~ bg_case, labeller = labeller(bg_case = bg_case_labs)) +
    scale_x_discrete(labels = count_case_labs) +
    labs(
      title = "Categorical-mark background contrast recovery",
      subtitle = "Dashed line is the true contrast",
      x = NULL,
      y = "Estimated contrast",
      fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  p_mark_contrast

}


############################################################
## EXPECTED BACKGROUND/TRIGGERED COUNT ERROR
############################################################

count_error_plot_data <- class_counts_tbl %>%
  mutate(
    class = factor(
      class,
      levels = c("background", "triggered"),
      labels = c("Background", "Triggered")
    ),
    bg_case = factor(
      bg_case,
      levels = c("constant", "linear_xy", "cov1", "mark_cat", "hotspot")
    ),
    count_case = factor(
      count_case,
      levels = c(
        "N1000_balanced_50_50",
        "N1000_bgdom_75_25",
        "N1000_trigdom_25_75"
      )
    )
  )

count_error_summary <- count_error_plot_data %>%
  group_by(bg_case, count_case, model, class) %>%
  summarise(
    n = n(),
    mean_true_count = mean(true_count, na.rm = TRUE),
    mean_expected_count_hat = mean(expected_count_hat, na.rm = TRUE),
    mean_difference = mean(difference, na.rm = TRUE),
    sd_difference = sd(difference, na.rm = TRUE),
    rmse_difference = sqrt(mean(difference^2, na.rm = TRUE)),
    q10_difference = quantile(difference, 0.10, na.rm = TRUE),
    q90_difference = quantile(difference, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

count_error_summary


############################################################
## BACKGROUND SHAPE RECOVERY: TRUE / CLASSIC / PARAMETRIC
############################################################

get_fit_value <- function(obj, model_name) {
  if (is.null(obj$fits[[model_name]])) return(NULL)
  if (!isTRUE(obj$fits[[model_name]]$ok)) return(NULL)
  obj$fits[[model_name]]$value
}

weighted_kde_grid <- function(x, y, weights, grid, h = NULL, chunk_size = 250L) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  
  ok <- is.finite(x) & is.finite(y) & is.finite(weights) & weights >= 0
  x <- x[ok]
  y <- y[ok]
  weights <- weights[ok]
  
  if (!length(x)) {
    return(rep(NA_real_, nrow(grid)))
  }
  
  if (is.null(h) || length(h) < 2 || any(!is.finite(h)) || any(h <= 0)) {
    hx <- stats::bw.nrd0(x)
    hy <- stats::bw.nrd0(y)
    h <- c(hx, hy)
  }
  
  hx <- as.numeric(h[1])
  hy <- as.numeric(h[2])
  
  if (!is.finite(hx) || hx <= 0) hx <- diff(range(x, na.rm = TRUE)) / 20
  if (!is.finite(hy) || hy <= 0) hy <- diff(range(y, na.rm = TRUE)) / 20
  
  gx <- as.numeric(grid$x)
  gy <- as.numeric(grid$y)
  
  dens <- numeric(length(gx))
  idx <- split(seq_along(x), ceiling(seq_along(x) / chunk_size))
  
  for (ii in idx) {
    dx <- outer(gx, x[ii], function(a, b) stats::dnorm((a - b) / hx) / hx)
    dy <- outer(gy, y[ii], function(a, b) stats::dnorm((a - b) / hy) / hy)
    dens <- dens + rowSums(dx * dy * matrix(weights[ii], nrow = length(gx), ncol = length(ii), byrow = TRUE))
  }
  
  dens
}

normalize_density_on_grid <- function(dens, grid) {
  dens <- as.numeric(dens)
  w <- grid$w
  
  dens[!is.finite(dens)] <- NA_real_
  
  integ <- sum(w * dens, na.rm = TRUE)
  
  if (!is.finite(integ) || integ <= 0) {
    return(rep(NA_real_, length(dens)))
  }
  
  dens / integ
}

center_log_shape <- function(dens, grid, eps = 1e-300) {
  logd <- log(pmax(dens, eps))
  logd - weighted.mean(logd, w = grid$w, na.rm = TRUE)
}

estimate_classic_background_grid <- function(fit, grid) {
  
  if (is.null(fit) || is.null(fit$cat)) {
    return(NULL)
  }
  
  df <- as.data.frame(fit$cat)
  
  x <- if ("xcat.work" %in% names(df)) df$xcat.work else if ("x_km" %in% names(df)) df$x_km else df$long
  y <- if ("ycat.work" %in% names(df)) df$ycat.work else if ("y_km" %in% names(df)) df$y_km else df$lat
  
  rho <- fit$rho.weights
  
  if (is.null(rho) || length(rho) != nrow(df)) {
    rho <- rep(1, nrow(df))
  }
  
  h <- fit$hdef
  
  dens_raw <- weighted_kde_grid(
    x = x,
    y = y,
    weights = rho,
    grid = grid,
    h = h
  )
  
  dens <- normalize_density_on_grid(dens_raw, grid)
  
  data.frame(
    x = grid$x,
    y = grid$y,
    dens = dens,
    log_shape = center_log_shape(dens, grid),
    stringsAsFactors = FALSE
  )
}

estimate_parametric_background_grid <- function(fit, grid, bg_case) {
  
  if (is.null(fit) || is.null(fit$model.bg) || is.null(fit$model.bg$mod_global)) {
    return(NULL)
  }
  
  mod <- fit$model.bg$mod_global
  
  if (bg_case == "mark_cat") {
    
    grid_A <- grid
    grid_B <- grid
    
    grid_A$mark_cat <- factor("A", levels = c("A", "B"))
    grid_B$mark_cat <- factor("B", levels = c("A", "B"))
    
    eta_A <- tryCatch(
      as.numeric(predict(mod, newdata = grid_A, type = "link")),
      error = function(e) rep(NA_real_, nrow(grid))
    )
    
    eta_B <- tryCatch(
      as.numeric(predict(mod, newdata = grid_B, type = "link")),
      error = function(e) rep(NA_real_, nrow(grid))
    )
    
    ## Ground/background spatial shape after integrating over the mark space.
    ## Since mark_cat has no spatial covariate, this should be spatially flat.
    lambda_shape <- exp(eta_A) + exp(eta_B)
    
  } else {
    
    eta <- tryCatch(
      as.numeric(predict(mod, newdata = grid, type = "link")),
      error = function(e) rep(NA_real_, nrow(grid))
    )
    
    lambda_shape <- exp(eta)
  }
  
  dens <- normalize_density_on_grid(lambda_shape, grid)
  
  data.frame(
    x = grid$x,
    y = grid$y,
    dens = dens,
    log_shape = center_log_shape(dens, grid),
    stringsAsFactors = FALSE
  )
}

true_background_grid_df <- function(grid) {
  
  dens <- grid$dens_bg_true_spatial
  
  data.frame(
    x = grid$x,
    y = grid$y,
    dens = dens,
    log_shape = grid$eta_bg_true_spatial -
      weighted.mean(grid$eta_bg_true_spatial, w = grid$w, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}



select_representative_main_files <- function(
    count_case_focus = "N1000_balanced_50_50",
    model_for_choice = "param_correct",
    metric_for_choice = "auc_bg"
) {
  
  candidate_metric <- class_metrics_tbl %>%
    filter(
      count_case == count_case_focus,
      model == model_for_choice,
      metric == metric_for_choice,
      is.finite(value)
    )
  
  if (!nrow(candidate_metric)) {
    stop("No classification metric found for representative selection.")
  }
  
  chosen <- candidate_metric %>%
    group_by(bg_case) %>%
    mutate(
      median_value = median(value, na.rm = TRUE),
      dist_to_median = abs(value - median_value)
    ) %>%
    slice_min(dist_to_median, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(bg_case, count_case, rep, value, median_value)
  
  chosen <- chosen %>%
    left_join(
      index_tbl %>% select(bg_case, count_case, rep, file),
      by = c("bg_case", "count_case", "rep")
    )
  
  chosen
}

representative_main <- select_representative_main_files(
  count_case_focus = "N1000_balanced_50_50",
  model_for_choice = "param_correct",
  metric_for_choice = "auc_bg"
)

representative_main


build_shape_recovery_dataset <- function(
    representative_tbl,
    grid_n_plot = GRID_N
) {
  
  out_list <- list()
  
  for (ii in seq_len(nrow(representative_tbl))) {
    
    bg_i <- representative_tbl$bg_case[ii]
    file_i <- representative_tbl$file[ii]
    
    obj_i <- readRDS(file_i)
    
    grid_i <- obj_i$grid_truth %||% make_grid_eval(bg_i, grid_n = grid_n_plot)
    
    fit_classic_i <- get_fit_value(obj_i, "classic")
    fit_param_i <- get_fit_value(obj_i, "param_correct")
    
    true_i <- true_background_grid_df(grid_i)
    classic_i <- estimate_classic_background_grid(fit_classic_i, grid_i)
    param_i <- estimate_parametric_background_grid(fit_param_i, grid_i, bg_case = bg_i)
    
    true_i$model_surface <- "True"
    classic_i$model_surface <- "Classic"
    param_i$model_surface <- "Parametric"
    
    tmp <- bind_rows(true_i, classic_i, param_i) %>%
      mutate(
        bg_case = bg_i,
        count_case = representative_tbl$count_case[ii],
        rep = representative_tbl$rep[ii]
      )
    
    out_list[[ii]] <- tmp
  }
  
  bind_rows(out_list)
}

shape_recovery_df <- build_shape_recovery_dataset(
  representative_tbl = representative_main,
  grid_n_plot = GRID_N
)

shape_recovery_df <- shape_recovery_df %>%
  mutate(
    bg_case = factor(
      bg_case,
      levels = c("constant", "linear_xy", "cov1", "mark_cat", "hotspot"),
      labels = c("Constant", "Linear x/y", "Covariate Z1", "Categorical mark", "Hotspot")
    ),
    model_surface = factor(
      model_surface,
      levels = c("True", "Classic", "Parametric")
    )
  )


## clipping per evitare che pochi picchi della KDE classica dominino la scala
lim_shape <- quantile(abs(shape_recovery_df$log_shape), 0.98, na.rm = TRUE)

if (!is.finite(lim_shape) || lim_shape <= 0) {
  lim_shape <- max(abs(shape_recovery_df$log_shape), na.rm = TRUE)
}

shape_recovery_df <- shape_recovery_df %>%
  mutate(
    log_shape_clip = pmax(pmin(log_shape, lim_shape), -lim_shape)
  )

p_shape_recovery <- ggplot(
  shape_recovery_df,
  aes(x = x, y = y, fill = log_shape_clip)
) +
  geom_raster() +
  coord_equal() +
  facet_grid(
    model_surface ~ bg_case
  ) +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-lim_shape, lim_shape),
    name = "Centered\nlog-shape"
  ) +
  labs(
    title = "Background shape recovery",
    subtitle = "One representative balanced 50/50 replication per background scenario",
    x = "x coordinate",
    y = "y coordinate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_shape_recovery

