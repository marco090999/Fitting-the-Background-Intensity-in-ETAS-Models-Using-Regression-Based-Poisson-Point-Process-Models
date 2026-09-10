#######################################################################
# ESPLORAZIONE DEI RISULTATI: COPERTURA NAIVE vs BOOTSTRAP
# (Reviewer 1, punto 3 -- errori standard/CI non condizionati)
#
# Da eseguire DOPO aver completato COV_MODE = "run" (versione
# sequenziale). Non rilancia nessun fit: legge solo i risultati gia'
# salvati e produce tabelle, diagnostiche e grafici pronti da
# condividere.
#
# USO:
#   setwd("<cartella con etas_coverage_bootstrap.R>")
#   source("explore_coverage_bootstrap_results.R")
#######################################################################

## ================================================================
## 0. CARICAMENTO SICURO DELL'AMBIENTE (se non gia' presente)
## ================================================================

if (!exists("collect_coverage_results", mode = "function", inherits = TRUE)) {
  cat("[INFO] Ambiente non trovato in sessione: ricarico in modalita' 'none'.\n")
  Sys.setenv(COV_MODE = "none")
  cov_script <- Sys.getenv("COV_SCRIPT_PATH", unset = "etas_coverage_bootstrap.R")
  if (!file.exists(cov_script)) {
    stop("Non trovo '", cov_script, "'. Imposta Sys.setenv(COV_SCRIPT_PATH = '/percorso/completo/file.R') ",
         "oppure esegui questo script dalla stessa cartella.")
  }
  source(cov_script, local = .GlobalEnv, echo = FALSE)
} else {
  cat("[INFO] Ambiente gia' presente in sessione: riuso le funzioni caricate.\n")
}

suppressPackageStartupMessages(library(ggplot2))
options(width = 140)


## ================================================================
## 1. STATO DI AVANZAMENTO: QUANTI JOB SONO STATI EFFETTIVAMENTE FATTI
## ================================================================
## Prima di guardare qualunque numero di copertura, verifica quante
## delle repliche esterne attese (fino a COV_M_OUTER per scenario) sono
## effettivamente presenti -- alcune potrebbero essere state scartate
## dal filtro di plausibilita' su theta_hat (c,d,q entro 5x dal vero
## valore), quindi il denominatore reale puo' essere minore di 30.

boot_files <- list.files(file.path(COV_OUT_ROOT, "bootstrap_fits"), pattern = "\\.rds$", full.names = TRUE)
cat("\n============================================================\n")
cat("STATO DI AVANZAMENTO\n")
cat("============================================================\n")
cat(sprintf("File di bootstrap trovati: %d\n", length(boot_files)))

job_status <- do.call(rbind, lapply(boot_files, function(f) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  boot_res <- obj$bootstrap
  n_boot_ok_vec <- if (!is.null(boot_res)) boot_res$n_boot_ok else NA
  data.frame(
    scenario = obj$metadata$scenario, outer_rep = obj$metadata$outer_rep,
    B_inner = obj$metadata$B_inner %||% NA,
    min_n_boot_ok = if (length(n_boot_ok_vec)) min(n_boot_ok_vec, na.rm = TRUE) else NA,
    max_n_boot_ok = if (length(n_boot_ok_vec)) max(n_boot_ok_vec, na.rm = TRUE) else NA,
    stringsAsFactors = FALSE
  )
}))

print(table(job_status$scenario))
cat("\n--- Distribuzione di n_boot_ok (fit bootstrap interni riusciti su B_inner) ---\n")
print(summary(job_status[, c("min_n_boot_ok", "max_n_boot_ok")]))

low_success <- job_status[job_status$min_n_boot_ok < 50, ]
if (nrow(low_success)) {
  cat("\n[ATTENZIONE] Repliche esterne con meno del 50% di fit bootstrap interni riusciti:\n")
  print(low_success, row.names = FALSE)
  cat("Per queste, la copertura bootstrap e' stimata su un campione ridotto -- interpretare con cautela.\n")
} else {
  cat("\n[OK] Tutte le repliche hanno almeno il 50% di fit bootstrap interni riusciti.\n")
}


## ================================================================
## 2. TABELLA PRINCIPALE: COPERTURA NAIVE vs BOOTSTRAP (Tabella S-style)
## ================================================================
## naive_coverage: ricalcolata su TUTTE le repliche disponibili nel main
## plan (non solo il sottoinsieme M_outer) -- il numero piu' solido.
## boot_coverage: sul sottoinsieme M_outer effettivamente processato.

cat("\n============================================================\n")
cat("TABELLA PRINCIPALE: COPERTURA NAIVE (CONDIZIONALE) vs BOOTSTRAP\n")
cat("============================================================\n")

tab_s12 <- build_table_S12()
print(tab_s12, row.names = FALSE, digits = 4)

## Intervallo di confidenza binomiale esatto sulla copertura bootstrap,
## dato il campione ridotto (M_outer tipicamente 30): serve a capire
## quanto ci si puo' fidare del singolo numero di copertura osservato.
cat("\n--- Intervallo di confidenza binomiale (Clopper-Pearson, 95%) sulla copertura bootstrap ---\n")
cat("(quantifica l'incertezza Monte Carlo dovuta al numero ridotto di repliche esterne)\n\n")
for (i in seq_len(nrow(tab_s12))) {
  row <- tab_s12[i, ]
  if (is.na(row$n_boot) || row$n_boot == 0) next
  n_covered <- round(row$boot_coverage * row$n_boot)
  ci <- binom.test(n_covered, row$n_boot)$conf.int
  cat(sprintf("%-12s %-8s: copertura=%.3f su n=%d  ->  IC 95%%: [%.3f, %.3f]\n",
             row$scenario, row$coefficient, row$boot_coverage, row$n_boot, ci[1], ci[2]))
}
cat("\nSe l'intervallo copre 0.95, il risultato e' compatibile con una copertura nominale corretta.\n",
    "Intervalli larghi riflettono semplicemente la scarsa numerosita' campionaria (M_outer),\n",
    "non un difetto del metodo.\n")


## ================================================================
## 3. CONFRONTO DELLE AMPIEZZE DEGLI INTERVALLI (naive vs bootstrap)
## ================================================================
## Il confronto piu' diretto per rispondere al revisore: il CI
## bootstrap e' sistematicamente piu' ampio di quello condizionale?

cat("\n============================================================\n")
cat("AMPIEZZA MEDIA DEGLI INTERVALLI: NAIVE vs BOOTSTRAP\n")
cat("============================================================\n")

width_summary <- tab_s12[, c("scenario", "coefficient", "naive_mean_width", "boot_mean_width")]
width_summary$ratio_boot_over_naive <- width_summary$boot_mean_width / width_summary$naive_mean_width
print(width_summary, row.names = FALSE, digits = 3)
cat("\n'ratio_boot_over_naive' > 1 indica che il CI bootstrap e' piu' ampio di quello\n",
    "condizionale -- coerente con l'ipotesi del revisore che il SE condizionale\n",
    "sottostimi l'incertezza reale, se il rapporto e' consistentemente > 1.\n")


## ================================================================
## 4. CONFRONTO APPAIATO PER REPLICA ESTERNA (naive vs bootstrap coverage)
## ================================================================
## Per ciascuna delle repliche esterne effettivamente bootstrappate,
## confronta se la copertura naive e quella bootstrap concordano o
## divergono sulla STESSA replica -- utile per capire se i fallimenti
## di copertura naive sono sistematici o isolati.

boot_res <- collect_coverage_results()

paired <- merge(
  boot_res$naive_subset[, c("scenario", "rep", "coefficient", "covered")],
  boot_res$bootstrap[, c("scenario", "outer_rep", "coefficient", "covered")],
  by.x = c("scenario", "rep", "coefficient"), by.y = c("scenario", "outer_rep", "coefficient"),
  suffixes = c("_naive", "_boot")
)

cat("\n============================================================\n")
cat("CONFRONTO APPAIATO: COPERTURA NAIVE vs BOOTSTRAP, STESSA REPLICA ESTERNA\n")
cat("============================================================\n")
if (nrow(paired)) {
  paired_summary <- aggregate(
    cbind(naive_covers = covered_naive, boot_covers = covered_boot) ~ scenario + coefficient,
    data = paired, FUN = function(x) mean(as.logical(x), na.rm = TRUE)
  )
  print(paired_summary, row.names = FALSE, digits = 3)

  ## casi in cui naive fallisce ma bootstrap copre correttamente (il
  ## pattern che ci si aspetta se il SE condizionale e' troppo stretto)
  paired$naive_fails_boot_ok <- (!paired$covered_naive) & paired$covered_boot
  n_pattern <- sum(paired$naive_fails_boot_ok, na.rm = TRUE)
  cat(sprintf("\nRepliche in cui il CI naive NON copre ma quello bootstrap SI': %d su %d confronti totali\n",
             n_pattern, nrow(paired)))
} else {
  cat("[NOTA] Nessuna riga in comune tra naive_subset e bootstrap: controlla la struttura di collect_coverage_results().\n")
}


## ================================================================
## 5. GRAFICI
## ================================================================

fig_dir <- file.path(COV_OUT_ROOT, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## 5a. Ampiezza degli intervalli: naive vs bootstrap, per coefficiente
if (nrow(boot_res$naive_subset) && nrow(boot_res$bootstrap)) {
  d_naive <- boot_res$naive_subset[, c("scenario", "coefficient", "width")]
  d_naive$type <- "Naive (condizionale)"
  d_boot <- boot_res$bootstrap
  d_boot$width <- d_boot$ci_hi - d_boot$ci_lo
  d_boot <- d_boot[, c("scenario", "coefficient", "width")]
  d_boot$type <- "Bootstrap"
  d_width <- rbind(d_naive, d_boot)
  d_width$label <- paste(d_width$scenario, d_width$coefficient, sep = " / ")

  p_width <- ggplot(d_width, aes(x = label, y = width, fill = type)) +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "Ampiezza degli intervalli di confidenza al 95%: naive vs bootstrap",
        subtitle = "Sullo stesso sottoinsieme di repliche esterne bootstrappate",
        x = NULL, y = "Ampiezza dell'intervallo", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  print(p_width)
  ggsave(file.path(fig_dir, "fig_ci_width_naive_vs_bootstrap.png"), p_width, width = 9, height = 5, dpi = 300)
}

## 5b. Copertura osservata (barre) con intervallo di confidenza
##     binomiale, naive (su tutte le repliche) vs bootstrap (su M_outer)
cov_plot_df <- do.call(rbind, lapply(seq_len(nrow(tab_s12)), function(i) {
  row <- tab_s12[i, ]
  out <- list()
  if (!is.na(row$n_naive) && row$n_naive > 0) {
    ci <- binom.test(round(row$naive_coverage * row$n_naive), row$n_naive)$conf.int
    out[[length(out) + 1]] <- data.frame(scenario = row$scenario, coefficient = row$coefficient,
                                         type = "Naive", coverage = row$naive_coverage,
                                         ci_lo = ci[1], ci_hi = ci[2], n = row$n_naive)
  }
  if (!is.na(row$n_boot) && row$n_boot > 0) {
    ci <- binom.test(round(row$boot_coverage * row$n_boot), row$n_boot)$conf.int
    out[[length(out) + 1]] <- data.frame(scenario = row$scenario, coefficient = row$coefficient,
                                         type = "Bootstrap", coverage = row$boot_coverage,
                                         ci_lo = ci[1], ci_hi = ci[2], n = row$n_boot)
  }
  do.call(rbind, out)
}))
cov_plot_df$label <- paste(cov_plot_df$scenario, cov_plot_df$coefficient, sep = " / ")

p_cov <- ggplot(cov_plot_df, aes(x = label, y = coverage, color = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.4), size = 0.6) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Copertura osservata degli intervalli al 95% nominale",
      subtitle = "Linea tratteggiata = copertura nominale (0.95); barre = IC binomiale di Clopper-Pearson",
      x = NULL, y = "Copertura osservata", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
print(p_cov)
ggsave(file.path(fig_dir, "fig_coverage_naive_vs_bootstrap.png"), p_cov, width = 9, height = 5, dpi = 300)


## ================================================================
## 6. SALVATAGGIO TABELLE
## ================================================================

utils::write.csv(tab_s12, file.path(COV_OUT_ROOT, "table_S12_final_summary.csv"), row.names = FALSE)
utils::write.csv(job_status, file.path(COV_OUT_ROOT, "table_job_status.csv"), row.names = FALSE)
if (exists("paired_summary")) {
  utils::write.csv(paired_summary, file.path(COV_OUT_ROOT, "table_paired_naive_vs_bootstrap.csv"), row.names = FALSE)
}

cat("\n[FATTO] Tabelle e grafici salvati in: ", normalizePath(COV_OUT_ROOT), "\n", sep = "")
cat("[FINE ESPLORAZIONE]\n")
