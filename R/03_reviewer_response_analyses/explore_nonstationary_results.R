#######################################################################
# ESPLORAZIONE DEI RISULTATI: ROBUSTEZZA ALLA NON-STAZIONARIETA' TEMPORALE
#
# Da eseguire DOPO aver completato NONSTAT_MODE = "run".
# Non rilancia nessuna simulazione/fit: legge solo gli .rds/.csv gia'
# prodotti e costruisce tabelle, grafici e un riassunto testuale
# orientato a rispondere ai punti 1 e 2 di Reviewer 1.
#
# USO:
#   setwd("<cartella con gli script originali e questo file>")
#   source("explore_nonstationary_results.R")
#
# Se la sessione corrente ha gia' in ambiente gli oggetti della pipeline
# (perche' hai appena eseguito NONSTAT_MODE = "run" nella stessa sessione)
# questo script li riusa direttamente. Altrimenti ricarica in modo sicuro
# (nessun autorun) lo script principale della pipeline.
#######################################################################

## ================================================================
## 0. CARICAMENTO SICURO DELL'AMBIENTE (se non gia' presente)
## ================================================================

if (!exists("collect_nonstat_results", mode = "function", inherits = TRUE)) {
  cat("[INFO] Ambiente non trovato in sessione: ricarico la pipeline in modalita' 'none'.\n")
  Sys.setenv(NONSTAT_MODE = "none")
  main_script <- Sys.getenv("NONSTAT_MAIN_SCRIPT", unset = "etas_nonstationary_background_robustness_check.R")
  if (!file.exists(main_script)) {
    stop("Non trovo '", main_script, "'. Imposta Sys.setenv(NONSTAT_MAIN_SCRIPT = '/percorso/completo/file.R') ",
         "oppure esegui questo script dalla stessa cartella.")
  }
  source(main_script, local = .GlobalEnv, echo = FALSE)
} else {
  cat("[INFO] Ambiente gia' presente in sessione: riuso le funzioni caricate.\n")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

options(width = 140)


## ================================================================
## 1. CARICAMENTO DEI RISULTATI
## ================================================================

res <- collect_nonstat_results()

cat("\n============================================================\n")
cat("RIEPILOGO DIMENSIONI DEI DATI RACCOLTI\n")
cat("============================================================\n")
for (nm in names(res)) {
  cat(sprintf("  %-16s : %5d righe\n", nm, nrow(res[[nm]])))
}


## ================================================================
## 2. STEP 1 -- STATO DEI FIT: SUCCESSO E STABILITA'
## ================================================================
## Da controllare SEMPRE per primo: se un modello fallisce in modo non
## casuale, il confronto finale e' su un campione distorto (sopravvivono
## solo i cataloghi "facili").

cat("\n============================================================\n")
cat("STEP 1 -- TASSO DI SUCCESSO DEI FIT (per modello e scenario)\n")
cat("============================================================\n")
success_tbl <- summarise_fit_success_rate(res)
print(success_tbl, row.names = FALSE)

if (any(success_tbl$success_rate < 1)) {
  cat("\n[ATTENZIONE] Almeno un modello/scenario ha un tasso di successo < 1.\n",
      "Le tabelle di bias sottostanti si basano SOLO sulle repliche riuscite:\n",
      "verificare se le repliche fallite sono associate a caratteristiche\n",
      "sistematiche (es. cataloghi piu' piccoli, piu' vicini al bordo della finestra).\n")
} else {
  cat("\n[OK] Tutti i modelli convergono su tutte le repliche in entrambi gli scenari.\n")
}

## Instabilita' numerica silenziosa: un fit puo' avere ok=TRUE ma essere
## degenere (es. k0 -> 0 o enorme, std_error mancante). Flag esplicito.
flag_unstable_trigger_fits <- function(trigger_bias_df, extreme_rel_error = 5) {
  d <- trigger_bias_df
  if (!nrow(d)) return(d)
  bad_by_row <- is.na(d$std_error) | (is.finite(d$rel_error) & abs(d$rel_error) > extreme_rel_error)
  unstable_ids <- unique(d[bad_by_row, c("model", "bg_case", "rep")])
  d$fit_key <- paste(d$model, d$bg_case, d$rep, sep = "__")
  unstable_key <- paste(unstable_ids$model, unstable_ids$bg_case, unstable_ids$rep, sep = "__")
  d$is_unstable <- d$fit_key %in% unstable_key
  d$fit_key <- NULL
  d
}

trig_flagged <- flag_unstable_trigger_fits(res$trigger_bias, extreme_rel_error = 5)

instab_summary <- aggregate(
  is_unstable ~ model + bg_case, data = unique(trig_flagged[, c("model", "bg_case", "rep", "is_unstable")]),
  FUN = function(x) mean(x, na.rm = TRUE)
)
names(instab_summary)[names(instab_summary) == "is_unstable"] <- "frac_repliche_instabili"

cat("\n--- Frazione di repliche con fit di triggering numericamente instabile ---\n")
cat("(std_error mancante su almeno un parametro, oppure |errore relativo| > 500%)\n")
print(instab_summary, row.names = FALSE)


## ================================================================
## 3. STEP 2 -- RECUPERO DELLA FORMA SPAZIALE DEL BACKGROUND (rMISEbg)
## ================================================================
## Risponde a: "il drift temporale finisce nel background spaziale?"
## Atteso, per costruzione (separabilita' mu(t)*f0(s)): NESSUNA differenza
## sistematica rispetto ai valori del piano stazionario originale.

cat("\n============================================================\n")
cat("STEP 2 -- RECUPERO DELLA FORMA SPAZIALE DEL BACKGROUND (rMISEbg)\n")
cat("============================================================\n")

if (nrow(res$bg_shape)) {
  bg_shape_summary <- do.call(rbind, lapply(split(res$bg_shape, list(res$bg_shape$bg_case, res$bg_shape$model)), function(sub) {
    if (!nrow(sub)) return(NULL)
    x <- sub$rmise_bg[is.finite(sub$rmise_bg)]
    data.frame(
      bg_case = sub$bg_case[1], model = sub$model[1], n = length(x),
      median = if (length(x)) median(x) else NA,
      q1 = if (length(x)) quantile(x, 0.25) else NA,
      q3 = if (length(x)) quantile(x, 0.75) else NA,
      n_na = sum(!is.finite(sub$rmise_bg))
    )
  }))
  bg_shape_summary <- bg_shape_summary[order(bg_shape_summary$bg_case, bg_shape_summary$model), ]
  rownames(bg_shape_summary) <- NULL
  print(bg_shape_summary, row.names = FALSE, digits = 4)

  cat("\nInterpretazione attesa: mediane rMISEbg comparabili, scenario per scenario,\n",
      "a quelle riportate nelle Tabelle S1/S3 del supplementary per il caso\n",
      "STAZIONARIO corrispondente (constant/cov1, regime Balanced N=1000).\n",
      "Se ETAS-P mostra un peggioramento marcato rispetto al caso stazionario,\n",
      "significa che il drift temporale sta distorcendo anche la stima spaziale\n",
      "(esito NON atteso sotto separabilita', da indagare ulteriormente).\n")
} else {
  cat("[ATTENZIONE] Nessun dato di rMISEbg disponibile.\n")
}

p_bg_shape <- plot_bg_shape_recovery(res)
if (!is.null(p_bg_shape)) print(p_bg_shape)


## ================================================================
## 4. STEP 3 -- BIAS DEI PARAMETRI DI TRIGGERING
## ================================================================
## Risponde a: "il drift finisce nel triggering (k0, c, p, d, q)?"
## Riportiamo il riassunto SIA su tutte le repliche SIA escludendo quelle
## flaggate come instabili (Step 1), per verificare la robustezza della
## conclusione a possibili ottimi degeneri isolati.

cat("\n============================================================\n")
cat("STEP 3 -- BIAS RELATIVO DEI PARAMETRI DI TRIGGERING\n")
cat("============================================================\n")

summarise_trigger_bias <- function(d) {
  d <- d[is.finite(d$rel_error), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  out <- do.call(rbind, lapply(split(d, list(d$bg_case, d$model, d$parameter)), function(sub) {
    if (!nrow(sub)) return(NULL)
    data.frame(
      bg_case = sub$bg_case[1], model = sub$model[1], parameter = sub$parameter[1],
      n = nrow(sub),
      median_rel_error = median(sub$rel_error),
      q1 = quantile(sub$rel_error, 0.25), q3 = quantile(sub$rel_error, 0.75)
    )
  }))
  out <- out[order(out$bg_case, out$parameter, out$model), ]
  rownames(out) <- NULL
  out
}

cat("\n--- (a) Tutte le repliche ---\n")
trig_summary_all <- summarise_trigger_bias(res$trigger_bias)
print(trig_summary_all, row.names = FALSE, digits = 3)

cat("\n--- (b) Escludendo le repliche flaggate come instabili (Step 1) ---\n")
trig_summary_stable <- summarise_trigger_bias(trig_flagged[!trig_flagged$is_unstable, ])
print(trig_summary_stable, row.names = FALSE, digits = 3)

cat("\nConfronta (a) e (b): se le mediane non cambiano in modo sostanziale,\n",
    "la conclusione e' robusta alla presenza di eventuali ottimi degeneri isolati.\n",
    "Se cambiano molto, il bias osservato in (a) e' guidato da poche repliche\n",
    "patologiche e va discusso separatamente dal pattern sistematico.\n")

p_trig <- plot_trigger_bias(res)
if (!is.null(p_trig)) print(p_trig)

## Versione robusta del grafico (solo repliche stabili), utile per la
## lettera di risposta se emergono outlier estremi nel grafico completo.
plot_trigger_bias_stable <- function(trig_flagged_df) {
  d <- trig_flagged_df[!trig_flagged_df$is_unstable & is.finite(trig_flagged_df$rel_error), ]
  if (!nrow(d)) return(NULL)
  d$model <- factor(d$model, levels = c("classic", "flp", "param"))
  d$bg_case <- factor(d$bg_case, levels = names(BGCASE_LABELS), labels = BGCASE_LABELS)
  ggplot(d, aes(x = parameter, y = rel_error, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
    facet_wrap(~ bg_case, nrow = 2) +
    scale_fill_manual(values = unname(MODEL_COLORS), labels = MODEL_LABELS, name = NULL) +
    labs(title = "Bias relativo dei parametri ETAS (repliche numericamente stabili)",
        subtitle = "Escluse le repliche con std_error mancante o |errore relativo| > 500% su almeno un parametro",
        x = NULL, y = "Errore relativo (stima - vero) / vero") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))
}

p_trig_stable <- plot_trigger_bias_stable(trig_flagged)
if (!is.null(p_trig_stable)) print(p_trig_stable)


## ================================================================
## 5. STEP 4 -- BIAS DEL COEFFICIENTE SPAZIALE Z1 (scenario cov1, solo ETAS-P)
## ================================================================
## Analogo diretto del coefficiente di distanza dalla faglia nell'applicazione
## reale: se il drift temporale "inquina" il coefficiente della covariata
## spaziale, e' il segnale piu' rilevante per l'interpretazione del catalogo
## italiano.

cat("\n============================================================\n")
cat("STEP 4 -- BIAS DEL COEFFICIENTE SPAZIALE Z1 (solo ETAS-P, scenario cov1)\n")
cat("============================================================\n")

if (nrow(res$bg_coef_bias)) {
  coef_summary <- do.call(rbind, lapply(split(res$bg_coef_bias, res$bg_coef_bias$coefficient), function(sub) {
    data.frame(
      coefficient = sub$coefficient[1], n = nrow(sub), true_value = sub$true[1],
      median_estimate = median(sub$estimate, na.rm = TRUE),
      median_error = median(sub$error, na.rm = TRUE),
      q1_error = quantile(sub$error, 0.25, na.rm = TRUE),
      q3_error = quantile(sub$error, 0.75, na.rm = TRUE)
    )
  }))
  print(coef_summary, row.names = FALSE, digits = 4)
} else {
  cat("[ATTENZIONE] Nessun dato per il coefficiente Z1 (verifica che bg_case=='cov1' abbia fit ETAS-P riusciti).\n")
}

p_coef <- plot_bg_coef_bias(res)
if (!is.null(p_coef)) print(p_coef)


## ================================================================
## 6. STEP 5 -- ACCURATEZZA DI CLASSIFICAZIONE BACKGROUND/TRIGGERED
## ================================================================

cat("\n============================================================\n")
cat("STEP 5 -- ACCURATEZZA DI CLASSIFICAZIONE (soglia 0.5 su rho_i)\n")
cat("============================================================\n")

if (nrow(res$classification) && all(c("metrics.metric", "metrics.value") %in% names(res$classification))) {
  acc_df <- res$classification[res$classification$metrics.metric == "accuracy_05", ]
  classif_summary <- aggregate(metrics.value ~ bg_case + model, data = acc_df,
                               FUN = function(x) median(x, na.rm = TRUE))
  names(classif_summary)[names(classif_summary) == "metrics.value"] <- "median_accuracy"
  print(classif_summary, row.names = FALSE)
  cat("\nAtteso: valori alti (~0.98+) e simili tra i tre modelli, come nel caso\n",
      "stazionario -- la classificazione binaria e' poco sensibile alla\n",
      "non-stazionarieta' temporale, che invece agisce sulla FORMA della\n",
      "densita' di background stimata (vedi Step 2, 3, 4).\n",
      "\nNOTA METODOLOGICA: un'alta accuratezza di classificazione, anche se\n",
      "uniforme nel tempo, NON implica che il modello stimato 'conosca' il\n",
      "trend di mu(t) -- implica solo che la classificazione si basa su\n",
      "evidenza spaziale/storica locale, robusta indipendentemente dal drift.\n",
      "Per questo la diagnostica temporale basata su rho_i aggregate per bin\n",
      "(Step 6) NON e' una prova di recupero del trend: si veda la nota li'.\n")
} else {
  cat("[NOTA] Struttura di res$classification diversa dal previsto:\n")
  if (nrow(res$classification)) print(head(res$classification))
}


## ================================================================
## 7. STEP 6 -- DIAGNOSTICA TEMPORALE: LA VERIFICA PIU' DIRETTA
## ================================================================
## Grafico gia' prodotto da build_all_nonstat_plots(): mu(t) vera vs massa
## di probabilita' posteriore di background stimata per bin, mediata sulle
## repliche. Qui aggiungiamo una sintesi NUMERICA (pendenza e correlazione)
## per quantificare quanto ciascun modello "vede" il trend, invece di
## affidarci solo all'ispezione visiva.

cat("\n============================================================\n")
cat("STEP 6 -- DIAGNOSTICA TEMPORALE: LETTURA CORRETTA\n")
cat("============================================================\n")
cat("\n[NOTA METODOLOGICA IMPORTANTE]\n",
    "Questa diagnostica somma le probabilita' posteriori rho_i per bin\n",
    "temporale. Dato che l'accuratezza di classificazione e' molto alta\n",
    "(vedi Step 5), la massa per bin finisce per riprodurre quasi\n",
    "meccanicamente il CONTEGGIO VERO di eventi di background realmente\n",
    "occorsi in quel bin (che varia nel tempo per costruzione del disegno\n",
    "sperimentale), indipendentemente da cosa 'creda' il parametro mu\n",
    "stimato dal modello (che resta uno scalare costante in tutti e tre i\n",
    "casi). Una correlazione alta e uniforme tra i tre modelli NON e'\n",
    "quindi prova che il modello 'recupera' il trend: e' prova che la\n",
    "CLASSIFICAZIONE binaria resta accurata nel tempo, cosa gia' mostrata\n",
    "dallo Step 5. La domanda 'dove va a finire il segnale nonstazionario'\n",
    "ha gia' una risposta diretta e corretta negli Step 2-4 (forma spaziale,\n",
    "parametri di triggering, coefficiente della covariata): questo Step 6\n",
    "va quindi letto come controllo di robustezza della classificazione\n",
    "sotto drift temporale, non come evidenza di recupero del trend.\n\n")

p_time <- plot_time_binned_recovery(res)
if (!is.null(p_time)) print(p_time)

## Per ogni replica e modello: regressione lineare fitted_mass_frac ~ bin,
## per estrarre una pendenza stimata confrontabile con la pendenza vera
## (identica per costruzione in ogni replica, dato bg_case/kappa fissati).
compute_temporal_recovery_slope <- function(time_binned_df) {
  d <- time_binned_df
  if (!nrow(d)) return(data.frame())
  out <- do.call(rbind, lapply(split(d, list(d$bg_case, d$model, d$rep)), function(sub) {
    if (!nrow(sub) || length(unique(sub$bin)) < 3) return(NULL)
    fit_hat <- try(lm(fitted_mass_frac ~ bin, data = sub), silent = TRUE)
    fit_true <- try(lm(true_mass_frac ~ bin, data = sub), silent = TRUE)
    if (inherits(fit_hat, "try-error") || inherits(fit_true, "try-error")) return(NULL)
    data.frame(
      bg_case = sub$bg_case[1], model = sub$model[1], rep = sub$rep[1],
      slope_fitted = unname(coef(fit_hat)[2]),
      slope_true = unname(coef(fit_true)[2]),
      cor_fitted_true = suppressWarnings(cor(sub$fitted_mass_frac, sub$true_mass_frac))
    )
  }))
  out
}

slope_tbl <- compute_temporal_recovery_slope(res$time_binned)

if (nrow(slope_tbl)) {
  cat("\n--- Pendenza stimata (fitted_mass_frac ~ bin) vs pendenza vera, per modello ---\n")
  slope_summary <- do.call(rbind, lapply(split(slope_tbl, list(slope_tbl$bg_case, slope_tbl$model)), function(sub) {
    if (!nrow(sub)) return(NULL)
    data.frame(
      bg_case = sub$bg_case[1], model = sub$model[1], n = nrow(sub),
      slope_true = sub$slope_true[1],
      median_slope_fitted = median(sub$slope_fitted, na.rm = TRUE),
      q1_slope_fitted = quantile(sub$slope_fitted, 0.25, na.rm = TRUE),
      q3_slope_fitted = quantile(sub$slope_fitted, 0.75, na.rm = TRUE),
      median_cor = median(sub$cor_fitted_true, na.rm = TRUE)
    )
  }))
  slope_summary <- slope_summary[order(slope_summary$bg_case, slope_summary$model), ]
  rownames(slope_summary) <- NULL
  print(slope_summary, row.names = FALSE, digits = 3)

  cat("\nLettura CORRETTA (aggiornata dopo i risultati osservati):\n",
      "- 'slope_true' e' la stessa per costruzione in ogni replica (dipende solo da kappa/n_bin).\n",
      "- Se 'median_slope_fitted' e 'median_cor' sono ALTI e SOSTANZIALMENTE\n",
      "  IDENTICI tra ETAS, ETAS-FLP ed ETAS-P (come tipicamente osservato),\n",
      "  questo NON significa che i modelli 'conoscano' il trend di mu(t):\n",
      "  significa che la massa di probabilita' posteriore di background per\n",
      "  bin riproduce il conteggio vero di eventi di background realmente\n",
      "  occorsi in quel bin, un artefatto della buona accuratezza di\n",
      "  classificazione (Step 5) che e' insensibile a come i tre modelli\n",
      "  parametrizzano mu (tutti con un singolo scalare costante).\n",
      "- Una differenza sistematica TRA i tre modelli sarebbe stata\n",
      "  informativa; l'assenza di differenza indica solo che questa specifica\n",
      "  diagnostica non discrimina tra modelli su questo aspetto.\n",
      "- La risposta a 'il modello vede il trend?' va cercata nella stima dei\n",
      "  PARAMETRI (Step 3) e della FORMA spaziale (Step 2), non nella massa\n",
      "  posteriore per bin: e' li' che i tre modelli, per costruzione\n",
      "  incapaci di rappresentare mu(t), possono davvero differire.\n")

  ## Boxplot delle pendenze stimate per modello, con riferimento alla pendenza vera
  d_plot <- slope_tbl
  d_plot$model <- factor(d_plot$model, levels = c("classic", "flp", "param"))
  d_plot$bg_case <- factor(d_plot$bg_case, levels = names(BGCASE_LABELS), labels = BGCASE_LABELS)
  true_ref <- unique(d_plot[, c("bg_case", "slope_true")])

  p_slope <- ggplot(d_plot, aes(x = model, y = slope_fitted, fill = model)) +
    geom_hline(data = true_ref, aes(yintercept = slope_true),
              linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_boxplot(outlier.alpha = 0.4) +
    facet_wrap(~ bg_case, nrow = 1) +
    scale_x_discrete(labels = MODEL_LABELS) +
    scale_fill_manual(values = unname(MODEL_COLORS), labels = MODEL_LABELS, guide = "none") +
    labs(title = "Pendenza della massa di background declusterizzata per bin temporale",
        subtitle = "Riflette l'accuratezza di classificazione (Step 5), non la conoscenza del trend da parte del modello stimato; vedi nota metodologica",
        x = NULL, y = "Pendenza (fitted_mass_frac ~ bin)") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
  print(p_slope)

  ggsave(file.path(NONSTAT_OUT_ROOT, "figures", "fig_E_temporal_slope_recovery.png"),
        p_slope, width = 9, height = 5, dpi = 300)
} else {
  cat("[ATTENZIONE] Nessun dato per la diagnostica temporale.\n")
}


## ================================================================
## 8. VISUALIZZAZIONE DEI FILE PNG GIA' SALVATI (controllo rapido)
## ================================================================

cat("\n============================================================\n")
cat("FILE GIA' SALVATI IN: ", file.path(NONSTAT_OUT_ROOT, "figures"), "\n")
cat("============================================================\n")
print(list.files(file.path(NONSTAT_OUT_ROOT, "figures"), full.names = FALSE))


## ================================================================
## 9. RIASSUNTO TESTUALE AUTOMATICO (bozza per la lettera di risposta)
## ================================================================

cat("\n============================================================\n")
cat("BOZZA DI RIASSUNTO PER LA LETTERA DI RISPOSTA AL REVISORE\n")
cat("============================================================\n\n")

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

for (bgc in unique(res$bg_shape$bg_case)) {

  cat(sprintf("--- Scenario: %s ---\n", bgc))

  bgs <- bg_shape_summary[bg_shape_summary$bg_case == bgc, ]
  if (nrow(bgs)) {
    cat(sprintf("  rMISEbg mediano: ETAS=%.4f, ETAS-FLP=%.4f, ETAS-P=%.4f\n",
                bgs$median[bgs$model == "classic"] %||% NA,
                bgs$median[bgs$model == "flp"] %||% NA,
                bgs$median[bgs$model == "param"] %||% NA))
  }

  ts <- trig_summary_stable[trig_summary_stable$bg_case == bgc & trig_summary_stable$model == "param", ]
  if (nrow(ts)) {
    worst <- ts[which.max(abs(ts$median_rel_error)), ]
    cat(sprintf("  ETAS-P, parametro di triggering con bias mediano piu' marcato: %s (%s)\n",
                worst$parameter, fmt_pct(worst$median_rel_error)))
  }

  ss <- slope_summary[slope_summary$bg_case == bgc, ]
  if (nrow(ss)) {
    for (i in seq_len(nrow(ss))) {
      cat(sprintf("  Modello %-8s: pendenza vera=%.4f, pendenza stimata (mediana)=%.4f, correlazione mediana=%.2f\n",
                  ss$model[i], ss$slope_true[i], ss$median_slope_fitted[i], ss$median_cor[i]))
    }
  }
  cat("\n")
}

cat("Ricorda di riportare anche:\n",
    "- il tasso di successo/instabilita' dei fit (Step 1), soprattutto se\n",
    "  ETAS-FLP ha richiesto retry sistematici sotto il trend, e la frazione\n",
    "  di repliche escluse come numericamente instabili nello Step 3;\n",
    "- il confronto DIRETTO tra i valori di rMISEbg e del bias su mu ottenuti\n",
    "  qui e i valori corrispondenti gia' riportati nelle vostre Tabelle S1/S3\n",
    "  e nella Sez. 4.4 per il caso STAZIONARIO: se sono comparabili (come\n",
    "  atteso qui), e' l'argomento piu' diretto per mostrare che il trend\n",
    "  temporale non introduce distorsioni NUOVE oltre a quelle gia' note e\n",
    "  documentate;\n",
    "- il bias del coefficiente Z1 (Step 4) come analogo diretto del\n",
    "  coefficiente di distanza dalla faglia nell'applicazione reale;\n",
    "- EVITARE di presentare la diagnostica temporale (Step 6) come prova che\n",
    "  i modelli 'recuperano' il trend: la sua lettura corretta e' che la\n",
    "  classificazione background/triggered resta robusta sotto drift\n",
    "  temporale non modellato, un punto di supporto secondario ma NON una\n",
    "  risposta alla domanda su dove finisce il segnale nonstazionario\n",
    "  (quella risposta sta negli Step 2-4).\n")

cat("\n[FINE ESPLORAZIONE]\n")
