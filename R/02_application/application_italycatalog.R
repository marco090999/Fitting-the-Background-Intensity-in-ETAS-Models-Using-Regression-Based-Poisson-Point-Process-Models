

################################################
#### STUDIO DI SIMULAZIONE ETAS PARAMETRICO ####
################################################

# -----------------------------------------------------------------------------------------------------------------------------------
# APPLICAZIONI SUL CATALOGO ITALIA
# 1) catalogo italia: confronto stime parametri triggering tra etas parametrico ed etas classico
# 2) catalogo italia: inclusione covariate esterne, confronto tra stime di covariate nel triggering tra etas parametrico e classico
# 3) catalogo italia: inclusione covariate, confronto tra covariate esterne nel background e nel triggering
# 4) catalogo italia: inclusione covariate sia nel background che nel triggering
# -----------------------------------------------------------------------------------------------------------------------------------


library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(viridis)
library(etasFLP)

data("italycatalog")

# 1. Trasformo il catalogo in oggetto sf
italycatalog_sf <- italycatalog %>%
  filter(
    !is.na(long),
    !is.na(lat),
    !is.na(magn1)
  ) %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )

# 2. Scarico la mappa del mondo e dell'Italia
world <- ne_countries(
  scale = 50,
  returnclass = "sf"
)

italy <- ne_countries(
  country = "Italy",
  scale = 50,
  returnclass = "sf"
)

# 3. Bounding box intorno all'Italia
italy_bbox <- c(
  xmin = 5.5,
  xmax = 19.5,
  ymin = 35.0,
  ymax = 47.8
)

p_italy_eq_clean <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey94",
    color = "white",
    linewidth = 0.25
  ) +
  geom_sf(
    data = italy,
    fill = "grey88",
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_sf(
    data = italycatalog_sf,
    aes(color = magn1),
    size = 2.2,
    alpha = 0.78
  ) +
  coord_sf(
    xlim = c(5.5, 19.5),
    ylim = c(35.0, 47.8),
    expand = FALSE
  ) +
  scale_color_viridis_c(
    option = "magma",
    direction = -1,
    name = "Magnitude"
  ) +
  guides(
    color = guide_colorbar(
      barheight = grid::unit(10, "cm"),   # altezza colorbar
      barwidth  = grid::unit(0.7, "cm"),  # larghezza colorbar
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = "Italian earthquake catalogue",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 17) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 28,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text = element_text(
      size = 21,
      color = "grey20"
    ),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.25
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 22
    ),
    legend.text = element_text(
      size = 20
    ),
    plot.margin = margin(10, 20, 10, 10)
  )

p_italy_eq_clean




# 1)

# ETAS CLASSICO
etas.class.p1 <- etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9,
                   mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1,q = 1.5, betacov = 0.7,
                   params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), w = replicate(nrow(catalog.withcov), 1),
                   hvarx = replicate(nrow(catalog.withcov), 1), hvary = replicate(nrow(catalog.withcov),1),
                   formula1 = "time ~  magnitude- 1", declustering = TRUE,
                   thinning = FALSE, flp = FALSE, ndeclust = 15, onlytime = FALSE,
                   is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
                   compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

summary(etas.class.p1)


# ETAS PARAMETRICO
etas.par.p1 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
  tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
  mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
  params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
  offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
  hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
  ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
  sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
  formula.bg = ~ s(x, y, bs = "tp", k = 50), process.type.bg = "s2d", spatial.cov.bg = FALSE, mult.bg = 4, ncube.bg = NULL,
  verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p1)
summary(etas.par.p1$model.bg$mod_global)
etas.par.p1$params

# 2)

# ETAS CLASSICO
etas.class.p2 <- etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9,
                           mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1,q = 1.5,
                           params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), w = replicate(nrow(cat.orig), 1),
                           hvarx = replicate(nrow(cat.orig), 1), hvary = replicate(nrow(cat.orig),1),
                           formula1 = "time ~ magnitude + nstaloc_rev + min_distance_rev + distmin - 1", declustering = TRUE,
                           thinning = FALSE, flp = FALSE, ndeclust = 15, onlytime = FALSE,
                           is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
                           compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

summary(etas.class.p2)
plot(etas.class.p2)

# ETAS PARAMETRICO
etas.par.p2 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
                             formula1 = "time ~ magnitude + nstaloc_rev + min_distance_rev + distmin - 1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, bs = "tp", k = 50), process.type.bg = "s2d", spatial.cov.bg = FALSE, mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p2)
plot(etas.par.p2)



# 3)

# ETAS PARAMETRICO
etas.par.p3 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                                 tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                                 mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                                 params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
                                 offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                                 hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                                 ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                                 sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                                 formula.bg = ~ s(x, y, bs = "tp", k = 50) + nstaloc_rev + min_distance_rev + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                                 type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                                 verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

etas.par.p3 <- esecov2_bkgdcov
summary(etas.par.p3)
summary(etas.par.p3$model.bg$mod_global)
plot(etas.par.p3)




# 4)

# ETAS PARAMETRICO
etas.par.p4 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
                             formula1 = "time ~ magnitude + nstaloc_rev + min_distance_rev + distmin - 1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, bs = "tp", k = 50) + nstaloc_rev + min_distance_rev + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p4)
summary(etas.par.p4$model.bg$mod_global)
plot(etas.par.p4)
plot(etas.par.p4$model.bg$mod_global)

save(etas.class.p1, etas.par.p1, etas.class.p2, etas.par.p2, etas.par.p3, etas.par.p4, file = "results_confronto_etas_class_par.RData")




# 5)

# ETAS PARAMETRICO - BEST MODEL
etas.par.p5 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, k = 20) + nstaloc_rev + min_distance_rev + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p5)
summary(etas.par.p5$model.bg$mod_global)
plot(etas.par.p5)


plot_prova_p1_class <- plot(etas.class.p1)
plot_prova_p5 <- plot.etasclass.par(etas.par.p5)



extract_residuals_etas_plot <- function(obj, model_name) {

  # Residui standardizzati modello totale
  res_total <- as.vector((obj$emp1 - obj$teo1) / sqrt(obj$teo1))

  # Residui standardizzati background
  res_back <- as.vector((obj$emp2 - obj$teo2) / sqrt(obj$teo2))

  out <- data.frame(
    model = model_name,
    component = rep(c("Total model", "Background"), each = length(res_total)),
    residual = c(res_total, res_back)
  )

  out <- out %>%
    filter(is.finite(residual))

  return(out)
}

res_df <- bind_rows(
  extract_residuals_etas_plot(plot_prova_p9, "Classical ETAS"),
  extract_residuals_etas_plot(plot_prova_p1_class.trig, "Parametric ETAS")
)

res_df$component <- factor(
  res_df$component,
  levels = c("Total model", "Background")
)

res_df$model <- factor(
  res_df$model,
  levels = c("Classical ETAS", "Parametric ETAS")
)


p_res_box_facet <- ggplot(res_df, aes(x = model, y = residual, fill = model)) +
  geom_boxplot(
    width = 0.6,
    outlier.alpha = 0.65
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ component, nrow = 1) +
  labs(
    x = NULL,
    y = "Standardized residuals",
    fill = NULL,
    title = "Standardized spatial residuals",
    subtitle = "Total model and background component"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 15),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

p_res_box_facet



################################################################################

# confronto modelli ETAS: classico vs. parametrico


#############################
### MODELLO ETAS CLASSICO ###
#############################

# ETAS CLASSICO
etas.class.p1 <- etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9,
                           mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1,q = 1.5, betacov = 0.7,
                           params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), w = replicate(nrow(catalog.withcov), 1),
                           hvarx = replicate(nrow(catalog.withcov), 1), hvary = replicate(nrow(catalog.withcov),1),
                           formula1 = "time ~  magnitude- 1", declustering = TRUE,
                           thinning = FALSE, flp = FALSE, ndeclust = 15, onlytime = FALSE,
                           is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
                           compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

summary(etas.class.p1)
plot_prova_p1_class <- plot(etas.class.p1)



#####################################
### MODELLO ETAS CLASSICO COVTRIG ###
#####################################

# ETAS CLASSICO CON COVARIATA NEL TRIGGERING
etas.class.trig.p1 <- etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9,
                           mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1,q = 1.5, betacov = 0.7,
                           params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), w = replicate(nrow(catalog.withcov), 1),
                           hvarx = replicate(nrow(catalog.withcov), 1), hvary = replicate(nrow(catalog.withcov),1),
                           formula1 = "time ~  magnitude+distmin- 1", declustering = TRUE,
                           thinning = FALSE, flp = FALSE, ndeclust = 15, onlytime = FALSE,
                           is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
                           compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

summary(etas.class.trig.p1)
plot_prova_p1_class.trig <- plot(etas.class.trig.p1)


etas.class.trig.p2 <- etasclass(cat.orig = catalog.withcov, magn.threshold = 2.5, magn.threshold.back = 3.9,
                                mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1,q = 1.5, betacov = 0.7,
                                params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), w = replicate(nrow(catalog.withcov), 1),
                                hvarx = replicate(nrow(catalog.withcov), 1), hvary = replicate(nrow(catalog.withcov),1),
                                formula1 = "time ~  magnitude+distmin- 1", declustering = TRUE,
                                thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
                                is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
                                compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)

summary(etas.class.trig.p2)
plot_prova_p2_class.trig <- plot(etas.class.trig.p2)




################################
### MODELLO ETAS PARAMETRICO ### (2)
################################

# ETAS PARAMETRICO

# 1) s(x, y, k = 40) + distmin # covariata non significativa
# 2) s(x, y, k = 30) + nstaloc_rev + min_distance_rev + distmin, migliore al momento (solo distmin significativa tra le covariate)
# 3) x*y + nstaloc_rev + min_distance_rev + distmin, non va bene, necessario usare le spline
# 4) s(x, y, k = 30) + distmin con distmin nel triggering (modello 8) fino ad ora il migliore

etas.par.p6 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, k = 30) + nstaloc_rev + min_distance_rev + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p6)
summary(etas.par.p6$model.bg$mod_global)
plot_prova_p6 <- plot.etasclass.par(etas.par.p6)


etas.par.p7 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ x*y + nstaloc_rev + min_distance_rev + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p7)
summary(etas.par.p7$model.bg$mod_global)
plot_prova_p7 <- plot.etasclass.par(etas.par.p7)



etas.par.p8 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude+distmin-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, k = 30) + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p8)
summary(etas.par.p8$model.bg$mod_global)
plot_prova_p8 <- plot.etasclass.par(etas.par.p8)



etas.par.p9 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude+distmin-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x, y, k = 25) + distmin, process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p9)
summary(etas.par.p9$model.bg$mod_global)
plot_prova_p9 <- plot.etasclass.par(etas.par.p9)


# MODELLO PARAMETRICO SENZA SPLINE
# 1) x*y + distmin + nstaloc_rev + min_distance_rev + x:distmin + y:distmin

etas.par.p10 <- etasclass.par(cat.orig = catalog.withcov, time.update = FALSE, magn.threshold = 2.5, magn.threshold.back = 3.9,
                             tmax = max(catalog.withcov$time), long.range = range(catalog.withcov$long), lat.range = range(catalog.withcov$lat),
                             mu = 0.3, k0 = 0.02, c = 0.015, p = 1.1, gamma = 0, d = 1, q = 1.5, betacov = 0.7,
                             params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), formula1 = "time~magnitude+distmin-1",
                             offset = 0, hdef = c(1, 1), w = replicate(nrow(catalog.withcov), 1), hvarx = replicate(nrow(catalog.withcov), 1),
                             hvary = replicate(nrow(catalog.withcov), 1), declustering = TRUE, thinning = FALSE, flp = FALSE, m1 = NULL,
                             ndeclust = 15, n.iterweight = 1, onlytime = FALSE, is.backconstant = FALSE, description = "", cat.back = NULL, back.smooth = 1,
                             sectoday = FALSE, longlat.to.km = TRUE, usenlm = TRUE, method = "BFGS", compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36,
                             formula.bg = ~ s(x,y, k = 25) + s(distmin, k = 10) + s(min_distance_rev, k = 10), process.type.bg = "s2d", spatial.cov.bg = TRUE,
                             type.cov.values.bg = list(nstaloc_rev = "interp", min_distance_rev = "interp", distmin = "interp"), mult.bg = 4, ncube.bg = NULL,
                             verbose.bg = TRUE, offset_k.bg = FALSE, grid.bg = FALSE, marked.bg = FALSE, mark.c.bg = FALSE, seed.bg = 2)

summary(etas.par.p10)
summary(etas.par.p10$model.bg$mod_global)
plot_prova_p10 <- plot.etasclass.par(etas.par.p10)


save(etas.class.p1, plot_prova_p1_class, etas.class.trig.p2, plot_prova_p2_class.trig, etas.par.p9, plot_prova_p9, file = "applicazione_ETAS_par.RData")



# IDEA ATTUALE #

# confrontare etas classico vs. etas con covariata nel triggering vs. etas parametrico con covariata nel triggering

# MODELLO CLASSICO
summary(etas.class.p1)
plot_prova_p1_class <- plot(etas.class.p1)

# MODELLO CLASSICO COV TRIG
summary(etas.class.trig.p1)
plot_prova_p1_class.trig <- plot(etas.class.trig.p1)

# MODELLO PARAMETRICO
summary(etas.par.p9)
summary(etas.par.p9$model.bg$mod_global)
plot_prova_p9 <- plot.etasclass.par(etas.par.p9)

extract_residuals_etas_plot <- function(obj, model_name) {

  res_total <- as.vector(
    (obj$emp1 - obj$teo1) / sqrt(pmax(obj$teo1, .Machine$double.eps))
  )

  res_back <- as.vector(
    (obj$emp2 - obj$teo2) / sqrt(pmax(obj$teo2, .Machine$double.eps))
  )

  out <- bind_rows(
    data.frame(
      model = model_name,
      component = "Total model",
      residual = res_total
    ),
    data.frame(
      model = model_name,
      component = "Background",
      residual = res_back
    )
  ) %>%
    filter(is.finite(residual))

  return(out)
}


res_summary <- res_df %>%
  group_by(model, component) %>%
  summarise(
    n = n(),
    min = min(residual, na.rm = TRUE),
    q1 = quantile(residual, 0.25, na.rm = TRUE),
    median = median(residual, na.rm = TRUE),
    mean = mean(residual, na.rm = TRUE),
    q3 = quantile(residual, 0.75, na.rm = TRUE),
    max = max(residual, na.rm = TRUE),
    sd = sd(residual, na.rm = TRUE),
    iqr = IQR(residual, na.rm = TRUE),
    .groups = "drop"
  )

res_summary

res_df <- bind_rows(
  extract_residuals_etas_plot(plot_prova_p1_class,      "ETAS-C"),
  extract_residuals_etas_plot(plot_prova_p2_class.trig, "ETAS-C-FLP"),
  extract_residuals_etas_plot(plot_prova_p9,            "ETAS-P")
)

res_df$component <- factor(
  res_df$component,
  levels = c("Total model", "Background")
)

res_df$model <- factor(
  res_df$model,
  levels = c("ETAS-C", "ETAS-C-FLP", "ETAS-P")
)

p_res_box_facet <- ggplot(res_df, aes(x = model, y = residual, fill = model)) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.55,
    outlier.size = 1.3
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(
    yintercept = c(-1.96, 1.96),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_wrap(~ component, nrow = 1) +
  labs(
    x = NULL,
    y = "Standardized residuals",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 28, color = "black"),

    axis.text.x = element_text(face = "bold", size = 28, color = "black"),
    axis.text.y = element_text(size = 28, color = "black"),
    axis.title.y = element_text(size = 28, color = "black"),

    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.18, "cm"),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),

    legend.position = "none"
  )

p_res_box_facet




italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")

extract_background_grid <- function(obj, model_name,
                                    long_range = range(catalog.withcov$long),
                                    lat_range  = range(catalog.withcov$lat)) {

  ngrid <- round(sqrt(length(obj$back.grid)))

  if (ngrid^2 != length(obj$back.grid)) {
    stop("back.grid does not have a square-grid length.")
  }

  # Usa x.grid e y.grid solo se hanno la lunghezza giusta
  if (!is.null(obj$x.grid) && length(obj$x.grid) == ngrid) {
    xg <- obj$x.grid
  } else {
    xg <- seq(long_range[1], long_range[2], length.out = ngrid)
  }

  if (!is.null(obj$y.grid) && length(obj$y.grid) == ngrid) {
    yg <- obj$y.grid
  } else {
    yg <- seq(lat_range[1], lat_range[2], length.out = ngrid)
  }

  grid_df <- expand.grid(
    long = xg,
    lat  = yg
  )

  grid_df$background <- as.vector(obj$back.grid)
  grid_df$model <- model_name

  grid_df
}


back_df <- bind_rows(
  extract_background_grid(plot_prova_p1_class,      "ETAS-C"),
  extract_background_grid(plot_prova_p1_class.trig, "ETAS-C-cov"),
  extract_background_grid(plot_prova_p9,            "ETAS-P")
)

back_df$model <- factor(
  back_df$model,
  levels = c("ETAS-C", "ETAS-C-cov", "ETAS-P")
)


p_background <- ggplot() +
  geom_tile(
    data = back_df,
    aes(x = long, y = lat, fill = background)
  ) +
  geom_sf(
    data = italy,
    fill = NA,
    color = "black",
    linewidth = 0.35
  ) +
  coord_sf(
    xlim = c(5.5, 19.8),
    ylim = c(36.0, 47.8),
    expand = FALSE
  ) +
  facet_wrap(~ model, nrow = 1) +
  scale_fill_viridis_c(
    name = expression(hat(lambda)[bkgd]),
    option = "plasma",
    trans = "sqrt",
    guide = guide_colorbar(
      barheight = unit(16, "cm"),
      barwidth  = unit(1, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 20, color = "black"),
    axis.text = element_text(size = 24, color = "black"),
    axis.title = element_text(size = 28, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_text(size = 28, color = "black"),
    legend.text = element_text(size = 24, color = "black"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  )

p_background



extract_triggered_grid <- function(obj, model_name,
                                   long_range = range(catalog.withcov$long),
                                   lat_range  = range(catalog.withcov$lat)) {

  ngrid <- round(sqrt(length(obj$trig.grid)))

  if (ngrid^2 != length(obj$trig.grid)) {
    stop("trig.grid does not have a square-grid length.")
  }

  if (!is.null(obj$x.grid) && length(obj$x.grid) == ngrid) {
    xg <- obj$x.grid
  } else {
    xg <- seq(long_range[1], long_range[2], length.out = ngrid)
  }

  if (!is.null(obj$y.grid) && length(obj$y.grid) == ngrid) {
    yg <- obj$y.grid
  } else {
    yg <- seq(lat_range[1], lat_range[2], length.out = ngrid)
  }

  grid_df <- expand.grid(
    long = xg,
    lat  = yg
  )

  grid_df$triggered <- as.vector(obj$trig.grid)
  grid_df$model <- model_name

  grid_df
}


trig_df <- bind_rows(
  extract_triggered_grid(plot_prova_p1_class,      "ETAS-C"),
  extract_triggered_grid(plot_prova_p1_class.trig, "ETAS-C-cov"),
  extract_triggered_grid(plot_prova_p9,            "ETAS-P")
)

trig_df$model <- factor(
  trig_df$model,
  levels = c("ETAS-C", "ETAS-C-cov", "ETAS-P")
)


p_triggered <- ggplot() +
  geom_tile(
    data = trig_df,
    aes(x = long, y = lat, fill = triggered)
  ) +
  geom_sf(
    data = italy,
    fill = NA,
    color = "black",
    linewidth = 0.35
  ) +
  coord_sf(
    xlim = c(5.5, 19.8),
    ylim = c(36.0, 47.8),
    expand = FALSE
  ) +
  facet_wrap(~ model, nrow = 1) +
  scale_fill_viridis_c(
    name = expression(hat(lambda)[trig]),
    option = "plasma",
    trans = "sqrt",
    guide = guide_colorbar(
      barheight = unit(16, "cm"),
      barwidth  = unit(1, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 20, color = "black"),
    axis.text = element_text(size = 24, color = "black"),
    axis.title = element_text(size = 28, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_text(size = 28, color = "black"),
    legend.text = element_text(size = 24, color = "black"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  )

p_triggered
