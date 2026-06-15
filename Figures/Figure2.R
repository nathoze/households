 
# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(scales)
  library(ggpubr)
  library(purrr)
  library(here)
})

source(here::here("src", "paths.R"))

cfg <- list(
  root = here::here("src"),
  fit_model = "quantitative_vl_no_contact_symptoms",
  h = 0.5,
  gi = 5,
  ds = 5,
  sar = 0.30,
  alpha = 1,
  n_samples_prior = 1e6,
  out_prefix = "Figure2",
  out_dir = here::here("Figures"),
  res_dpi = 300
)

source(here::here("src", "main.R"))
source(here::here("src", "visualization", "MCMC_predictions.R"))
source(here::here("src", "visualization", "MCMC_compare_input_estimates.R"))
source(here::here("src", "models", "transmission probability conditional.R"))
source(here::here("src", "helpers.R"))

# Backward-compatible color names used below
if (!exists("colors")) {
  colors <- list()
}
if (is.null(colors$data_constant)) {
  colors$data_constant <- if (!is.null(colors$input)) colors$input else "#5C9B58"
}
if (is.null(colors$data_exp)) {
  colors$data_exp <- if (!is.null(colors$fit)) colors$fit else "#FF8A5B"
}


# ======================================================================
# Results-system helpers
# ======================================================================
# Expected results structure:
#
# results/
#   chains/
#     symptoms_only/
#     constant/
#     qualitative_vl_no_contact_symptoms/
#     quantitative_vl_no_contact_symptoms/
#     quantitative_vl_contact_symptoms/
#   ct/
#     symptoms_only/
#     constant/
#     qualitative_vl_no_contact_symptoms/
#     quantitative_vl_no_contact_symptoms/
#     quantitative_vl_contact_symptoms/
#
# Figure 2 main panel uses:
#   quantitative_vl_no_contact_symptoms
#
# Figure S2 compares:
#   symptoms_only
#   qualitative_vl_no_contact_symptoms
#   quantitative_vl_no_contact_symptoms
#
# If your paths.R uses the alternative name
# "quantitative_vl_with_contact_symptoms", the helper below also supports it.

resolve_model_name <- function(model) {
  aliases <- c(
    quantitative_vl_with_contact_symptoms = "quantitative_vl_contact_symptoms"
  )

  if (model %in% names(aliases)) {
    return(unname(aliases[[model]]))
  }

  model
}

load_fit_objects <- function(model, gi, ds, sar, alpha, n_hh, index) {
  model <- resolve_model_name(model)

  paths <- load_household_fit(
    model = model,
    gi = gi,
    ds = ds,
    sar = sar,
    alpha = alpha,
    n_hh = n_hh,
    index = index
  )

  missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing_paths) > 0) {
    stop(
      "Missing file(s) for model '", model, "':
",
      paste(missing_paths, collapse = "
"),
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())

  load(paths$data_path,   envir = env)  # expected: xdata, parameters, ...
  load(paths$chains_path, envir = env)  # expected: Chains
  load(paths$ct_path,     envir = env)  # expected: CT

  list(
    xdata = env$xdata,
    parameters = env$parameters,
    Chains = env$Chains,
    CT = env$CT,
    paths = paths
  )
}

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

 
 
# ---- Chargement des données  ----
gi <- cfg$gi
ds <- cfg$ds
sar <- cfg$sar
alpha_in <- cfg$alpha
h <- cfg$h
I <- 12
n_hh <- 100

fit_main <- load_fit_objects(
  model = "quantitative_vl_no_contact_symptoms",
  gi = gi,
  ds = ds,
  sar = sar,
  alpha = alpha_in,
  n_hh = n_hh,
  index = I
)

xdata <- fit_main$xdata
parameters <- fit_main$parameters
Chains <- fit_main$Chains
CT <- fit_main$CT



# ---- Plots A/B (viral traj + graph) ----
hh_id <- 79
ct <- CT[[hh_id]]

 
# Main Figure 2 uses the quantitative viral-load model without symptom information for contacts.
# To match that analysis dataset, symptoms are removed for infected contacts only.
contact.id <- which(xdata$col_infected > 1)
xdata$symptomatic[contact.id] <- 0


# A: viral load trajectories (ta fonction)
g1 <- MCMC_plot_predicted_viral_load_hhindex(Chains = Chains, xdata = xdata, Xlim = 30, hh.index = hh_id)

p1 <- g1 +
  geom_vline(xintercept = c(0, 3, 7, 14, 21, 28), linetype = "dashed") +
  xlab("Days post-inclusion")


# ---- Input distributions (priors simulés) vs Fit (postérieurs MCMC) ----

n.samples <- cfg$n_samples_prior
 m <- xdata$all.HH$m[1]
# Input (priors)
T_P_input <- parameters$T_P_star * exp(rnorm(n.samples, sd = parameters$T_P_eta))
T_U_input <- parameters$T_U_star * exp(rnorm(n.samples, sd = parameters$T_U_eta))
V_P_input <- parameters$V_P_star * exp(rnorm(n.samples, sd = parameters$V_P_eta))

# GI_input <- (1 + 4 / V_P_input) * T_P_input
## -- NEW GI INPUT WEIGHTED BY THE PROBABILITY TO TRANSMIT

 
parameters$m  <-m
parameters$fixed_param = FALSE

sim <- simulate_transmissions_input(
  parameters = parameters,
  N = 20000,
  t_grid = seq(0, 35, by = 0.05)
)


GI_input = sim$summary$t_weighted 
SAR_input = sim$summary$P

input.gi = sim$summary$t_conditional


incubation_input <- parameters$Delta_S_star * exp(rnorm(n.samples, sd = parameters$Delta_S_eta))


VL =NULL
for(hh_id in 1:100){
  VL =rbind(VL,get_VL(Chains,xdata = xdata,id=hh_id))
}

sar.fit = pinf_from_vl(T_P =  VL$T_P,T_U = VL$T_U,V_P = VL$V_P,m = VL$m,model = "exp")
print(c(mean(sar.fit), median(sar.fit), quantile975(sar.fit)))



make_network_plot <- function(ct, node_radius = 0.08) {
  EDGES <- matrix(0, nrow = 8, ncol = 8)
  EDGES[1:4, 1:4] <- ct$input.chain
  EDGES[5:8, 5:8] <- ct$P
  
  nodes <- data.frame(
    id = 1:8,
    x = c(1, 2, 1.2, 3, 3.5, 4.5, 3.7, 5.5),
    y = c(1.6, 1.6, 1.1, 1.1, 1.6, 1.6, 1.1, 1.1),
    fill = rep(c("#c45441", "#1212ca", "#128a18", "#ffb20d"), 2),
    label = rep(1:4, 2)
  )
  
  edges <- expand.grid(from = nodes$id, to = nodes$id) |>
    subset(((from > 4 & to > 4) | (from <= 4 & to <= 4)) & from != to) |>
    as_tibble()
  
  # poids : on garde ton choix EDGES[to, from]
  edges <- edges |>
    mutate(
      weight = purrr::pmap_dbl(list(to, from), ~ EDGES[..1, ..2])
    ) |>
    filter(weight > 0.01) |>
    mutate(
      size = scales::rescale(weight, to = c(0.5, 3)),
      x_start = nodes$x[from],
      y_start = nodes$y[from],
      x_end0  = nodes$x[to],
      y_end0  = nodes$y[to],
      x_end   = x_end0 - (x_end0 - x_start) * node_radius,
      y_end   = y_end0 - (y_end0 - y_start) * node_radius,
      x_mid   = (x_start + x_end) / 2,
      y_mid   = (y_start + y_end) / 2,
      label   = paste0(round(weight * 100), "%")
    )
  
  ggplot() +
    geom_segment(
      data = edges,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end, size = size),
      arrow = arrow(length = unit(0.4, "inches"), type = "closed"),
      color = "grey"
    ) +
    geom_text(data = edges, aes(x = x_mid, y = y_mid, label = label),
              size = 8, fontface = "bold", color = "black") +
    geom_point(data = nodes, aes(x = x, y = y, fill = fill),
               shape = 21, size = 16, color = "black") +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = 14, fontface = "bold", color = "lightgrey") +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

p2 <- make_network_plot(ct)
 
## fit.gi 
fit.incubation = c()
for(i in 1:xdata$n_infected_total){
  fit.incubation=c(fit.incubation, as.numeric(Chains$TI_tmp_positive[,i]) )
}
input.gi = fit.gi = c()
for(i in 1:xdata$n_households){
  ct = CT[[i]]
  # input.gi= c(input.gi, ct$input.generation.times)
  fit.gi= c(fit.gi, ct$generation.times)
}

   
data.comparison <- bind_rows(
  tibble(x = fit.gi,           type = "fit",   figure = "GI"),
  tibble(x = GI_input,         type = "input", figure = "GI"),
  tibble(x = fit.incubation,   type = "fit",   figure = "Incubation"),
  tibble(x = incubation_input, type = "input", figure = "Incubation"),
  tibble(x = sar.fit,          type = "fit",   figure = "SAR"),
  tibble(x = SAR_input,        type = "input", figure = "SAR")
)

# ---- Plot helper: histogram overlay "input vs fit" ----
# ---- Density overlay "input vs fit" ----
plot_hist_compare <- function(df, fig, xlab,
                              fill_colors = c("input" = colors$data_constant,
                                              "fit"   = colors$data_exp),
                              alpha = 0.5) {
  
  df %>%
    filter(tolower(figure) == tolower(fig)) %>%
    mutate(type = as.factor(type)) %>%
    ggplot(aes(x = x, fill = type)) +
    geom_density(alpha = alpha, show.legend = TRUE) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 16, angle = 0, vjust = 0),
      axis.text.y = element_text(size = 16),
      text = element_text(size = 16),
      strip.background = element_rect(fill = TRUE)
    ) +
    xlab(xlab) +
    ylab("Proportion") +
    scale_fill_manual(values = fill_colors)
}

histogram.GI <- plot_hist_compare(data.comparison, "GI", "Generation interval (days)")
histogram.Incubation <- plot_hist_compare(data.comparison, "Incubation", "Incubation period (days)")
histogram.SAR <- plot_hist_compare(data.comparison, "SAR", "Transmission probability")

# ---- right_column : remplace p3..p6 ----
right_column <- ggarrange(
  histogram.GI,
  histogram.Incubation,
  histogram.SAR,
  ncol = 1, nrow = 3,
  labels = c("C", "D", "E")
)

# ---- Graph "B" : réseau d'infection (factorisé) ----
 
# ---- Layout final ----

final_plot <- ggarrange(
  ggarrange(p1, p2, ncol = 1, labels = c("A", "B")),
  right_column,
  ncol = 2,
  widths = c(1, 1)
)

print(final_plot)

# ---- Export ----
A <- 2

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(cfg$out_dir, "Figure2.pdf")

ggpubr::ggexport(
  final_plot,
  filename = pdf_file,
  width = 8 * A,
  height = 4 * A
)
 
# ======================================================================
# Figure S2 — Inferred transmission chains for alternative study designs
# ======================================================================

plot_transmission_chain_for_model <- function(model,
                                              panel_title = NULL,
                                              hh_id,
                                              gi_value,
                                              ds_value,
                                              sar_value,
                                              alpha_value,
                                              n_hh_value,
                                              index_value) {
  fit <- load_fit_objects(
    model = model,
    gi = gi_value,
    ds = ds_value,
    sar = sar_value,
    alpha = alpha_value,
    n_hh = n_hh_value,
    index = index_value
  )

  CT <- fit$CT

  if (!hh_id %in% seq_along(CT)) {
    stop("Household ", hh_id, " is not available in CT for model: ", model)
  }

  p <- make_network_plot(CT[[hh_id]])

  if (!is.null(panel_title)) {
    p <- p + ggtitle(panel_title)
  }

  p
}

S2_clinical <- plot_transmission_chain_for_model(
  model = "symptoms_only",
  panel_title = "Symptoms only",
  hh_id = 79,
  gi_value = gi,
  ds_value = ds,
  sar_value = sar,
  alpha_value = alpha_in,
  n_hh_value = n_hh,
  index_value = I
)

S2_qualitative <- plot_transmission_chain_for_model(
  model = "qualitative_vl_no_contact_symptoms",
  panel_title = "Qualitative VL, no contact symptoms",
  hh_id = 79,
  gi_value = gi,
  ds_value = ds,
  sar_value = sar,
  alpha_value = alpha_in,
  n_hh_value = n_hh,
  index_value = I
)

S2_quantitative <- plot_transmission_chain_for_model(
  model = "quantitative_vl_no_contact_symptoms",
  panel_title = "Quantitative VL, no contact symptoms",
  hh_id = 79,
  gi_value = gi,
  ds_value = ds,
  sar_value = sar,
  alpha_value = alpha_in,
  n_hh_value = n_hh,
  index_value = I
)

Figure_S2 <- ggpubr::ggarrange(
  S2_clinical,
  S2_qualitative,
  S2_quantitative,
  ncol = 3,
  labels = c("A", "B", "C")
)

print(Figure_S2)

# ---- Export Figure S2 ----
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

ggpubr::ggexport(
  Figure_S2,
  filename = file.path(cfg$out_dir, "Figure_S2_transmission_chains.pdf"),
  width = 15,
  height = 5
)
 