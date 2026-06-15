# ======================================================================
# Figure S1 — Parameter recovery for one simulation
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(here)
  library(HDInterval)
})

source(here::here("src", "paths.R"))
source(here::here("src", "main.R"))
source(here::here("src", "visualization", "MCMC_predictions.R"))
source(here::here("src", "visualization", "MCMC_compare_input_estimates.R"))
source(here::here("src", "helpers.R"))

cfg <- list(
  model = "quantitative_vl_with_contact_symptoms",
  gi = 5,
  ds = 5,
  sar = 0.30,
  alpha = 1,
  n_hh = 100,
  index = 12,
  out_dir = here::here("Figures")
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

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
      "Missing file(s) for model '", model, "':\n",
      paste(missing_paths, collapse = "\n"),
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  load(paths$data_path,   envir = env)
  load(paths$chains_path, envir = env)
  load(paths$ct_path,     envir = env)

  list(
    xdata = env$xdata,
    parameters = env$parameters,
    Chains = env$Chains,
    CT = env$CT,
    paths = paths
  )
}

fit <- load_fit_objects(
  model = cfg$model,
  gi = cfg$gi,
  ds = cfg$ds,
  sar = cfg$sar,
  alpha = cfg$alpha,
  n_hh = cfg$n_hh,
  index = cfg$index
)

xdata <- fit$xdata
parameters <- fit$parameters
Chains <- fit$Chains
CT <- fit$CT

plot.parameters <- MCMC_compare_input_estimates_residuals(
  Chains = Chains,
  xdata = xdata,
  lims = c(4, 4, 4, 4, 6)
)

p1 <- plot.parameters[[1]] + xlab("Infection time - input") + ylab("Bias")
p2 <- plot.parameters[[2]] + xlab("Time to peak - input") + ylab("Bias")
p3 <- plot.parameters[[3]] + xlab("Time to clearance - input") + ylab("Bias")
p4 <- plot.parameters[[4]] + xlab("Peak value - input") + ylab("Bias")

FigureS1 <- ggpubr::ggarrange(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  p4 + theme(legend.position = "none"),
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D")
)

print(FigureS1)

summary.stats <- MCMC_compare_summary_statistics(Chains = Chains, xdata = xdata)
name <- c("tinfection", "T_P", "T_U", "V_P")

for (i in c(2, 3, 4, 1)) {
  message(name[i])
  S <- summary.stats[[i]]
  print(S %>% summarise(coverage = mean(covered, na.rm = TRUE)))
  print(
    S %>%
      summarise(
        mean_bias = mean(bias, na.rm = TRUE),
        lower = q025(bias),
        upper = q975(bias)
      )
  )
}

summary.stats.index.contact <- MCMC_compare_summary_statistics_index_contact(
  Chains = Chains,
  xdata = xdata
)

for (i in c(2, 3, 4, 1)) {
  message(name[i])
  S <- summary.stats.index.contact[[i]]
  S_index <- S %>% filter(status == "index")
  S_contact <- S %>% filter(status == "contact")

  print(tibble::tibble(
    parameter = name[i],
    coverage_index = mean(S_index$covered, na.rm = TRUE),
    coverage_contact = mean(S_contact$covered, na.rm = TRUE)
  ))
}

ggpubr::ggexport(
  FigureS1,
  filename = file.path(cfg$out_dir, "FigureS1.pdf"),
  width = 8,
  height = 6
)
