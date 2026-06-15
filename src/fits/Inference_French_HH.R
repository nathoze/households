# ======================================================================
# Run inference and save outputs in the results/ directory 
# ======================================================================


rm(list = ls())

suppressPackageStartupMessages({
  library(here)
})

source(here::here("src", "main.R"))
source(here::here("src", "paths.R"))

set.seed(1954)

# ======================================================================
# Global settings
# ======================================================================

fit_epid_model <- "exp"

gi <- 5
ds <- 5
alpha <- 1
vp <- 6
n_hh <- 100

sar_values <- c(0.15, 0.3, 0.45)
index_values <- 1:100

iter_warmup <- 1500

overwrite_chains <- FALSE
overwrite_ct <- FALSE
reconstruct_ct <- TRUE

# For testing, uncomment:
# sar_values <- c(0.15)
# index_values <- 1:2

# ======================================================================
# Path helpers
# ======================================================================

make_household_filename <- function(gi, ds, sar, alpha, n_hh, index) {
  sprintf(
    "Households_26082025_gi_%s_incubation_%s_sar_%s_alpha_%s_NHH_%s_index_%s.rds",
    gi,
    ds,
    as.character(sar),
    alpha,
    n_hh,
    index
  )
}

data_file_path <- function(filename) {
  file.path(here::here("data", "simulated_data"), filename)
}

results_file_path <- function(kind = c("chains", "ct"), model, filename) {
  kind <- match.arg(kind)

  out_dir <- here::here("results", kind, model)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  file.path(out_dir, filename)
}

check_input_data <- function(path) {
  if (!file.exists(path)) {
    stop("Missing data file:\n", path, call. = FALSE)
  }

  invisible(path)
}

# ======================================================================
# Model-specific data modification
# ======================================================================

remove_contact_symptoms <- function(xdata) {
  contact_ids <- which(xdata$col_infected > 1)
  xdata$symptomatic[contact_ids] <- 0
  xdata
}

# ======================================================================
# Model-specific inference and reconstruction
# ======================================================================

run_model_inference <- function(model, xdata) {
  if (model == "quantitative_vl_with_contact_symptoms") {
    return(
      run_inference(
        xdata,
        init,
        rstan_model,
        iter_warmup = iter_warmup
      )
    )
  }

  if (model == "quantitative_vl_no_contact_symptoms") {
    return(
      run_inference(
        xdata,
        init,
        rstan_model = rstan_model_exp,
        iter_warmup = iter_warmup
      )
    )
  }

  if (model == "qualitative_vl_no_contact_symptoms") {
    return(
      run_inference(
        xdata,
        init,
        rstan_model = rstan_model_positivity_3,
        iter_warmup = iter_warmup
      )
    )
  }

  if (model == "symptoms_only") {
    return(
      run_inference_symptoms(
        xdata,
        init,
        rstan_model = rstan_model_symptoms_gamma,
        iter_warmup = iter_warmup
      )
    )
  }

  if (model == "constant") {
    return(
      run_inference(
        xdata,
        init,
        rstan_model_constant,
        iter_warmup = iter_warmup
      )
    )
  }

  stop("Unknown model: ", model, call. = FALSE)
}

reconstruct_one_household <- function(model, xdata, Chains, hh_index) {
  if (model == "quantitative_vl_with_contact_symptoms") {
    return(
      chain_transmission_reconstruction(
        xdata,
        Chains = Chains,
        hh_index,
        model = fit_epid_model
      )
    )
  }

  if (model == "quantitative_vl_no_contact_symptoms") {
    return(
      reconstruct_transmission_chain(
        xdata,
        Chains = Chains,
        hh_index,
        model = fit_epid_model
      )
    )
  }

  if (model == "qualitative_vl_no_contact_symptoms") {
    return(
      reconstruct_transmission_chain(
        xdata,
        Chains = Chains,
        hh_index,
        model = fit_epid_model
      )
    )
  }

  if (model == "symptoms_only") {
    return(
      reconstruct_transmission_chain_symptoms(
        xdata,
        Chains = Chains,
        hh_index
      )
    )
  }

  if (model == "constant") {
    return(
      chain_transmission_reconstruction(
        xdata,
        Chains = Chains,
        hh_index,
        model = "constant"
      )
    )
  }

  stop("Unknown model: ", model, call. = FALSE)
}

# ======================================================================
# Core runner
# ======================================================================

run_one_fit <- function(model, gi, ds, sar, alpha, n_hh, index) {
  filename <- make_household_filename(
    gi = gi,
    ds = ds,
    sar = sar,
    alpha = alpha,
    n_hh = n_hh,
    index = index
  )

  data_path <- data_file_path(filename)
  chains_path <- results_file_path("chains", model, filename)
  ct_path <- results_file_path("ct", model, filename)

  check_input_data(data_path)

  load(data_path) # expected: xdata, parameters, init, ...

  if (model %in% c(
    "quantitative_vl_no_contact_symptoms",
    "qualitative_vl_no_contact_symptoms"
  )) {
    xdata <- remove_contact_symptoms(xdata)
  }

  message("------------------------------------------------------------")
  message("Model: ", model)
  message("SAR: ", sar, " | index: ", index)
  message("Chains: ", chains_path)
  message("CT: ", ct_path)

  if (!file.exists(chains_path) || overwrite_chains) {
    message("Running inference...")

    res <- run_model_inference(model = model, xdata = xdata)

    Chains <- res$Chains
    DIC <- res$DIC
    pD <- res$pD

    message("DIC: ", DIC)

    save(
      Chains,
      DIC,
      pD,
      file = chains_path
    )
  } else {
    message("Chains already exist; loading existing file.")
    load(chains_path) # expected: Chains, DIC, pD
  }

  if (reconstruct_ct && (!file.exists(ct_path) || overwrite_ct)) {
    message("Reconstructing transmission chains...")

    CT <- vector("list", xdata$n_households)

    for (hh_index in seq_len(xdata$n_households)) {
      message("  household ", hh_index, "/", xdata$n_households)

      CT[[hh_index]] <- reconstruct_one_household(
        model = model,
        xdata = xdata,
        Chains = Chains,
        hh_index = hh_index
      )

      # Save progressively so long jobs are not lost if interrupted.
      save(CT, file = ct_path)
    }
  } else if (file.exists(ct_path)) {
    message("CT already exists; skipping reconstruction.")
  }

  invisible(
    data.frame(
      model = model,
      gi = gi,
      ds = ds,
      sar = sar,
      alpha = alpha,
      n_hh = n_hh,
      index = index,
      data_path = data_path,
      chains_path = chains_path,
      ct_path = ct_path,
      chains_exists = file.exists(chains_path),
      ct_exists = file.exists(ct_path)
    )
  )
}

run_model_grid <- function(model,
                           sar_values = c(0.15, 0.3, 0.45),
                           index_values = 1:100,
                           gi = 5,
                           ds = 5,
                           alpha = 1,
                           n_hh = 100) {
  rows <- list()
  k <- 0L

  for (sar in sar_values) {
    for (index in index_values) {
      k <- k + 1L

      rows[[k]] <- run_one_fit(
        model = model,
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n_hh,
        index = index
      )
    }
  }

  do.call(rbind, rows)
}

# ======================================================================
# Models to run
# ======================================================================

models_to_run <- c(
  "quantitative_vl_with_contact_symptoms",
  "symptoms_only",
  "constant",
  "qualitative_vl_no_contact_symptoms",
  "quantitative_vl_no_contact_symptoms"
)

# Keep constant at SAR 0.45 by default, matching the old script.
# Change to sar_values if you want all three SAR values.
sar_values_by_model <- list(
  quantitative_vl_with_contact_symptoms = sar_values,
  symptoms_only = sar_values,
  constant = c(0.45),
  qualitative_vl_no_contact_symptoms = sar_values,
  quantitative_vl_no_contact_symptoms = sar_values
)

# ======================================================================
# Run
# ======================================================================

manifest_rows <- list()
m <- 0L

for (model in models_to_run) {
  m <- m + 1L

  manifest_rows[[m]] <- run_model_grid(
    model = model,
    sar_values = sar_values_by_model[[model]],
    index_values = index_values,
    gi = gi,
    ds = ds,
    alpha = alpha,
    n_hh = n_hh
  )
}
message("Inference script completed.")

# manifest <- do.call(rbind, manifest_rows)

# dir.create(here::here("results"), recursive = TRUE, showWarnings = FALSE)
# 
# write.csv(
#   manifest,
#   file = here::here("results", "manifest_inference_outputs.csv"),
#   row.names = FALSE
# )

 # message("Manifest written to: ", here::here("results", "manifest_inference_outputs.csv"))
