 
# ======================================================================
# Paths for household simulation project
# ======================================================================

household_filename <- function(gi, ds, sar, alpha, n_hh, index) {
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

load_household_fit <- function(model,
                               gi, ds, sar, alpha, n_hh, index,
                               root = here::here()) {
  
  model <- match.arg(
    model,
    allowed_models <- c(
      "symptoms_only",
      "constant",
      "qualitative_vl_no_contact_symptoms",
      "quantitative_vl_no_contact_symptoms",
      "quantitative_vl_with_contact_symptoms"
    )  )
  
  filename <- household_filename(
    gi = gi,
    ds = ds,
    sar = sar,
    alpha = alpha,
    n_hh = n_hh,
    index = index
  )
  
  list(
    data_path = file.path(
      root,
      "data",
      "simulated_data",
      filename
    ),
    chains_path = file.path(
      root,
      "results",
      "chains",
      model,
      filename
    ),
    ct_path = file.path(
      root,
      "results",
      "ct",
      model,
      filename
    )
  )
}

load_fit_objects <- function(model,
                             gi, ds, sar, alpha, n_hh, index,
                             root = here::here()) {
  paths <- load_household_fit(
    model = model,
    gi = gi,
    ds = ds,
    sar = sar,
    alpha = alpha,
    n_hh = n_hh,
    index = index,
    root = root
  )
  
  missing <- c()
  
  if (!file.exists(paths$data_path)) {
    missing <- c(missing, paths$data_path)
  }
  
  if (!file.exists(paths$chains_path)) {
    missing <- c(missing, paths$chains_path)
  }
  
  if (!file.exists(paths$ct_path)) {
    missing <- c(missing, paths$ct_path)
  }
  
  if (length(missing) > 0) {
    stop(
      "Missing file(s) for model '", model, "':\n",
      paste(missing, collapse = "\n")
    )
  }
  
  env <- new.env(parent = emptyenv())
  
  load(paths$data_path, envir = env)
  load(paths$chains_path, envir = env)
  load(paths$ct_path, envir = env)
  
  list(
    xdata = env$xdata,
    parameters = env$parameters,
    Chains = env$Chains,
    CT = env$CT,
    paths = paths
  )
}