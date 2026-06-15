 
# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)     # ggarrange
  library(scales)
  library(transport)  # wasserstein1d
})

library(here)

source(here::here("src", "paths.R"))

# ---- Config ----
cfg <- list(
  fit_epid_model = "exp",
  h              = 0.5,
  out_prefix     = "Figure3",
  out_dir        = here::here("Figures"),
  res_dpi        = 300
)

# ---- Source project code ----
source(here::here("src", "models", "functions_household_transmission.R"))
source(here::here("src", "main.R"))
source(here::here("src", "visualization", "MCMC_predictions.R"))
source(here::here("src", "models", "transmission probability conditional.R"))
source(here::here("src", "helpers.R"))

# ---- Clean loader for the reorganized results directory ----
load_fit_objects <- function(model, gi, ds, sar, alpha, n_hh, index) {
  model_candidates <- switch(
    model,
    "quantitative_vl_contact_symptoms" = c(
      "quantitative_vl_contact_symptoms",
      "quantitative_vl_with_contact_symptoms"
    ),
    model
  )

  last_error <- NULL

  for (model_candidate in model_candidates) {
    paths <- tryCatch(
      load_household_fit(
        model = model_candidate,
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n_hh,
        index = index
      ),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (is.null(paths)) next

    missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
    if (length(missing_paths) > 0) {
      last_error <- simpleError(
        paste0(
          "Missing file(s) for model '", model_candidate, "':\n",
          paste(missing_paths, collapse = "\n")
        )
      )
      next
    }

    env <- new.env(parent = emptyenv())

    load(paths$data_path,   envir = env)  # expected: xdata, parameters, ...
    load(paths$chains_path, envir = env)  # expected: Chains
    load(paths$ct_path,     envir = env)  # expected: CT

    return(
      list(
        xdata = env$xdata,
        parameters = env$parameters,
        Chains = env$Chains,
        CT = env$CT,
        paths = paths,
        model_used = model_candidate
      )
    )
  }

  stop(
    "Could not load model '", model, "'. Tried: ",
    paste(model_candidates, collapse = ", "),
    "\nLast error:\n",
    conditionMessage(last_error)
  )
}

 

# Transmission probability for a given set of VL parameters
pinf_from_params <- function(T_P, T_U, V_P, m, model, h = cfg$h) {
  if (model == "exp") {
    1 - exp(-m * 10^(2*h) / (h * V_P * log(10)) * (T_P + T_U) * (10^(h * V_P) - 1))
  } else if (model == "constant") {
    1 - exp(-m * (T_P + T_U))
  } else stop("Unknown model: ", model)
}

sar_household <- function(pinf0, hh.size, alpha) {
  1 - (1 - pinf0 / (hh.size - 1)^alpha)^(hh.size - 1)
}

# Viral load parameter draws from Chains (your “model_1” version)
get_VL_model_1 <- function(Chains, xdata, id) {
  data.frame(
    T_P   = exp(Chains$log_T_P_mu + Chains$log_T_P_sd * Chains$T_P_eta[, id]) * xdata$T_P_star,
    T_U   = exp(Chains$log_T_U_mu + Chains$log_T_U_sd * Chains$T_U_eta[, id]) * xdata$T_U_star,
    V_P   = exp(Chains$log_V_P_mu + Chains$log_V_P_sd * Chains$V_P_eta[, id]) * xdata$V_P_star,
    m     = 10^(Chains$logm),
    alpha = Chains$alpha
  )
}

# Add SAR per row (pinf0 for individual link)
get_SAR_Chains <- function(VL, model) {
  SAR <- pinf_from_params(VL$T_P, VL$T_U, VL$V_P, VL$m, model = model, h = cfg$h)
  VL %>% mutate(SAR = SAR)
}

# Add household SAR using alpha + hh.size
get_SAR_Chains_hh_size <- function(VL, model, hh.size) {
  pinf0 <- pinf_from_params(VL$T_P, VL$T_U, VL$V_P, VL$m, model = model, h = cfg$h)
  SAR   <- sar_household(pinf0, hh.size = hh.size, alpha = VL$alpha)
  VL %>% mutate(SAR = SAR)
}

compare_distributions <- function(x1, x2, comparison.type = c("Mean","KS","WS")) {
  comparison.type <- match.arg(comparison.type)
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  if (comparison.type == "WS") return(transport::wasserstein1d(x1, x2))
  if (comparison.type == "KS") return(as.numeric(stats::ks.test(x1, x2)$statistic))
  mean(x1) - mean(x2)
}
hdi_contains <- function(v, X, mass = 0.95) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) stop("Need >=2 finite values.")
  
  if (!requireNamespace("HDInterval", quietly = TRUE)) {
    stop("Package 'HDInterval' required: install.packages('HDInterval').")
  }
  
  dens <- density(v)
  hdi_obj <- HDInterval::hdi(dens, credMass = mass)
  
  if (is.null(dim(hdi_obj))) {
    # vector
    lower <- hdi_obj[1]
    upper <- hdi_obj[2]
    return(isTRUE(X >= lower & X <= upper))
  } else {
    # matrix (k x 2)
    return(any(X >= hdi_obj[, 1] & X <= hdi_obj[, 2]))
  }
}

# ======================================================================
# Core evaluation function  
# ======================================================================
compare_fit_input_model <- function(xdata, CT_input, ch, parameters, type,
                                    fit.epid.model = cfg$fit_epid_model) {
  
  # ---- Reconstruction scores per household ----
  df1 <- df2 <- NULL
  for (i in seq_len(xdata$n_households)) {
    ct <- CT_input[[i]]
    df1 <- bind_rows(df1, data.frame(
      index = i,
      size_households = xdata$size_households[i],
      n.infected = xdata$n_infected[i],
      good.reconstruction = mean(ct$number.generations == ct$input.number.generations)
    ))
    df2 <- bind_rows(df2, data.frame(
      index = i,
      size_households = xdata$size_households[i],
      n.infected = xdata$n_infected[i],
      good.reconstruction = mean(ct$difference.edges == 0)
    ))
  }
  
  M1 <- mean(df1 %>% filter(n.infected > 2) %>% pull(good.reconstruction), na.rm = TRUE)
  M2 <- mean(df2 %>% filter(n.infected > 2) %>% pull(good.reconstruction), na.rm = TRUE)
  
  # ---- Generation interval: pooled over households ----
  input.gi <- fit.gi <- c()
  for (i in seq_len(xdata$n_households)) {
    ct <- CT_input[[i]]
    input.gi <- c(input.gi, ct$input.generation.times)
    fit.gi   <- c(fit.gi,   ct$generation.times)
  }
  
  # ---- Incubation: pooled over infected ----
  input.incubation <- fit.incubation <- c()
  if (type != "Symptoms") {
    for (i in seq_len(xdata$n_infected_total)) {
      input.incubation <- c(input.incubation, xdata$all.HH$incubation[i])
      fit.incubation   <- c(fit.incubation, as.numeric(ch$TI_tmp_positive[, i]))
    }
  } else {
    for (i in seq_len(xdata$n_infected_total)) {
      input.incubation <- c(input.incubation, xdata$all.HH$incubation[i])
      fit.incubation   <- c(fit.incubation, as.numeric(ch$incubation))
    }
  }
  
  # ---- Alpha bias (alpha.input hard-coded as in your script) ----
  alpha.input <- 1
  bias.alpha.mean  <- mean(ch$alpha - alpha.input, na.rm = TRUE)
  bias.alpha.lower <- q025(ch$alpha - alpha.input)
  bias.alpha.upper <- q975(ch$alpha - alpha.input)
  
  # ---- SAR / Pinf ----
  sar.fit <- sar.input <- c()
  
  if (type == "Symptoms") {
    # Symptoms model uses q = 1-exp(-m) and household transform
    for (k in seq_len(xdata$n_households)) {
      q <- 1 - exp(-ch$m)
      Pinf = q
      # hh.size <- xdata$size_households[k]
      # Pinf <- 1 - (1 - q / (hh.size - 1)^ch$alpha)^(hh.size - 1)
      # sar.input <- c(sar.input, 1 - (1 - xdata$all.HH$SAR[k] / (hh.size - 1)^alpha.input)^(hh.size - 1))
      
      sar.input <- c(sar.input,  xdata$all.HH$SAR[k] )
      
      sar.fit   <- c(sar.fit,   Pinf)
    }
  } else {
    # Quantitative/virological model uses VL -> pinf0 -> hh SAR
    index.of.index <- which(xdata$col_infected == 1)
    for (id.household in seq_len(xdata$n_households)) {
      id <- index.of.index[id.household]
      hh.size <- xdata$size_households[id.household]
      
      VL <- get_VL_model_1(ch, xdata, id)
      
      pinf0 <- pinf_from_params(VL$T_P, VL$T_U, VL$V_P, VL$m, model = fit.epid.model, h = cfg$h)
      
      # subsample for speed (as your script)
      if (nrow(VL) > 500) VL <- VL[sample.int(nrow(VL), 500), ]
      
      sar.input <- c(sar.input,  xdata$all.HH$SAR[id.household]  )
      sar.fit   <- c(sar.fit,   pinf0)
    }
  }
  
  # ---- Parameter errors (only for non-symptoms) ----
  T_P.error <- T_U.error <- V_P.error <- c()
  if (type != "Symptoms") {
    row_infected <- xdata$row_infected
    col_infected <- xdata$col_infected
    T_P.error <- numeric(xdata$n_infected_total)
    T_U.error <- numeric(xdata$n_infected_total)
    V_P.error <- numeric(xdata$n_infected_total)
    
    for (k in seq_len(xdata$n_infected_total)) {
      T_P.error[k] <- mean(ch$T_P[, row_infected[k], col_infected[k]]) - xdata$all.HH$T_P[k]
      T_U.error[k] <- mean(ch$T_U[, row_infected[k], col_infected[k]]) - xdata$all.HH$T_U[k]
      V_P.error[k] <- mean(ch$V_P[, row_infected[k], col_infected[k]]) - xdata$all.HH$V_P[k]
    }
  }
  
  # ---- Assemble outputs (keep your names for downstream code) ----
  out <- data.frame(
    type = type,
    
    mean.good.number.generations = M1,
    mean.good.reconstruction     = M2,
    
    bias.alpha.mean  = bias.alpha.mean,
    bias.alpha.lower = bias.alpha.lower,
    bias.alpha.upper = bias.alpha.upper,
    
    M1.GI   = mean(fit.gi, na.rm = TRUE),
    M1.DS   = mean(fit.incubation, na.rm = TRUE),
    M1.Pinf = mean(sar.fit, na.rm = TRUE),
    
    M2.GI   = var(fit.gi, na.rm = TRUE),
    M2.DS   = var(fit.incubation, na.rm = TRUE),
    M2.Pinf = var(sar.fit, na.rm = TRUE),
    
    coverage.GI   = hdi_contains(fit.gi, parameters$GI.input),
    coverage.DS   = hdi_contains(fit.incubation, parameters$DS.input),
    coverage.Pinf = hdi_contains(sar.fit, parameters$pinf.input),
    
    T_P.error.mean  = if (length(T_P.error)) mean(T_P.error, na.rm = TRUE) else NA_real_,
    T_P.error.lower = if (length(T_P.error)) q025(T_P.error) else NA_real_,
    T_P.error.upper = if (length(T_P.error)) q975(T_P.error) else NA_real_,
    
    T_U.error.mean  = if (length(T_U.error)) mean(T_U.error, na.rm = TRUE) else NA_real_,
    T_U.error.lower = if (length(T_U.error)) q025(T_U.error) else NA_real_,
    T_U.error.upper = if (length(T_U.error)) q975(T_U.error) else NA_real_,
    
    V_P.error.mean  = if (length(V_P.error)) mean(V_P.error, na.rm = TRUE) else NA_real_,
    V_P.error.lower = if (length(V_P.error)) q025(V_P.error) else NA_real_,
    V_P.error.upper = if (length(V_P.error)) q975(V_P.error) else NA_real_
  )
  
  out
}

# ======================================================================
# RUN GRID
# ======================================================================

fit.epid.model <- cfg$fit_epid_model
alpha <- 1
vp <- 6

gi <- 5
ds <- 5
n.h.values <- c(100)
sar.values <- c(0.15, 0.3, 0.45)
I.values <- 1:100

# Offsets for plotting (used multiple places)
pinf_offset <- data.frame(
  pinf.input = factor(c("0.15","0.3","0.45"), levels = c("0.15","0.3","0.45")),
  offset = c(0.2, 0.1, 0)
)

type_pos <- data.frame(
  type = c("Virological+Clinical","Qualitative", "Clinical"),
  type_numeric = c(3,2, 1)
)


m_map <- tibble::tibble(
  pinf.input = factor(c("0.15","0.3","0.45"), levels = c("0.15","0.3","0.45")),
  m = c(1.4e-05, 3.082862e-05, 5.3e-05)
)

compute_reference_from_sim <- function(parameters, m_map,
                                       N = 20000,
                                       t_grid = seq(0, 35, by = 0.05),
                                       fixed_param = FALSE) {
  
  ref_list <- vector("list", nrow(m_map))
  
  for (i in seq_len(nrow(m_map))) {
    p <- parameters
    p$m <- m_map$m[i]
    p$fixed_param <- fixed_param
    
    sim <- simulate_transmissions_input(
      parameters = p,
      N = N,
      t_grid = t_grid
    )
    
    GI <- sim$summary$t_weighted
    SAR <- sim$summary$P
    
    ref_list[[i]] <- data.frame(
      pinf.input = m_map$pinf.input[i],
      gi_mean = mean(GI, na.rm = TRUE),
      gi_median = median(GI, na.rm = TRUE),
      gi_sd   = stats::sd(GI, na.rm = TRUE),
      pinf_mean = mean(SAR, na.rm = TRUE),
      pinf_median = median(SAR, na.rm = TRUE),
      pinf_sd   = stats::sd(SAR, na.rm = TRUE)
    )
  }
  
  dplyr::bind_rows(ref_list)
}

rows <- list()
krow <- 0L

for (n.h in n.h.values) {
  for (sar in sar.values) {
    for (I in I.values) {

      # ----------------------------------------------------------------
      # Quantitative viral load, no symptoms for infected contacts
      # Main Figure 3 + Supplementary comparator
      # ----------------------------------------------------------------
      fit.quant_no_contact_symptoms <- load_fit_objects(
        model = "quantitative_vl_no_contact_symptoms",
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n.h,
        index = I
      )

      xdata <- fit.quant_no_contact_symptoms$xdata
      parameters <- fit.quant_no_contact_symptoms$parameters
      Chains.quant_no_contact_symptoms <- fit.quant_no_contact_symptoms$Chains
      CT.quant_no_contact_symptoms <- fit.quant_no_contact_symptoms$CT

      # Same simulated data, but symptom information is removed for
      # infected contacts. Index cases keep their original status.
      contact_ids <- which(xdata$col_infected > 1)
      xdata_no_contact_symptoms <- xdata
      xdata_no_contact_symptoms$symptomatic[contact_ids] <- 0

      # ----------------------------------------------------------------
      # Quantitative viral load, with symptoms for infected contacts
      # Used in the Supplementary figure
      # ----------------------------------------------------------------
      fit.quant_contact_symptoms <- load_fit_objects(
        model = "quantitative_vl_contact_symptoms",
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n.h,
        index = I
      )

      Chains.quant_contact_symptoms <- fit.quant_contact_symptoms$Chains
      CT.quant_contact_symptoms <- fit.quant_contact_symptoms$CT

      # ----------------------------------------------------------------
      # Qualitative viral load, no symptoms for infected contacts
      # Main Figure 3
      # ----------------------------------------------------------------
      fit.qual_no_contact_symptoms <- load_fit_objects(
        model = "qualitative_vl_no_contact_symptoms",
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n.h,
        index = I
      )

      Chains.qual_no_contact_symptoms <- fit.qual_no_contact_symptoms$Chains
      CT.qual_no_contact_symptoms <- fit.qual_no_contact_symptoms$CT

      # ----------------------------------------------------------------
      # Symptoms only
      # Main Figure 3
      # ----------------------------------------------------------------
      fit.symptoms_only <- load_fit_objects(
        model = "symptoms_only",
        gi = gi,
        ds = ds,
        sar = sar,
        alpha = alpha,
        n_hh = n.h,
        index = I
      )

      Chains.symptoms_only <- fit.symptoms_only$Chains
      CT.symptoms_only <- fit.symptoms_only$CT

      # Ensure input references exist for coverage calls
      parameters$GI.input   <- gi
      parameters$DS.input   <- ds
      parameters$pinf.input <- sar

      d.quant_no_contact_symptoms <- compare_fit_input_model(
        xdata = xdata_no_contact_symptoms,
        CT_input = CT.quant_no_contact_symptoms,
        ch = Chains.quant_no_contact_symptoms,
        parameters = parameters,
        type = "Quantitative VL",
        fit.epid.model = fit.epid.model
      )

      d.quant_contact_symptoms <- compare_fit_input_model(
        xdata = xdata,
        CT_input = CT.quant_contact_symptoms,
        ch = Chains.quant_contact_symptoms,
        parameters = parameters,
        type = "Quantitative VL + contact symptoms",
        fit.epid.model = fit.epid.model
      )

      d.qual_no_contact_symptoms <- compare_fit_input_model(
        xdata = xdata_no_contact_symptoms,
        CT_input = CT.qual_no_contact_symptoms,
        ch = Chains.qual_no_contact_symptoms,
        parameters = parameters,
        type = "Qualitative VL",
        fit.epid.model = fit.epid.model
      )

      d.symptoms_only <- compare_fit_input_model(
        xdata = xdata,
        CT_input = CT.symptoms_only,
        ch = Chains.symptoms_only,
        parameters = parameters,
        type = "Symptoms",
        fit.epid.model = fit.epid.model
      )

      df0 <- data.frame(gi = gi, ds = ds, pinf.input = sar, n.h = n.h, I = I)

      krow <- krow + 1L
      rows[[krow]] <- bind_rows(
        bind_cols(df0, d.quant_no_contact_symptoms) %>%
          mutate(model = "quantitative_vl_no_contact_symptoms",
                 analysis_set = "main_and_supplementary"),

        bind_cols(df0, d.quant_contact_symptoms) %>%
          mutate(model = "quantitative_vl_contact_symptoms",
                 analysis_set = "supplementary"),

        bind_cols(df0, d.qual_no_contact_symptoms) %>%
          mutate(model = "qualitative_vl_no_contact_symptoms",
                 analysis_set = "main"),

        bind_cols(df0, d.symptoms_only) %>%
          mutate(model = "symptoms_only",
                 analysis_set = "main")
      )
    }
  }
}

df <- bind_rows(rows)
# Input ----

ref_sim <- compute_reference_from_sim(
  parameters = parameters,
  m_map = m_map,
  N = 20000,
  t_grid = seq(0, 35, by = 0.05),
  fixed_param = FALSE
)

# Objects expected by the plotting code:
mean_input_gi <- ref_sim %>% dplyr::transmute(pinf.input, input = gi_mean)
median_input_gi <- ref_sim %>% transmute(pinf.input, input = gi_median)
sd_input_gi   <- ref_sim %>% transmute(pinf.input, sdGI = gi_sd)

mean_input_pinf <- ref_sim %>% transmute(pinf.input, input = pinf_mean)
median_input_pinf <- ref_sim %>% transmute(pinf.input, input  = pinf_median)
sd_input_pinf   <- ref_sim %>% transmute(pinf.input, sdGI = pinf_sd)

# ---- Incubation SD: independent of pinf ----
n.samples <- 20000
incubation_input_distribution <- parameters$Delta_S_star * exp(rnorm(n.samples, sd = parameters$Delta_S_eta))
sdDS_true <- stats::sd(incubation_input_distribution, na.rm = TRUE)
mean_input_incubation <- m_map %>% transmute( pinf.input, incubation = mean(incubation_input_distribution))
median_input_incubation <- m_map %>% transmute(pinf.input, incubation = median(incubation_input_distribution))
sd_input_incubation <- m_map %>% transmute(pinf.input, sdDS = sdDS_true)


# ======================================================================
# FIGURES — shared plotting helpers
# ======================================================================

prepare_plot_df <- function(df_plot, type_levels) {
  df_plot %>%
    mutate(
      type = dplyr::recode(type, "Symptoms" = "Symptoms only"),
      type = factor(type, levels = type_levels, ordered = TRUE),
      pinf.input = factor(pinf.input, levels = c("0.15", "0.3", "0.45")),
      n.h = factor(n.h, levels = as.character(n.h.values), ordered = TRUE)
    )
}

make_type_pos <- function(type_levels) {
  data.frame(
    type = type_levels,
    type_numeric = rev(seq_along(type_levels))
  )
}

summarise_metric <- function(df2, variable) {
  df2 %>%
    group_by(type, pinf.input, n.h) %>%
    dplyr::summarise(
      Mean  = mean(.data[[variable]], na.rm = TRUE),
      lower = q025(.data[[variable]]),
      upper = q975(.data[[variable]]),
      .groups = "drop"
    )
}

summarise_sd_metric <- function(df2, variable) {
  df2 %>%
    group_by(type, pinf.input, n.h) %>%
    dplyr::summarise(
      Mean  = median(sqrt(.data[[variable]]), na.rm = TRUE),
      lower = q025(sqrt(.data[[variable]])),
      upper = q975(sqrt(.data[[variable]])),
      .groups = "drop"
    )
}

plot_estimate_panel <- function(dat, reference_df, reference_col,
                                xlab, type_pos_local, xlim_values = NULL) {
  p <- dat %>%
    left_join(reference_df, by = "pinf.input") %>%
    left_join(type_pos_local, by = "type") %>%
    left_join(pinf_offset, by = "pinf.input") %>%
    mutate(y_adjusted = type_numeric + offset) %>%
    ggplot(aes(x = Mean, y = y_adjusted, color = pinf.input)) +
    geom_segment(aes(x = lower, xend = upper, y = y_adjusted, yend = y_adjusted)) +
    geom_point(aes(x = .data[[reference_col]], y = y_adjusted),
               size = 2, stroke = 1, color = "red", shape = 5) +
    geom_point(size = 2, stroke = 1) +
    scale_y_continuous(breaks = type_pos_local$type_numeric,
                       labels = type_pos_local$type) +
    labs(x = xlab, y = "Data type", color = "") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.text.y = element_text(size = 0),
      legend.position = "none"
    ) +
    scale_color_brewer(palette = "Dark2")

  if (!is.null(xlim_values)) {
    p <- p + xlim(xlim_values)
  }

  p
}

plot_chain_panel <- function(dat, xlab, type_pos_local) {
  dat %>%
    left_join(type_pos_local, by = "type") %>%
    left_join(pinf_offset, by = "pinf.input") %>%
    mutate(y_adjusted = type_numeric + offset) %>%
    ggplot(aes(x = Mean, y = y_adjusted, color = pinf.input)) +
    geom_point(size = 2, stroke = 1) +
    geom_segment(aes(x = lower, xend = upper, y = y_adjusted, yend = y_adjusted)) +
    scale_y_continuous(breaks = type_pos_local$type_numeric,
                       labels = type_pos_local$type) +
    labs(x = xlab, y = "Data type", color = "") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.text.y = element_text(size = 0),
      legend.position = "none"
    ) +
    scale_color_brewer(palette = "Dark2") +
    xlim(c(0, 100))
}

make_metric_plots <- function(df_plot, type_levels) {
  df2 <- prepare_plot_df(df_plot, type_levels)
  type_pos_local <- make_type_pos(type_levels)

  df.GI <- summarise_metric(df2, "M1.GI")
  df.DS <- summarise_metric(df2, "M1.DS")
  df.Pinf <- summarise_metric(df2, "M1.Pinf")

  df.GI.M2 <- summarise_sd_metric(df2, "M2.GI")
  df.DS.M2 <- summarise_sd_metric(df2, "M2.DS")
  df.Pinf.M2 <- summarise_sd_metric(df2, "M2.Pinf")

  df.Nb.Generations <- df2 %>%
    filter(!is.na(mean.good.number.generations)) %>%
    group_by(type, pinf.input, n.h) %>%
    dplyr::summarise(
      Mean  = 100 * mean(mean.good.number.generations, na.rm = TRUE),
      lower = 100 * q025(mean.good.number.generations),
      upper = 100 * q975(mean.good.number.generations),
      .groups = "drop"
    )

  df.Good.Reconstructions <- df2 %>%
    filter(!is.na(mean.good.reconstruction)) %>%
    group_by(type, pinf.input, n.h) %>%
    dplyr::summarise(
      Mean  = 100 * mean(mean.good.reconstruction, na.rm = TRUE),
      lower = 100 * q025(mean.good.reconstruction),
      upper = 100 * q975(mean.good.reconstruction),
      .groups = "drop"
    )

  list(
    G.GI = plot_estimate_panel(
      df.GI,
      reference_df = mean_input_gi,
      reference_col = "input",
      xlab = "Mean Generation interval (days)",
      type_pos_local = type_pos_local
    ),
    G.DS = plot_estimate_panel(
      df.DS,
      reference_df = mean_input_incubation,
      reference_col = "incubation",
      xlab = "Mean Incubation period (days)",
      type_pos_local = type_pos_local
    ),
    G.Pinf = plot_estimate_panel(
      df.Pinf,
      reference_df = mean_input_pinf,
      reference_col = "input",
      xlab = "Mean transmission probability",
      type_pos_local = type_pos_local,
      xlim_values = c(0, NA)
    ),
    G.GI.M2 = plot_estimate_panel(
      df.GI.M2,
      reference_df = sd_input_gi,
      reference_col = "sdGI",
      xlab = "sd Generation interval (days)",
      type_pos_local = type_pos_local,
      xlim_values = c(0, NA)
    ),
    G.DS.M2 = plot_estimate_panel(
      df.DS.M2,
      reference_df = sd_input_incubation,
      reference_col = "sdDS",
      xlab = "sd Incubation period (days)",
      type_pos_local = type_pos_local
    ),
    G.Pinf.M2 = plot_estimate_panel(
      df.Pinf.M2,
      reference_df = sd_input_pinf,
      reference_col = "sdGI",
      xlab = "sd Transmission probability",
      type_pos_local = type_pos_local
    ),
    G.Good.Reconstructions = plot_chain_panel(
      df.Good.Reconstructions,
      xlab = "Proportion of correct transmission chains (%)",
      type_pos_local = type_pos_local
    ),
    G.Nb.Generations = plot_chain_panel(
      df.Nb.Generations,
      xlab = "Proportion of correct number of generations (%)",
      type_pos_local = type_pos_local
    )
  )
}

# ======================================================================
# EXPORTS
# ======================================================================
A <- 1.1

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

export_plot <- function(plot, filename_base, width, height) {
  ggpubr::ggexport(
    plot,
    filename = file.path(cfg$out_dir, paste0(filename_base, ".pdf")),
    width = width,
    height = height
  )
}

# ======================================================================
# MAIN FIGURE 3
# Uses the three main designs:
#   Quantitative VL, no symptoms for infected contacts
#   Qualitative VL, no symptoms for infected contacts
#   Symptoms only
# ======================================================================

df_main <- df %>%
  filter(model %in% c(
    "quantitative_vl_no_contact_symptoms",
    "qualitative_vl_no_contact_symptoms",
    "symptoms_only"
  ))

main_type_levels <- c(
  "Quantitative VL",
  "Qualitative VL",
  "Symptoms only"
)

main_plots <- make_metric_plots(df_main, main_type_levels)

G.Good.Reconstructions.legend <- main_plots$G.Good.Reconstructions +
  theme(legend.position = "right")

export_plot(
  G.Good.Reconstructions.legend,
  filename_base = "Figure3_legend",
  width = 6 * 0.9,
  height = 2 * 0.9
)

Figure3 <- ggarrange(
  main_plots$G.GI + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  main_plots$G.GI.M2,
  main_plots$G.DS + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  main_plots$G.DS.M2,
  main_plots$G.Pinf + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  main_plots$G.Pinf.M2,
  main_plots$G.Good.Reconstructions,
  main_plots$G.Nb.Generations,
  labels = c("A", "D", "B", "E", "C", "F", "G", "H"),
  nrow = 4,
  ncol = 2
)

print(Figure3)

export_plot(
  Figure3,
  filename_base = "Figure3",
  width = 8 * A,
  height = 8 * A
)

# ======================================================================
# SUPPLEMENTARY FIGURE
# Compares quantitative VL with versus without symptoms for contacts.
# Both are loaded within the same simulation loop above.
# ======================================================================

df_supplementary <- df %>%
  filter(model %in% c(
    "quantitative_vl_no_contact_symptoms",
    "quantitative_vl_contact_symptoms"
  ))

supplementary_type_levels <- c(
  "Quantitative VL",
  "Quantitative VL + contact symptoms"
)

supplementary_plots <- make_metric_plots(df_supplementary, supplementary_type_levels)

Figure3_supplementary_metrics <- ggarrange(
  supplementary_plots$G.GI + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  supplementary_plots$G.DS + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  supplementary_plots$G.Pinf + theme(strip.background = element_blank(), axis.text.y = element_text(size = 0)),
  supplementary_plots$G.Good.Reconstructions,
  supplementary_plots$G.GI.M2,
  supplementary_plots$G.DS.M2,
  supplementary_plots$G.Pinf.M2,
  supplementary_plots$G.Nb.Generations,
  labels = c("A", "B", "C", "G", "D", "E", "F", "H"),
  nrow = 2,
  ncol = 4
)

print(Figure3_supplementary_metrics)

export_plot(
  Figure3_supplementary_metrics,
  filename_base = "Figure3_supplementary_metrics",
  width = 10 * A,
  height = 5 * A
)
