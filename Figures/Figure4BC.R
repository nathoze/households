# ======================================================================
# Figure 4BC — Time to reach 95% cumulative infectiousness
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(here)
})

source(here::here("src", "paths.R"))
source(here::here("src", "models", "functions_household_transmission.R"))
source(here::here("src", "main.R"))
source(here::here("src", "visualization", "MCMC_predictions.R"))
source(here::here("src", "models", "transmission probability conditional.R"))
source(here::here("src", "helpers.R"))

cfg <- list(
  h = 0.5,
  out_dir = here::here("Figures")
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

h <- cfg$h
alpha <- 1
gi <- 5
ds <- 5
n.h.values <- c(100)
sar.values <- c(0.15, 0.30, 0.45)
I.values <- 1:100
t_grid <- seq(0, 30, by = 0.05)

m_map <- tibble(
  pinf.input = factor(c("0.15", "0.30", "0.45"), levels = c("0.15", "0.30", "0.45")),
  m_true = c(1.4e-05, 3.082862e-05, 5.3e-05)
)

pinf_offset <- tibble(
  pinf.input = factor(c("0.15", "0.30", "0.45"), levels = c("0.15", "0.30", "0.45")),
  offset = c(0.20, 0.10, 0.00)
)

type_pos <- tibble(
  type = c("Quantitative VL", "Qualitative VL", "Symptoms only"),
  type_numeric = c(3, 2, 1)
)

models_to_load <- tibble(
  type = c("Quantitative VL", "Qualitative VL", "Symptoms only"),
  model = c(
    "quantitative_vl_no_contact_symptoms",
    "qualitative_vl_no_contact_symptoms",
    "symptoms_only"
  )
)

q025 <- function(x) unname(quantile(x, 0.025, na.rm = TRUE))
q975 <- function(x) unname(quantile(x, 0.975, na.rm = TRUE))

resolve_model_name <- function(model) {
  aliases <- c(
    quantitative_vl_with_contact_symptoms = "quantitative_vl_contact_symptoms"
  )
  if (model %in% names(aliases)) return(unname(aliases[[model]]))
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

get_t95_infectiousness <- function(time, lambda) {
  ok <- is.finite(time) & is.finite(lambda)
  time <- time[ok]
  lambda <- lambda[ok]

  if (length(time) < 2L || sum(lambda, na.rm = TRUE) <= 0) return(NA_real_)

  dt <- c(0, diff(time))
  cumulative <- cumsum(lambda * dt)
  cumulative <- cumulative / max(cumulative, na.rm = TRUE)

  idx <- which(cumulative >= 0.95)[1]
  if (is.na(idx)) return(NA_real_)
  time[idx]
}

get_VL_model_1 <- function(Chains, xdata, id) {
  data.frame(
    T_P = exp(Chains$log_T_P_mu + Chains$log_T_P_sd * Chains$T_P_eta[, id]) * xdata$T_P_star,
    T_U = exp(Chains$log_T_U_mu + Chains$log_T_U_sd * Chains$T_U_eta[, id]) * xdata$T_U_star,
    V_P = exp(Chains$log_V_P_mu + Chains$log_V_P_sd * Chains$V_P_eta[, id]) * xdata$V_P_star,
    m = 10^(Chains$logm)
  )
}

get_t95_from_params <- function(T_P, T_U, V_P, m, t_grid, h) {
  VL_t <- model_gen_triangle(
    T_P = T_P,
    T_U = T_U,
    V_P = V_P,
    t = t_grid,
    tinfection = 0,
    normalize.time = FALSE
  )

  lambda_t <- m * 10^(h * VL_t$logV)
  get_t95_infectiousness(time = VL_t$t, lambda = lambda_t)
}

get_t95_posterior_median <- function(Chains, xdata, id, t_grid, h) {
  VL <- get_VL_model_1(Chains, xdata, id)

  get_t95_from_params(
    T_P = median(VL$T_P, na.rm = TRUE),
    T_U = median(VL$T_U, na.rm = TRUE),
    V_P = median(VL$V_P, na.rm = TRUE),
    m = median(VL$m, na.rm = TRUE),
    t_grid = t_grid,
    h = h
  )
}

get_t95_gamma_infectiousness <- function(shape, inverse_scale, t_grid) {
  lambda_t <- dgamma(t_grid, shape = shape, rate = inverse_scale)
  get_t95_infectiousness(time = t_grid, lambda = lambda_t)
}

get_t95_symptom_posterior_median <- function(Chains, t_grid) {
  shape_med <- median(Chains$gamma_infection_shape, na.rm = TRUE)
  get_t95_gamma_infectiousness(shape = shape_med, inverse_scale = 0.3, t_grid = t_grid)
}

get_param_sd <- function(parameters, base_name) {
  candidates <- c(
    paste0(base_name, "_eta"),
    paste0(base_name, ".eta"),
    paste0(base_name, "_sd"),
    paste0(base_name, ".sd")
  )

  for (nm in candidates) {
    if (!is.null(parameters[[nm]])) return(parameters[[nm]])
  }

  stop("Could not find an input SD parameter for ", base_name)
}

sample_input_vl_parameters <- function(parameters, n_samples) {
  T_P_sd <- get_param_sd(parameters, "T_P")
  T_U_sd <- get_param_sd(parameters, "T_U")
  V_P_sd <- get_param_sd(parameters, "V_P")

  tibble(
    T_P = parameters$T_P_star * exp(rnorm(n_samples, mean = 0, sd = T_P_sd)),
    T_U = parameters$T_U_star * exp(rnorm(n_samples, mean = 0, sd = T_U_sd)),
    V_P = parameters$V_P_star * exp(rnorm(n_samples, mean = 0, sd = V_P_sd))
  )
}

compute_input_t95_reference <- function(parameters, m_map, n_samples = 20000, t_grid, h) {
  sampled_vl <- sample_input_vl_parameters(parameters, n_samples)
  ref_list <- vector("list", nrow(m_map))

  for (i in seq_len(nrow(m_map))) {
    m_i <- m_map$m_true[i]

    t95_i <- vapply(
      seq_len(nrow(sampled_vl)),
      function(j) {
        get_t95_from_params(
          T_P = sampled_vl$T_P[j],
          T_U = sampled_vl$T_U[j],
          V_P = sampled_vl$V_P[j],
          m = m_i,
          t_grid = t_grid,
          h = h
        )
      },
      numeric(1)
    )

    ref_list[[i]] <- tibble(
      pinf.input = m_map$pinf.input[i],
      input.mean.t95 = mean(t95_i, na.rm = TRUE),
      input.var.t95 = var(t95_i, na.rm = TRUE)
    )
  }

  bind_rows(ref_list) %>%
    mutate(pinf.input = factor(pinf.input, levels = levels(m_map$pinf.input)))
}

safe_var <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(0)
  var(x)
}

rows <- list()
krow <- 0L
last_parameters <- NULL

for (n.h in n.h.values) {
  for (sar in sar.values) {
    for (I in I.values) {
      for (j in seq_len(nrow(models_to_load))) {
        type_name <- models_to_load$type[j]
        model_name <- models_to_load$model[j]

        fit <- load_fit_objects(
          model = model_name,
          gi = gi,
          ds = ds,
          sar = sar,
          alpha = alpha,
          n_hh = n.h,
          index = I
        )

        xdata <- fit$xdata
        parameters <- fit$parameters
        Chains <- fit$Chains
        last_parameters <- parameters

        if (model_name == "symptoms_only") {
          krow <- krow + 1L
          rows[[krow]] <- tibble(
            type = type_name,
            model = model_name,
            pinf.input = sprintf("%.2f", sar),
            n.h = n.h,
            I = I,
            id = NA_integer_,
            t95 = get_t95_symptom_posterior_median(Chains = Chains, t_grid = t_grid)
          )
        } else {
          for (id in seq_len(xdata$n_infected_total)) {
            krow <- krow + 1L
            rows[[krow]] <- tibble(
              type = type_name,
              model = model_name,
              pinf.input = sprintf("%.2f", sar),
              n.h = n.h,
              I = I,
              id = id,
              t95 = get_t95_posterior_median(
                Chains = Chains,
                xdata = xdata,
                id = id,
                t_grid = t_grid,
                h = h
              )
            )
          }
        }
      }
    }
  }
}

if (length(rows) == 0L) stop("No T95 values were extracted. Check paths.R and results/ folders.")

df_t95 <- bind_rows(rows) %>%
  mutate(
    type = factor(type, levels = type_pos$type, ordered = TRUE),
    pinf.input = factor(pinf.input, levels = c("0.15", "0.30", "0.45")),
    n.h = factor(n.h, levels = as.character(n.h.values), ordered = TRUE)
  )

ref_t95 <- compute_input_t95_reference(
  parameters = last_parameters,
  m_map = m_map,
  n_samples = 20000,
  t_grid = t_grid,
  h = h
)

df_t95_summary <- df_t95 %>%
  group_by(type, pinf.input, n.h, I) %>%
  summarise(
    mean.t95 = mean(t95, na.rm = TRUE),
    var.t95 = safe_var(t95),
    .groups = "drop"
  )

df_t95_mean <- df_t95_summary %>%
  group_by(type, pinf.input, n.h) %>%
  summarise(
    Mean = mean(mean.t95, na.rm = TRUE),
    lower = q025(mean.t95),
    upper = q975(mean.t95),
    .groups = "drop"
  )

df_t95_var <- df_t95_summary %>%
  group_by(type, pinf.input, n.h) %>%
  summarise(
    Mean = mean(var.t95, na.rm = TRUE),
    lower = q025(var.t95),
    upper = q975(var.t95),
    .groups = "drop"
  )

G.T95.mean <- df_t95_mean %>%
  left_join(ref_t95, by = "pinf.input") %>%
  left_join(type_pos, by = "type") %>%
  left_join(pinf_offset, by = "pinf.input") %>%
  mutate(y_adjusted = type_numeric + offset) %>%
  ggplot(aes(x = Mean, y = y_adjusted, color = pinf.input)) +
  geom_segment(aes(x = lower, xend = upper, y = y_adjusted, yend = y_adjusted)) +
  geom_point(aes(x = input.mean.t95, y = y_adjusted), size = 2, stroke = 1, color = "red", shape = 5) +
  geom_point(size = 2, stroke = 1) +
  scale_y_continuous(breaks = type_pos$type_numeric, labels = type_pos$type) +
  labs(
    x = "Mean time to reach 95% cumulative infectiousness (days)",
    y = "Data type",
    color = ""
  ) +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.y = element_text(size = 0), legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  xlim(c(0, NA))

G.T95.var <- df_t95_var %>%
  left_join(ref_t95, by = "pinf.input") %>%
  left_join(type_pos, by = "type") %>%
  left_join(pinf_offset, by = "pinf.input") %>%
  mutate(y_adjusted = type_numeric + offset) %>%
  ggplot(aes(x = Mean, y = y_adjusted, color = pinf.input)) +
  geom_segment(aes(x = lower, xend = upper, y = y_adjusted, yend = y_adjusted)) +
  geom_point(aes(x = input.var.t95, y = y_adjusted), size = 2, stroke = 1, color = "red", shape = 5) +
  geom_point(size = 2, stroke = 1) +
  scale_y_continuous(breaks = type_pos$type_numeric, labels = type_pos$type) +
  labs(
    x = "Variance of time to reach 95% cumulative infectiousness (days^2)",
    y = "Data type",
    color = ""
  ) +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.y = element_text(size = 0), legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  xlim(c(0, NA))

Figure4BC <- ggarrange(
  G.T95.mean,
  G.T95.var,
  labels = c("B", "C"),
  nrow = 1,
  ncol = 2
)

print(Figure4BC)

A <- 0.9

ggpubr::ggexport(
  Figure4BC,
  filename = file.path(cfg$out_dir, "Figure4BC.pdf"),
  width = 6 * A,
  height = 3 * A
)

saveRDS(
  list(
    df_t95 = df_t95,
    df_t95_summary = df_t95_summary,
    df_t95_mean = df_t95_mean,
    df_t95_var = df_t95_var,
    ref_t95 = ref_t95,
    Figure4BC = Figure4BC
  ),
  file = file.path(cfg$out_dir, "Figure4BC_objects.rds")
)
