# helpers

colors <- list(
  within_host = "#20B4B6",
  fit = "#FF8A5B",
  data_exp = "#5C9B58",
  data_constant = "#EA526F",
  data_constant_fr = "#e41e43",
  data_constant_d = "#A9424F",
  data_exp_fr = "#6b8a69",
  data_exp_d = "#1C5B38",
  fr = "#DDDDDD",
  d = "#444444",
  und = "#FCEADE",
  symptoms = "#EA526F",
  full = "#5C9B58"
)
q025 <- function(x) unname(stats::quantile(x, probs = 0.025, na.rm = TRUE))
q975 <- function(x) unname(stats::quantile(x, probs = 0.975, na.rm = TRUE))

# SAR "pinf" (infection prob. per contact) pour exp/constant
pinf_from_vl <- function(T_P, T_U, V_P, m, model, h = 0.5) {
  if (model == "exp") {
    1 - exp(-m * 10^(2*h) / (h * V_P * log(10)) * (T_P + T_U) * (10^(h * V_P) - 1))
  } else if (model == "constant") {
    1 - exp(-m * (T_P + T_U))
  } else {
    stop("Unknown model: ", model)
  }
}

# SAR household-size (comme ton get_SAR_Chains_hh_size)
sar_from_pinf_hh <- function(pinf0, hh_size, alpha) {
  1 - (1 - pinf0 / (hh_size - 1)^alpha)^(hh_size - 1)
}

get_VL <- function(Chains, xdata, id) {
  data.frame(
    T_P   = exp(Chains$log_T_P_mu + Chains$log_T_P_sd * Chains$T_P_eta[, id]) * xdata$T_P_star,
    T_U   = exp(Chains$log_T_U_mu + Chains$log_T_U_sd * Chains$T_U_eta[, id]) * xdata$T_U_star,
    V_P   = exp(Chains$log_V_P_mu + Chains$log_V_P_sd * Chains$V_P_eta[, id]) * xdata$V_P_star,
    m     = 10^(Chains$logm),
    alpha = Chains$alpha
  )
}

add_SAR <- function(VL, model, h = 0.5) {
  pinf0 <- pinf_from_vl(VL$T_P, VL$T_U, VL$V_P, VL$m, model = model, h = h)
  VL %>% mutate(SAR = pinf0)
}

add_SAR_hh <- function(VL, model, hh_size, h = 0.5) {
  pinf0 <- pinf_from_vl(VL$T_P, VL$T_U, VL$V_P, VL$m, model = model, h = h)
  SAR <- sar_from_pinf_hh(pinf0, hh_size = hh_size, alpha = VL$alpha)
  VL %>% mutate(SAR = SAR)
}

save_plot <- function(p, filename, width, height, dpi = 300) {
  ggplot2::ggsave(filename = filename, plot = p, width = width, height = height, dpi = dpi)
}

quantile.times <- function(cumsumq, t, N = 20000) {
  # cumsumq should be increasing and end near 1
  Times <- numeric(N)
  
  for (n in seq_len(N)) {
    w <- which(runif(1) < cumsumq)[1]
    Times[n] <- t[w]
  }
  
  list(
    q025 = quantile025(Times),
    q05  = quantile05(Times),
    q95  = quantile95(Times),
    q975 = quantile975(Times)
  )
}
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}