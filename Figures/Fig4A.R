# ======================================================================
# Figure 4A — Recovery of transmission parameter m across study designs
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(ggpubr)
  library(here)
})

source(here::here("src", "paths.R"))

# ---- Config ----
cfg <- list(
  out_dir = here::here("Figures"),
  res_dpi = 300
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Simulation grid ----
gi <- 5
ds <- 5
alpha <- 1
n.h.values <- c(100)
sar.values <- c(0.15, 0.30, 0.45)
I.values <- 1:100

# ---- True input values of m ----
m_map <- tibble(
  pinf.input = factor(
    c("0.15", "0.30", "0.45"),
    levels = c("0.15", "0.30", "0.45")
  ),
  m_true = c(1.4e-05, 3.082862e-05, 5.3e-05)
)

# ---- Models to load ----
# Figure 4A compares the two viral-load study designs that estimate m:
#   1. quantitative viral load, no symptom information for contacts
#   2. qualitative viral load, no symptom information for contacts
#
# "Symptoms only" is kept as an empty reference row because this model
# does not estimate the viral-load transmission parameter m.

fits_to_load <- tibble(
  design = c(
    "Quantitative VL",
    "Qualitative VL"
  ),
  model = c(
    "quantitative_vl_no_contact_symptoms",
    "qualitative_vl_no_contact_symptoms"
  )
)

design_levels <- c(
  "Quantitative VL",
  "Qualitative VL",
  "Symptoms only"
)

# Optional alias if another script still uses this variant.
resolve_model_name <- function(model) {
  aliases <- c(
    quantitative_vl_with_contact_symptoms = "quantitative_vl_contact_symptoms"
  )

  if (model %in% names(aliases)) {
    return(unname(aliases[[model]]))
  }

  model
}

load_chains_object <- function(model, gi, ds, sar, alpha, n_hh, index) {
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

  if (!file.exists(paths$chains_path)) {
    stop(
      "Missing chains file for model '", model, "':\n",
      paths$chains_path,
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  load(paths$chains_path, envir = env)

  if (!exists("Chains", envir = env)) {
    stop("Object 'Chains' not found in: ", paths$chains_path, call. = FALSE)
  }

  env$Chains
}

# ======================================================================
# Extract posterior summaries of m
# ======================================================================

m_rows <- list()
k <- 0L

for (n.h in n.h.values) {
  for (sar in sar.values) {
    for (I in I.values) {

      for (j in seq_len(nrow(fits_to_load))) {

        design_name <- fits_to_load$design[j]
        model_name <- fits_to_load$model[j]

        Chains <- load_chains_object(
          model = model_name,
          gi = gi,
          ds = ds,
          sar = sar,
          alpha = alpha,
          n_hh = n.h,
          index = I
        )

        m_post <- if ("logm" %in% names(Chains)) {
          10^(Chains$logm)
        } else if ("m" %in% names(Chains)) {
          Chains$m
        } else {
          warning("No m or logm found for model: ", model_name)
          next
        }

        k <- k + 1L
        m_rows[[k]] <- tibble(
          design = design_name,
          model = model_name,
          pinf.input = sprintf("%.2f", sar),
          gi = gi,
          ds = ds,
          n.h = n.h,
          I = I,
          m_mean = mean(m_post, na.rm = TRUE),
          m_median = median(m_post, na.rm = TRUE),
          m_lower = quantile(m_post, 0.025, na.rm = TRUE),
          m_upper = quantile(m_post, 0.975, na.rm = TRUE)
        )
      }
    }
  }
}

if (length(m_rows) == 0L) {
  stop("No posterior m summaries were extracted. Check paths.R and results/ folders.")
}

# Same colors as previous figure
pinf_cols <- c(
  "0.15" = "#1b9e77",
  "0.30" = "#d95f02",
  "0.45" = "#7570b3"
)

df_m <- bind_rows(m_rows) %>%
  mutate(
    design = factor(design, levels = design_levels),
    pinf.input = factor(
      pinf.input,
      levels = c("0.15", "0.30", "0.45")
    )
  )

# Add blank Symptoms-only level explicitly
df_blank <- expand.grid(
  design = factor("Symptoms only", levels = design_levels),
  pinf.input = factor(
    c("0.15", "0.30", "0.45"),
    levels = levels(df_m$pinf.input)
  )
)

# True input points for all study-design rows
input_points <- expand.grid(
  design = design_levels,
  pinf.input = c("0.15", "0.30", "0.45")
) %>%
  as_tibble() %>%
  mutate(
    design = factor(design, levels = design_levels),
    pinf.input = factor(pinf.input, levels = levels(df_m$pinf.input))
  ) %>%
  left_join(m_map, by = "pinf.input")

# ======================================================================
# Plot
# ======================================================================

G.m <- ggplot(df_m, aes(x = design, y = m_median, fill = pinf.input)) +

  geom_blank(
    data = df_blank,
    aes(x = design, y = NA_real_, fill = pinf.input),
    inherit.aes = FALSE
  ) +

  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +

  geom_jitter(
    aes(color = pinf.input),
    size = 1.1,
    alpha = 0.35,
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    )
  ) +

  geom_point(
    data = input_points,
    aes(x = design, y = m_true, group = pinf.input),
    inherit.aes = FALSE,
    color = "red",
    shape = 5,
    size = 2,
    position = position_dodge(width = 0.75)
  ) +

  scale_x_discrete(limits = rev(design_levels)) +
  scale_fill_manual(values = pinf_cols) +
  scale_color_manual(values = pinf_cols) +
  scale_y_log10(labels = label_scientific()) +

  labs(
    x = "Study design",
    y = expression("Estimated transmission parameter " * m),
    fill = expression(P[tr]),
    color = expression(P[tr])
  ) +

  coord_flip() +

  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

print(G.m)

# ======================================================================
# Export
# ======================================================================

A <- 0.9

ggpubr::ggexport(
  G.m,
  filename = file.path(cfg$out_dir, "Figure4A.pdf"),
  width = 5 * A,
  height = 3 * A
)
