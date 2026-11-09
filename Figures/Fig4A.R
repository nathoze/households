# ======================================================================
# Display: rows = Ptr, colors = study design
# ======================================================================

design_cols <- c(
  "Quantitative VL" = "#1b9e77",
  "Symptoms only"   = "#d95f02",
  "Qualitative VL"  = "#7570b3"
)

design_labels <- expression(D^{VL}, D^{TSO}, D^{PCR})

ptr_positions <- tibble(
  pinf.input = factor(
    c("0.15", "0.30", "0.45"),
    levels = c("0.15", "0.30", "0.45")
  ),
  ptr_y = c(3, 2, 1)
)

design_positions <- tibble(
  design = factor(design_levels, levels = design_levels),
  design_offset = c(0.18, 0, -0.18)
)

df_m <- bind_rows(m_rows) %>%
  mutate(
    design = factor(design, levels = design_levels),
    pinf.input = factor(
      pinf.input,
      levels = levels(ptr_positions$pinf.input)
    )
  ) %>%
  left_join(ptr_positions, by = "pinf.input") %>%
  left_join(design_positions, by = "design") %>%
  mutate(y_design = ptr_y + design_offset)

# One input reference per Ptr, shared by all three designs.
input_ticks <- m_map %>%
  left_join(ptr_positions, by = "pinf.input")

# ======================================================================
# Plot
# ======================================================================

G.m <- ggplot(
  df_m,
  aes(
    x = m_median,
    y = y_design,
    fill = design,
    color = design,
    group = interaction(pinf.input, design)
  )
) +
  geom_boxplot(
    orientation = "y",
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_point(
    size = 1.1,
    alpha = 0.35,
    position = position_jitter(
      width = 0,
      height = 0.035,
      seed = 42
    ),
    # Include every design in the legend, including Symptoms only.
    show.legend = TRUE
  ) +
  geom_segment(
    data = input_ticks,
    aes(
      x = m_true,
      xend = m_true,
      y = ptr_y - 0.29,
      yend = ptr_y + 0.29
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  scale_x_log10(labels = label_scientific()) +
  scale_y_continuous(
    breaks = c(3, 2, 1),
    labels = c("0.15", "0.30", "0.45"),
    limits = c(0.5, 3.5),
    minor_breaks = NULL
  ) +
  scale_fill_manual(
    values = design_cols,
    limits = design_levels,
    guide = "none"
  ) +
  scale_color_manual(
    values = design_cols,
    limits = design_levels,
    breaks = design_levels,
    labels = design_labels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 2.5, shape = 16)
    )
  ) +
  labs(
    x = expression("Estimated transmission parameter " * m),
    y = expression(P[tr]),
    color = "Study design"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(
      size = 16,
      angle = 0,
      vjust = 0.5
    ),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

print(G.m)

# ======================================================================
# Export
# ======================================================================

ggsave(
  file.path(cfg$out_dir, "Figure4A.pdf"),
  G.m,
  width = 5,
  height = 6
)

ggsave(
  file.path(cfg$out_dir, "Figure4A.png"),
  G.m,
  width = 5,
  height = 6,
  dpi = cfg$res_dpi
)