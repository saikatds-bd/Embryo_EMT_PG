library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(patchwork)

gags <- readxl::read_xlsx("GAGs_in_adult_tissue.xlsx")

gags_long <- gags %>%
  pivot_longer(
    cols = 6:86,
    names_to = "Cell Type",
    values_to = "Log2Expression"
  )

subgroup_cols <- c(
  "Initiation"     = "sienna1",
  "Polymerization" = "#1D4E89",
  "Sulfation"      = "#56B4E9",
  "Epimerization"  = "blue",
  "Regulation"     = "#F4A261",
  "Modification"   = "red"
)

subgroup_order <- c(
  "Initiation",
  "Polymerization",
  "Sulfation",
  "Epimerization",
  "Regulation",
  "Modification"
)



gene_order <- gags_long %>%
  group_by(Row.names, Subgroup) %>%
  summarise(
    median_expr = median(Log2Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) %>%
  arrange(Subgroup, desc(median_expr)) %>%
  pull(Row.names)

gags_long_plot <- gags_long %>%
  mutate(
    Subgroup = factor(Subgroup, levels = subgroup_order),
    Row.names = factor(Row.names, levels = gene_order)
  )

# ------------------------------------------------------------
# Subgroup-level distribution summary
# This will be shown as a background band in each facet
# ------------------------------------------------------------

subgroup_summary <- gags_long_plot %>%
  group_by(Subgroup) %>%
  summarise(
    q10 = quantile(Log2Expression, 0.10, na.rm = TRUE),
    q25 = quantile(Log2Expression, 0.25, na.rm = TRUE),
    median = median(Log2Expression, na.rm = TRUE),
    q75 = quantile(Log2Expression, 0.75, na.rm = TRUE),
    q90 = quantile(Log2Expression, 0.90, na.rm = TRUE),
    mean = mean(Log2Expression, na.rm = TRUE),
    sd = sd(Log2Expression, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_gag_gene_boxplot <- ggplot(
  gags_long_plot,
  aes(x = Row.names, y = Log2Expression)
) +
  # Wider background range: 10th to 90th percentile
  geom_rect(
    data = subgroup_summary,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = q10,
      ymax = q90,
      fill = Subgroup
    ),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  # Stronger background band: IQR, Q1 to Q3
  geom_rect(
    data = subgroup_summary,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = q25,
      ymax = q75,
      fill = Subgroup
    ),
    inherit.aes = FALSE,
    alpha = 0.16
  ) +
  # Subgroup median line
  geom_hline(
    data = subgroup_summary,
    aes(yintercept = median),
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey25"
  ) +
  # Gene-level boxplots
  geom_boxplot(
    aes(fill = Subgroup),
    outlier.shape = NA,
    width = 0.65,
    linewidth = 0.3,
    alpha = 0.95
  ) +
  # Cell-type dots
  geom_jitter(
    color = "black",
    width = 0.18,
    size = 0.35,
    alpha = 0.45
  ) +
  facet_grid(
    . ~ Subgroup,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = subgroup_cols) +
  labs(
    x = "GAG biosynthetic enzyme",
    y = "Log2 expression",
    fill = "Subgroup"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.35
    ),
    strip.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    panel.spacing.x = unit(0.35, "lines"),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

p_gag_gene_boxplot

ggsave("Distribution of GAGs in adult tissue.pdf", height = 4, width = 8)

