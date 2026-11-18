# Panel B: CD4 CTL proportion, R vs NR, two-sided Mann-Whitney U; bar, violin and dot versions

source("scripts/00_setup.R")
library(ggpubr)
library(ggbeeswarm)

samples <- read_csv(file.path(dir_tab, "risk_samples.csv"), show_col_types = FALSE) %>%
  filter(sample_id %in% read_csv(file.path(dir_tab, paste0("bayesprism_fractions", run_tag, ".csv")),
                                 col_select = sample_id, show_col_types = FALSE)$sample_id)
n <- table(samples$response)

df <- read_csv(file.path(dir_tab, paste0("bayesprism_fractions", run_tag, ".csv")), show_col_types = FALSE) %>%
  select(sample_id, freq = `CD4 CTL`) %>%
  inner_join(samples, by = "sample_id") %>%
  mutate(response = factor(response, levels = c("R", "NR")))
stopifnot(nrow(df) == nrow(samples), all(df$freq > 0))
write_csv(arrange(df, response, subject), file.path(dir_tab, paste0("fig_cd4ctl", run_tag, ".csv")))
print(wilcox.test(freq ~ response, data = df))

log_breaks <- as.vector(outer(c(1, 3), 10^(floor(log10(min(df$freq*100))):ceiling(log10(max(df$freq*100))))))

cd4ctl_plot <- function(summary_layer) {
  set.seed(1)
  df %>%
  ggplot(aes(x = response, y = freq*100)) +
    summary_layer +
    geom_beeswarm(size = 2.5, aes(color = response), cex = 2,
                  shape = 16, corral = "random") +
    scale_y_log10(breaks = log_breaks, labels = format(log_breaks, scientific = FALSE, drop0trailing = TRUE, trim = TRUE),
                  expand = expansion(mult = c(0.03, 0.12))) +
    scale_fill_manual(values = response_colors, name = NULL,
                      labels = paste0(names(response_colors), " (n = ", n[names(response_colors)], ")")) +
    scale_color_manual(values = response_colors, name = NULL,
                       labels = paste0(names(response_colors), " (n = ", n[names(response_colors)], ")")) +
    ylab("CD4 CTL (%)") +
    xlab("") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.location = "panel") +
    stat_compare_means(label = "p.format", p.digits = 2, size = 7,
                       label.x = 1.5, hjust = 0.5,
                       label.y = log10(max(df$freq*100) * 1.8)) +
    no_background +
    font_arial +
    ggtheme()
}

p <- cd4ctl_plot(stat_summary(fun = "median", fun.min = function(y) -Inf, fun.max = "median",
                         geom = "crossbar",
                         aes(fill = response), alpha = 0.3,
                         width = 0.75,
                         color = "grey30"))
save_svg(p, file.path(dir_fig, paste0("fig_cd4ctl_bar", run_tag, ".svg")), width = 14, height = 11)

p <- cd4ctl_plot(list(geom_violin(alpha = 0.3, aes(fill = response), scale = "width"),
                 stat_summary(fun = "median",
                              geom = "crossbar",
                              width = 0.75,
                              color = "grey30")))
save_svg(p, file.path(dir_fig, paste0("fig_cd4ctl_violin", run_tag, ".svg")), width = 14, height = 11)

p <- cd4ctl_plot(stat_summary(fun = "median",
                         geom = "crossbar",
                         width = 0.75,
                         color = "grey30"))
save_svg(p, file.path(dir_fig, paste0("fig_cd4ctl_dots", run_tag, ".svg")), width = 14, height = 11)
