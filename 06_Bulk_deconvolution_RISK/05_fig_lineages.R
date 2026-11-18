# Panel A: main population proportions per labelled patient

source("scripts/00_setup.R")
library(ggh4x)

samples <- read_csv(file.path(dir_tab, "risk_samples.csv"), show_col_types = FALSE) %>%
  filter(sample_id %in% read_csv(file.path(dir_tab, paste0("bayesprism_fractions", run_tag, ".csv")),
                                 col_select = sample_id, show_col_types = FALSE)$sample_id)
lineage <- read_csv(file.path(dir_tab, "cell_type_lineage.csv"), show_col_types = FALSE)
n <- table(samples$response)

df <- read_csv(file.path(dir_tab, paste0("bayesprism_fractions", run_tag, ".csv")), show_col_types = FALSE) %>%
  pivot_longer(-sample_id, names_to = "type", values_to = "freq") %>%
  inner_join(samples, by = "sample_id") %>%
  left_join(lineage, by = "type") %>%
  group_by(sample_id, subject, response, main) %>%
  summarise(freq = sum(freq), .groups = "drop") %>%
  mutate(response = factor(response, levels = c("R", "NR")),
         main = factor(main, levels = lineages),
         subject = factor(subject, levels = unique(subject[order(response, subject)])))
stopifnot(!anyNA(df$main), all(abs(tapply(df$freq, df$sample_id, sum) - 1) < 1e-6))
write_csv(arrange(df, response, subject, main), file.path(dir_tab, paste0("fig_lineages", run_tag, ".csv")))

gap <- unit(5.5, "pt")

p <- df %>%
ggplot(aes(x = subject, y = freq*100, fill = main)) +
  geom_bar(stat = "identity") +
  facet_grid2(cols = vars(response), scales = "free_x", space = "free_x",
              labeller = as_labeller(c(R  = paste0("R (n = ", n[["R"]], ")"),
                                       NR = paste0("NR (n = ", n[["NR"]], ")"))),
              strip = strip_themed(background_x = elem_list_rect(
                fill = response_strip_colors[levels(df$response)], colour = NA))) +
  scale_fill_manual(values = lineage_colors, breaks = lineages, labels = lineage_labels, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Estimated proportion (%)") +
  xlab("Participants") +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA),
        strip.text = element_text(size = 20),
        strip.placement = "outside",
        strip.switch.pad.grid = gap,
        panel.spacing.x = gap,
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black", linewidth = 0.8),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.8),
        axis.ticks.length.y = unit(6, "pt"),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.x = element_text(margin = margin(t = 9))) +
  no_background +
  font_arial +
  ggtheme()

save_svg(p, file.path(dir_fig, paste0("fig_lineages", run_tag, ".svg")), width = 36, height = 12)
