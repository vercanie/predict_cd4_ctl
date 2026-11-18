# Paths, parameters and plot style shared by scripts 02-06

library(dplyr)
library(tidyr)
library(readr)
library(Matrix)
library(ggplot2)

# Paths ----------------------------------------------------------------------------------
proj <- "/path/to/project"
cr   <- file.path(proj, "clean_run")

raw <- list(
  risk_runs   = file.path(cr, "scripts/risk_runs.tsv"),
  risk_labels = file.path(proj, "data/raw/RISK/bulk_data.rd"),
  hgnc        = file.path(proj, "data/raw/annot/hgnc_complete_set.txt"),
  fig4_map    = file.path(cr, "scripts/fig4_rename_map.csv")
)
dir_atlas  <- file.path(cr, "output/data/atlas")
dir_salmon <- file.path(cr, "output/data/salmon")
dir_data  <- file.path(cr, "output/data")
dir_tab   <- file.path(cr, "output/tables")
dir_fig   <- file.path(cr, "output/figures")

# Parameters -----------------------------------------------------------------------------
min_cells_profile <- 10
min_profiles_type <- 20
sig_n_per_type    <- 60
sig_min_pct       <- 0.5
noise_genes <- "^(HB[ABDEGMQZ][0-9]?$|RP[LS][0-9]|MRP[LS]|MT-|MALAT1$|NEAT1$|XIST$|RPS4Y|DDX3Y|UTY$|KDM5D$|EIF1AY$)"
ncores <- as.integer(Sys.getenv("NCORES", "4"))
# sensitivity analysis: MIN_MAPPING=30
min_mapping <- as.numeric(Sys.getenv("MIN_MAPPING", "0"))
run_tag <- if (min_mapping > 0) paste0("_map", min_mapping) else ""

fig4_type <- function(fine, rename = read_csv(raw$fig4_map, show_col_types = FALSE)) {
  t <- coalesce(rename$fig4_name[match(fine, rename$object_annotations_fine)], fine)
  t <- sub("^Homaostatic", "Homeostatic", t)   # typo in the atlas
  if_else(t %in% c("Proliferating T cells", "Proliferating"), "Proliferating T", t)
}

# Plot style -----------------------------------------------------------------------------
ggtheme <- function() {
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    text = element_text(size = 20, colour = "black"),
    legend.text = element_text(size = 20),
    legend.key.size =  unit(10, units = "points")
  )
}

lineages <- c("B cells", "Endothelial cells", "Epithelial cells", "Mast cells", "Myeloid cells",
              "Plasma cells", "Proliferating cells", "Stromal cells", "T cells")
lineage_labels <- c("B cells", "Endothelial c.", "Epithelial cells", "Mast cells", "Myeloid cells",
                    "Plasma cells", "Proliferating c.", "Stromal cells", "T cells")
lineage_colors <- setNames(scales::hue_pal()(9), lineages)

systemfonts::register_font("Arial", plain = file.path(proj, "data/fonts/Arimo.ttf"))
font_arial <- theme(text = element_text(family = "Arial"), geom = element_geom(family = "Arial"))

save_svg <- function(plot, file, width, height) {
  ggsave(file, plot, width = width, height = height, units = "cm", bg = "transparent",
         fix_text_size = FALSE)
  svg <- readLines(file)
  writeLines(gsub('font-family: "Arimo"', 'font-family: "Arial"', svg, fixed = TRUE), file)
}

no_background <- theme(plot.background = element_blank(), panel.background = element_blank(),
                       legend.background = element_blank(), legend.key = element_blank())

response_colors <- c("R" = "green3", "NR" = "red")
response_strip_colors <- c("R" = "#B3F0B3", "NR" = "#FFB3B3")
