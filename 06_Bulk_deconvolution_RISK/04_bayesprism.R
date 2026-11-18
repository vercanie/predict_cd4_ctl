# BayesPrism, all biopsies jointly; final cell-type theta

source("scripts/00_setup.R")
library(BayesPrism)

bulk <- readRDS(file.path(dir_data, "risk_counts.rds"))
qc   <- read_csv(file.path(dir_tab, "risk_salmon_qc.csv"), show_col_types = FALSE)
bulk <- bulk[, qc$sample_id[qc$percent_mapped >= min_mapping]]
message("biopsies deconvolved: ", ncol(bulk), " (salmon mapping rate >= ", min_mapping, "%)")
ref  <- readRDS(file.path(dir_data, "reference.rds"))
genes <- colnames(ref$reference)

prism <- new.prism(reference = ref$reference, mixture = t(bulk[genes, ]),
                   input.type = "count.matrix",
                   cell.type.labels = ref$cell_type, cell.state.labels = ref$cell_type,
                   key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
res <- run.prism(prism, n.cores = ncores)

theta <- get.fraction(res, which.theta = "final", state.or.type = "type")
stopifnot(nrow(theta) == ncol(bulk), all(abs(rowSums(theta) - 1) < 1e-6))
theta %>% as.data.frame() %>% tibble::rownames_to_column("sample_id") %>%
  write_csv(file.path(dir_tab, paste0("bayesprism_fractions", run_tag, ".csv")))
