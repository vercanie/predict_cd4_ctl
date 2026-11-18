# BayesPrism reference: mean counts per cell per (cell type, donor) profile, marker genes

source("scripts/00_setup.R")

risk_genes <- rownames(readRDS(file.path(dir_data, "risk_counts.rds")))
genes <- readLines(file.path(dir_atlas, "genes.tsv"))
profiles <- read_csv(file.path(dir_atlas, "profiles.csv"), show_col_types = FALSE) %>%
  mutate(type = fig4_type(annotations_fine), key = paste(type, donor, sep = "---"))
counts <- readMM(file.path(dir_atlas, "profiles_counts.mtx"))
dimnames(counts) <- list(genes, profiles$profile)

# Main population per cell type -------------------------------------------------------------
lineage <- profiles %>%
  group_by(type, main) %>% summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  group_by(type) %>% slice_max(n_cells, n = 1, with_ties = FALSE) %>% ungroup() %>%
  # "Proliferating T" merges two fine labels filed under different main populations
  mutate(main = if_else(type == "Proliferating T", "Proliferating cells", main)) %>%
  select(type, main)
stopifnot(all(lineage$main %in% lineages))
write_csv(lineage, file.path(dir_tab, "cell_type_lineage.csv"))

# Sum profiles per (Fig. 4 type, donor) -----------------------------------------------------
keys <- unique(profiles$key)
G <- sparseMatrix(i = match(profiles$key, keys), j = seq_len(nrow(profiles)), x = 1,
                  dims = c(length(keys), nrow(profiles)))
counts <- as(counts, "CsparseMatrix") %*% t(G)
colnames(counts) <- keys
ref_profiles <- profiles %>%
  group_by(key, type, donor) %>% summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  slice(match(keys, key))

# Filter profiles, cell types and genes -----------------------------------------------------
keep_types <- ref_profiles %>% filter(n_cells >= min_cells_profile) %>%
  count(type) %>% filter(n >= min_profiles_type) %>% pull(type)
keep <- ref_profiles$n_cells >= min_cells_profile & ref_profiles$type %in% keep_types
ref_profiles <- ref_profiles[keep, ]
counts <- counts[!grepl(noise_genes, rownames(counts)), keep]
counts <- counts[intersect(rownames(counts), risk_genes), ]
message("reference: ", length(keep_types), " cell types, ", ncol(counts), " profiles, ",
        n_distinct(ref_profiles$donor), " donors, ", nrow(counts), " genes shared with RISK")

# Mean counts per cell ----------------------------------------------------------------------
ref <- counts %*% Diagonal(x = 1 / ref_profiles$n_cells)
colnames(ref) <- colnames(counts)

# Marker genes ------------------------------------------------------------------------------
cp10k <- t(t(ref) / colSums(ref)) * 1e4
types <- sort(unique(ref_profiles$type))
mean_type <- sapply(types, function(t) rowMeans(cp10k[, ref_profiles$type == t, drop = FALSE]))
pct_type  <- sapply(types, function(t) rowMeans(ref[, ref_profiles$type == t, drop = FALSE] > 0))
markers <- unique(unlist(lapply(types, function(t) {
  other <- apply(mean_type[, types != t, drop = FALSE], 1, max)
  lfc <- log2((mean_type[, t] + 0.01) / (other + 0.01))
  ok <- pct_type[, t] >= sig_min_pct & lfc > 1
  head(names(sort(lfc[ok], decreasing = TRUE)), sig_n_per_type)
})))
message("marker genes: ", length(markers))

saveRDS(list(reference = t(as.matrix(ref[markers, ])), cell_type = ref_profiles$type),
        file.path(dir_data, "reference.rds"))
