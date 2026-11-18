# RISK bulk: tximport gene counts, symbol matching to the atlas, response labels, salmon QC

source("scripts/00_setup.R")

runs <- read_tsv(raw$risk_runs, show_col_types = FALSE)
stopifnot(nrow(runs) == 199, !anyDuplicated(runs$sample_id))

# Response labels -------------------------------------------------------------------------
load(raw$risk_labels)
samples <- bulk$design %>%
  tibble::rownames_to_column("sample_id") %>%
  filter(dataset == "RISK", !is.na(response), sample_id %in% runs$sample_id) %>%
  transmute(sample_id, subject = as.character(subject), sex = as.character(sex),
            status = as.character(status), response = as.character(response))
stopifnot(!anyDuplicated(samples$subject), all(samples$response %in% c("R", "NR")))
samples %>% count(response, status) %>% print()
write_csv(samples, file.path(dir_tab, "risk_samples.csv"))

# Salmon -> gene counts (tximport) --------------------------------------------------------
files <- setNames(file.path(dir_salmon, runs$sample_id, "quant.sf"), runs$sample_id)
stopifnot(all(file.exists(files)))
tx2gene <- read_tsv(file.path(dir_data, "salmon_index/tx2gene.tsv"), col_names = c("tx", "gene_id", "gene"),
                    show_col_types = FALSE) %>% select(tx, gene)
txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
counts <- txi$counts

qc <- bind_rows(lapply(runs$sample_id, function(s) {
  m <- jsonlite::read_json(file.path(dir_salmon, s, "aux_info/meta_info.json"))
  tibble(sample_id = s, reads = m$num_processed, percent_mapped = m$percent_mapped)
}))
write_csv(qc, file.path(dir_tab, "risk_salmon_qc.csv"))
message("salmon mapping rate (%): ", paste(round(quantile(qc$percent_mapped, c(0, .5, 1)), 1), collapse = " / "),
        " (min / median / max)")
message("reads assigned per sample (M): ",
        paste(round(quantile(colSums(counts), c(0, .5, 1)) / 1e6, 1), collapse = " / "))

# Symbol matching to the atlas ------------------------------------------------------------
atlas_genes <- readLines(file.path(dir_atlas, "genes.tsv"))
hgnc <- read_tsv(raw$hgnc, col_select = c(symbol, prev_symbol, alias_symbol),
                 show_col_types = FALSE)
synonyms <- bind_rows(
  hgnc %>% separate_longer_delim(prev_symbol, "|") %>% select(a = symbol, b = prev_symbol),
  hgnc %>% separate_longer_delim(alias_symbol, "|") %>% select(a = symbol, b = alias_symbol)) %>%
  filter(!is.na(b), b != "")
synonyms <- bind_rows(synonyms, rename(synonyms, a = b, b = a)) %>% distinct()
symbol_map <- synonyms %>%
  transmute(old = a, new = b) %>%
  filter(old %in% setdiff(rownames(counts), atlas_genes),
         new %in% atlas_genes, !new %in% rownames(counts)) %>%
  group_by(old) %>% filter(n() == 1) %>% ungroup() %>%   # one-to-one matches only
  group_by(new) %>% filter(n() == 1) %>% ungroup()
message("RISK symbols renamed to atlas symbols: ", nrow(symbol_map))
rownames(counts) <- coalesce(symbol_map$new[match(rownames(counts), symbol_map$old)], rownames(counts))
stopifnot(!anyDuplicated(rownames(counts)))

message("RISK counts: ", nrow(counts), " genes x ", ncol(counts), " samples; ",
        sum(rownames(counts) %in% atlas_genes), " genes shared with the atlas")
saveRDS(counts, file.path(dir_data, "risk_counts.rds"))
