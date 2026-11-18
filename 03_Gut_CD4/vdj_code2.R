#! /usr/bin/Rscript
create.VDJ.metadata2 <- function(input.vdj){
  
  input.vdj <- input.vdj %>% dplyr::arrange(chain, cdr3, productive, reads, v_gene, j_gene, d_gene, c_gene) %>%
    dplyr::filter(cdr3 != 'None' & cdr3 != '') %>% 
    dplyr::filter(!duplicated(.[c('chain', 'barcode', 'cdr3')], fromLast = TRUE)) %>%
    dplyr::filter(barcode %in% md_sub$raw_barcode)
  
  # remove all non-productive VDJs
  input.vdj <- input.vdj %>% dplyr::filter(productive) %>% 
    dplyr::arrange(chain, barcode, desc(reads), desc(umis))
  
  # reformat data
  # first we'll treat columns we can't/don't wan't to just aggregate
  # Multi results will be aggregated
  vdj.data <- input.vdj %>% dplyr::select(barcode, chain, v_gene, d_gene, j_gene, cdr3, cdr3_nt)
  vdj.summary <- vdj.data %>% group_by(barcode, chain) %>% summarise(count = n())
  vdj.summary.A <- vdj.summary %>% dplyr::filter(chain == 'TRA') %>% dplyr::select(-chain) %>%
    dplyr::rename(count_TRA = count)
  vdj.summary.B <- vdj.summary %>% dplyr::filter(chain == 'TRB') %>% dplyr::select(-chain) %>% 
    dplyr::rename(count_TRB = count)
  vdj.summary.Multi <- vdj.summary %>% dplyr::filter(chain == 'Multi') %>% dplyr::select(-chain) %>%
    dplyr::rename(count_Multi = count)
  vdj.data.A1 <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRA') %>%  dplyr::filter(row_number()==1) %>%
    dplyr::select(-d_gene, -chain) %>%
    dplyr::rename(cdr3_A1 = cdr3, cdr3_A1_nt = cdr3_nt, 
                  v_gene_A1 = v_gene, j_gene_A1 = j_gene)
  vdj.data.A2 <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRA') %>% dplyr::filter(row_number()==2) %>%
    dplyr::select(-d_gene, -chain) %>% 
    dplyr::rename(cdr3_A2 = cdr3, cdr3_A2_nt = cdr3_nt, 
                  v_gene_A2= v_gene, j_gene_A2 = j_gene)
  vdj.data.B <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRB') %>% dplyr::filter(row_number()==1) %>%
    dplyr::select(-chain) %>% 
    dplyr::rename(cdr3_B = cdr3, cdr3_B_nt = cdr3_nt, v_gene_B= v_gene,
                  d_gene_B = d_gene, j_gene_B = j_gene)
  vdj.data.Multi <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'Multi') %>% 
    summarise(v_gene = paste(v_gene, collapse = ','), d_gene = paste(d_gene, collapse = ','), 
              j_gene = paste(j_gene, collapse = ','), cdr3 = paste(cdr3, collapse = ','),
              cdr3_nt = paste(cdr3_nt, collapse = ',')) %>%
    dplyr::rename(cdr3_Multi = cdr3, cdr3_Multi_nt = cdr3_nt, v_gene_Multi = v_gene,
                  d_gene_Multi = d_gene, j_gene_Multi = j_gene)
  
  # descriptors for each cell
  VDJ.cell.descr <- input.vdj %>% dplyr::select(barcode, is_cell) %>% 
    group_by(barcode) %>% dplyr::filter(row_number()==1)
  
  # descriptors for TCR that are not that instrumental
  # filter to keep only first two TRA and single TRB (arranged as previously)
  VDJ.TCR.descr <- input.vdj %>% 
    group_by(barcode, chain) %>%
    dplyr::filter((row_number() <= 2 & chain == 'TRA') | (row_number() <= 1 & chain == 'TRB') | 
                    (row_number() <= 1 & chain == 'Multi')) %>% ungroup %>%
    dplyr::select(barcode, chain, contig_id, high_confidence, productive, length, reads, umis) %>%
    group_by(barcode) %>% 
    summarise(chain = paste(chain, collapse = ','), contig_id = paste(contig_id, collapse = ','), 
              high_confidence = paste(high_confidence, collapse = ','), productive = paste(productive, collapse = ','),
              length = paste(length, collapse = ','), reads = paste(reads, collapse = ','),
              umis = paste(umis, collapse = ','))
  
  VDJ.meta.data <- VDJ.cell.descr %>% left_join(VDJ.TCR.descr, by = "barcode") %>% 
    ungroup %>% dplyr::arrange(barcode) %>% left_join(vdj.summary.A, by = "barcode") %>%
    left_join(vdj.summary.B, by = "barcode") %>%
    left_join(vdj.summary.Multi, by = "barcode") %>%
    left_join(vdj.data.A1, by = "barcode") %>% left_join(vdj.data.A2, by = "barcode") %>%
    left_join(vdj.data.B, by = "barcode") %>% left_join(vdj.data.Multi, by = "barcode")
  
  VDJ.meta.data <- as.data.frame(VDJ.meta.data)
  rownames(VDJ.meta.data) <- VDJ.meta.data$barcode
  VDJ.meta.data <- dplyr::select(VDJ.meta.data, is_cell, chain, contig_id, high_confidence, productive, length, 
                                 reads, umis, count_TRA, count_TRB, count_Multi, cdr3_A1, cdr3_A1_nt, 
                                 cdr3_A2, cdr3_A2_nt, cdr3_B, cdr3_B_nt, cdr3_Multi, cdr3_Multi_nt, v_gene_A1, j_gene_A1,
                                 v_gene_A2, j_gene_A2, v_gene_B, d_gene_B, j_gene_B, v_gene_Multi, d_gene_Multi, 
                                 j_gene_Multi)
  
  return(VDJ.meta.data)
}