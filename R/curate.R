# decide how to deal with antisense counts!
resolveAmbiguity <- function(scSet,
                             n_cells = 5000,
                             min_cor = 0.2)
{
  message ('Resolving TE-Gene ambiguity')
  set.seed(scSet@seed)
  tes         <- tes(scSet)
  genes       <- genes(scSet)
  counts      <- counts(scSet)
  counts_gene <- counts[which(type == 'gene'), ]
  counts_te   <- counts[which(type == 'te'), ]
  
  # add exonic info to counts
  counts_te$exonic <- counts_te$id_unique %in% tes[tes$exonic]$id_unique
  
  # summarize non-exonic counts per familiy and cell
  fam_sc_counts <- counts_te[which(!exonic), .(n = sum(n)), by = c('name', 'barcode')]
  
  # summarize exonic counts per locus and cell
  loc_sc_counts <- counts_te[which(exonic), .(n = sum(n), name = name[1]), by = c('id_unique', 'barcode')]

  # create mat of fam and loci counts
  smat <- Reputils::longToSparse(rbindlist(list(fam_sc_counts, loc_sc_counts[, 1:3]), use.names = FALSE))

  # sample cells if too many
  if (length(scSet@cells) > n_cells)
  {
    smat <- smat[, sample(colnames(smat), n_cells)]
  }
  
  # all vs all correlation
  mat_cor <- tgstat::tgs_cor(t(as.matrix(smat)), spearman = TRUE)
  
  # loci vs fam correlation
  loci    <- grep('|', rownames(mat_cor), fixed = TRUE)
  mat_cor <- mat_cor[loci, -loci]
  
  # get max correlated fam per locus
  max_cor_df <- 
    as.data.table(melt(mat_cor, varnames = c('id_unique', 'name')))[, .(fam_max = name[value == max(value, na.rm = TRUE)], cor = max(value, na.rm = TRUE)), by = 'id_unique']
  max_cor_df <- merge(max_cor_df, unique(loc_sc_counts[, c('id_unique', 'name')])) # add repname to locus
  
  # label and return loci correlated to their fam
  max_cor_df[, cor2fam := name == fam_max & cor >= min_cor]
  loci_rescued <- max_cor_df[which(cor2fam), ]$id_unique
   
  # kill gene counts likely coming from TE
  genes2throw   <- sort(unique(unlist(strsplit(sapply(strsplit(loci_rescued, '|', fixed = TRUE), function(x) x[3]), '_', fixed = TRUE))))
  counts_gene_f <- counts_gene[!which(name %in% genes2throw), ]
  
  # kill non-correlated exonic TE loci
  counts_te_f   <- counts_te[which(!exonic | id_unique %in% loci_rescued), !'exonic']
 
  #message ('Rescued ', length(loci_rescued), ' TE loci in favour of ', paste0(genes2throw, collapse = ' '))
  message ('Rescued ', paste0(loci_rescued, collapse = ' '), ' TE loci')
 
  # plot locus correlation
  p_loc_fam_cor <-
    max_cor_df %>% 
    ggplot(aes(cor, fill = cor2fam)) + 
      geom_density(alpha = 0.25)
 
  # update scSet
  scSet@counts            <- rbindlist(list(counts_gene_f, counts_te_f))
  scSet@plots$loc_fam_cor <- p_loc_fam_cor
  
  print(p_loc_fam_cor)
  return(scSet)
}