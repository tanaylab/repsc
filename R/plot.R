plotCellSize <- function(scSet)
{
  counts <- scSet@counts
  
  cell_size <- counts[, .(n_tot = sum(n)), by = c('barcode', 'type')]
  
  p <-
    cell_size %>% 
      ggplot(aes(n_tot, fill = type)) + 
        geom_histogram(bins = 75, alpha = 0.5, position = 'identity') + 
        scale_x_log10()
        
  return(p)
}

plotCPN <- function(tes,
                    name = NULL)
{
  repname <- name
  cpn_df <- cpn(tes %>% filter(name == repname))
  
  cpn_df %>%
    ggplot(aes(pos_con, bps)) +
      geom_bar(stat = 'identity')
}   

plotEnr <- function(scSet)
{
  counts <- scSet@counts
  
  ggdata <- 
    counts[, .(enr = sum(n)), by = c('barcode', 'peak', 'type')] %>% 
      spread(peak, enr, fill = 0, sep = '_') %>% 
      as_tibble %>%
      mutate(perc_peak = peak_TRUE / (peak_FALSE + peak_TRUE) * 100)
  
  p <-
    ggdata %>%
      ggplot(aes(perc_peak, fill = type)) +
        geom_histogram(bins = 75, position = 'identity', alpha = 0.25)
  
  return(p)
}               

plotTECov <- function(scSet, 
                      pattern = NULL,
                      center = TRUE,
                      min_expr = 10)
{
  counts    <- scSet@coutns
  counts_te <- counts[which(type == 'te'),] 
  tes       <- scSet@tes
  
  if (!is.null(pattern))
  {
    counts_te <- counts_te[grep(pattern, name), ]
  }
  
  # get TE copynumber
  cpn_df <- cpn(tes)
  
  # aggregate counts across cells and TE loci
  counts_aggr <- aggr(counts_te)
  
  # calc density
  con_signal <- 
    merge(cpn_df, counts_aggr, all.x = TRUE, by = c('name', 'pos_con'))[which(is.na(n_tot)), n_tot := 0
        ][, dens := n_tot / bps
          ][which(is.na(dens)), dens := 0]
  
  # exclude NA pos_con
  con_signal <- con_signal[!which(is.na(pos_con)),]
  
  # centering and transformations
  con_signal[, maxi := dens == max(dens), by = 'name']
  con_signal[, dens_scaled := dens / dens[maxi][1], by = 'name']
  con_signal[, dens_log := cut(dens, breaks = c(0, 0.001 * 2^(1:16)), include.lowest = TRUE) , by = 'name']
  con_signal[, umis_log := cut(n_tot, breaks = c(0, 2^(1:20)), include.lowest = TRUE, labels = FALSE) , by = 'name']
  if (center)
  {
    con_signal[, pos_con_cent := as.numeric(pos_con) - as.numeric(pos_con[maxi][1]), by = 'name']
  }
  
  # add repclass annotation
  con_signal <- merge(con_signal, unique(as.data.table(tes)[, c('name', 'repclass')]), all.x = TRUE, by = 'name')
  
  # plot heatmap
  con_signal %>%
    ggplot(aes(pos_con_cent, name, fill = dens_log)) +
      geom_raster() +
      xlim(-50, 50) + 
      scale_fill_manual(values = colorRampPalette(c(rev(RColorBrewer::brewer.pal(11, 'Spectral'))))(length(levels(con_signal$dens_log)) + 1)) + 
      facet_wrap(~repclass, scales = 'free_y') +
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            panel.grid.major = element_blank())
}

plotGeneCov <- function(scSet,
                        highest_n = 2500)
{
  protocol <- scSet@protocol
  counts   <- scSet@counts
  
  counts_gene <- counts[which(type == 'gene'), ]
  
  ggdata <- 
    counts_gene[, .(tot = sum(n)), by = c('name', 'pos_con')
      ][, gene_tot := sum(tot), by = 'name'
        ]
        
  top_genes <- 
    ggdata %>% 
    select(name, gene_tot) %>% 
    distinct %>% 
    top_n(highest_n, gene_tot) %>% 
    pull(name)
  
  
  if (protocol == 'fiveprime' | protocol == 'threeprime')
  {
    # scale UMIs per gene
    ggdata[, maxi := tot == max(tot), by = 'name']
    ggdata[, tot_scaled := tot / tot[maxi][1], by = 'name']
    
    # center on max bin
    ggdata[, pos_con_cent := as.numeric(pos_con) - as.numeric(pos_con[maxi][1]), by = 'name']
  }
  
  p <- 
    ggdata %>%
      filter(name %in% top_genes) %>%
      ggplot(aes(pos_con_cent, name, fill = log10(tot))) +
        geom_raster() +
        xlim(-100, 100) +
        scale_fill_distiller(palette = 'Spectral') +
        theme(axis.text        = element_blank(),
              axis.ticks       = element_blank(),
              axis.title       = element_blank(),
              panel.background = element_blank())
}

plotMitoRibo <- function(scSet)
{
  ribo_prefixes <- c('Human' = '^RP[LS]', 'Mouse' = '^Rp[ls]') 
  genes         <- scSet@genes
  organism      <- scSet@genome@common_name
  counts_gene   <- scSet@counts[which(type == 'gene'), ]
  ribo_prefix   <- ribo_prefixes[organism]
  
  mt_genes <- 
    genes %>%
      filter(gene_type %in% c('Mt_tRNA', 'Mt_rRNA')) %>%
      as_tibble %>%
      pull(gene_name) %>%
      unique
      
  ribo_genes <-
    genes %>%
      filter(grepl(ribo_prefix, gene_name)) %>%
      as_tibble %>%
      pull(gene_name) %>%
      unique
   
  # tag mito/ribo genes
  counts_gene <- 
    counts_gene[, ':=' (mito = name %in% mt_genes, ribo = name %in% ribo_genes)]
  
  # summarize percentage  
  counts_aggr <-
    counts_gene[, .(n_tot = sum(n), n_ribo = sum(n[ribo]), n_mito = sum(n[mito])), by = c('barcode')
      ][, .(perc_ribo = n_ribo / n_tot * 100, perc_mito = n_mito / n_tot * 100)]
  
  p <- 
    counts_aggr %>% 
      ggplot(aes(perc_ribo, perc_mito)) + 
        geom_point(size = 0.25) + 
        scale_y_log10() + 
        scale_x_log10()
  
  return(p)
}

plotMSA <- function(msa, 
                    cluster = TRUE,
                    type = 'base',
                    max_gap = 0.95,
                    threads = 1,
                    ds = NULL)
{
  if (!is.null(ds))
  {
    msa <- sample(msa, ds)
  }
  
  msa <- DNAMultipleAlignment(msa)
  
  # remove gapped columns
  gap_cols <- which(consensusMatrix(msa, as.prob = T)['-',] > max_gap)
  colmask(msa) <- IRanges(start = gap_cols, end = gap_cols)
  
  msa <- DNAStringSet(msa)
  
  # Remove short loci after gap removal
  msa <- msa[width(DECIPHER::RemoveGaps(msa)) >= 6]
  
  if (cluster)
  {
    message ('Clustering loci')
    # msa <- DNAStringSet(msa)
    # d <- DECIPHER::DistanceMatrix(msa, processors = 70)
    # clusts <- IdClusters(myDistMatrix = d, method = "UPGMA", processors = 70)
    # msa <- msa[order(clusts[,1])]
    
    d <- tgs_dist(oligonucleotideFrequency(msa, width = 6))
    clusts <- hclust(d, 'ward.D2')

    msa <- msa[clusts$order]
  }
    
  if (type == 'base')
  {
    #msa_bin <- as.DNAbin(msa)
    
    ape::image.DNAbin(as.DNAbin(msa), col = c("#5DA731", "#D63317", "#DFD93C", "#2F4C9B", "lightgrey", "white"), show.labels = FALSE)
  } else {
    msa_mat <- as.matrix(msa)
    msa_dt <- as.data.table(msa_mat)
    msa_dt$id_unique <- rownames(msa)
    msa_long <- melt(msa_dt, id.vars = 'id_unique')
    
    # plot ggplot
    msa_long %>% 
      ggplot(aes(variable, id_unique, fill = value)) + 
        geom_raster() +
        scale_fill_manual(values = c("white", "#5DA731", "#D63317", "#DFD93C", "#2F4C9B")) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              panel.background = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
  }
}

plotPercUnique <- function(scSet,
                           type = c('cells', 'family'))
{
  counts <- scSet@counts
  
  if (type == 'cells')
  {
    ggdata <- counts[, .(n = sum(n), n_all = sum(n_all)), by = 'barcode'][, perc_unique := n / n_all * 100]
    p <-
      ggdata %>%
        ggplot(aes(perc_unique)) + 
          geom_histogram(bins = 75)
  }
  
  if (type == 'family')
  {
    counts_te <- counts[which(type == 'te'), ]
    ggdata <- counts_te[, .(tot = sum(n_all), perc = sum(n) / sum(n_all) * 100), by = 'name']
   
    p <- 
      ggdata %>%
        mutate(label = ifelse(perc < 50, name, '')) %>%
      ggplot(aes(tot, perc, label = label)) + 
        geom_point(size = 0.25) + 
        scale_x_log10() #+
        #ggrepel::geom_text_repel()
  }

  return(p)
}

plotTELoad <- function(scSet)
{
  counts <- scSet@counts
  
  te_load <- compTELoad(counts)
  
  p <-
    te_load %>%
      ggplot(aes(te_load)) +
        geom_histogram(bins = 75) +
        scale_x_log10()
        
  return(p)
}

plotQC <- function(counts,
                   genes,
                   tes)
{
  counts_gene <- counts[which(type == 'gene'), ]
  counts_gene$repclass <- NA
  counts_te   <- counts[which(type == 'te'), ]
  counts_te   <- merge(counts_te, unique(as.data.table(tes)[, c('name', 'repclass')]), by = 'name')
  counts      <- rbindlist(list(counts_gene, counts_te))
  
  
  p_size          <- plotCellSize(scSet)
  p_load          <- plotTELoad(scSet)
  p_mitoribo      <- plotMitoRibo(scSet)
  p_unique_cells  <- plotPercUnique(scSet, type = 'cells')
  p_unique_fams   <- plotPercUnique(scSet, type = 'family')
  p_gene_cov      <- plotGeneCov(scSet)
  p_te_cov        <- plotTECov(counts)
  
  cowplot::plot_grid(p_size, 
                     p_load, 
                     p_unique_cells, 
                     p_gene_cov,
                     p_unique_fams)

}