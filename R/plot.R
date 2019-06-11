plotCPN <- function(tes,
                    name = NULL)
{
  repname <- name
  cpn_df <- cpn(tes %>% filter(name == repname))
  
  cpn_df %>%
    ggplot(aes(pos_con, bps)) +
      geom_bar(stat = 'identity')
}                  

plotTECov <- function(counts, 
                      tes, 
                      pattern = NULL,
                      center = TRUE,
                      min_expr = 10)
{
  counts_te <- counts[which(type == 'te'),] 
  
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

plotGeneCov <- function(counts,
                        top_n = 100,
                        max_size = 1000)
{
  counts_gene <- counts[which(type == 'gene'), ]
  
  ggdata <- 
    counts_gene[, .(tot = sum(n)), by = c('name', 'pos_con')
      ][, gene_tot := sum(tot), by = 'name'
        ]
        
  top_genes <- 
    ggdata %>% 
    select(name, gene_tot) %>% 
    distinct %>% 
    top_n(top_n, gene_tot) %>% 
    pull(name)
  
  ggdata %>%
    filter(name %in% top_genes) %>%
    ggplot(aes(pos_con, name, fill = log10(tot))) +
      geom_raster() +
      #xlim(0, max_size) +
      scale_fill_distiller(palette = 'Spectral') +
      theme(axis.text        = element_blank(),
            axis.ticks       = element_blank(),
            axis.title       = element_blank(),
            panel.background = element_blank())
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

plotPercUnique <- function(counts)
{
  counts_te <- counts[which(type == 'te'), ]
  
  ggdata <- counts_te[, .(tot = sum(n_all), perc = sum(n) / sum(n_all) * 100), by = 'name']

  p <- 
    ggdata %>%
      mutate(label = ifelse(perc < 50, name, '')) %>%
    ggplot(aes(tot, perc, label = label)) + 
      geom_point() + 
      scale_x_log10() +
      ggrepel::geom_text_repel()

}