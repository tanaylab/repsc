plotCellSize <- function(scSet,
                         what = c('hist', 'rank'))
{
  cstats <- scSet@cstats
  
  if (sum(dim(cstats)) == 0)
  {
    cstats <- cstats(scSet)    
  }
 
  if (what == 'hist')
  {  
    p <-
      cstats %>% 
        ggplot(aes(x = size_genic, y = meta, fill = meta)) + 
          ggridges::geom_density_ridges(alpha = 0.95) +
          #geom_histogram(bins = 75, alpha = 0.25) + 
          #scale_x_log10() + 
          #facet_wrap(~pass_size, scales = 'free_y') +
          scale_fill_viridis_d() + 
          theme(legend.position = 'none')
  }
  
  if (what == 'rank')
  {
    p <-
      cstats[which(size_genic > 10), ] %>%
        ggplot(aes(rank, size_genic, col = meta)) + 
          geom_point(size = 0.1) + 
          scale_x_log10() + 
          scale_y_log10() +
          #geom_hline(yintercept = size_t) +
          #scale_size_manual(values = c(8, 0.01)) + 
          scale_colour_viridis_d() + 
          theme(legend.position = 'none') 
  }
  
  return(p)
}

plotCorMat <- function(scSet, 
                       ds = NULL,
                       what = 'features',
                       nclust = 8,
                       min_cor = 0.1,
                       min_expr = 100)
{
  set.seed(scSet@seed)
  scSet     <- suppressWarnings(annoClass(scSet))
  gstats    <- scSet@gstats
  counts    <- Repsc::counts(scSet)
  counts_ds <- scSet@counts_ds
  feats     <- gstats[which(type == 'te' & feature), ]$id_unique 
    
  if (what == 'features')
  {
    # aggregate feature counts per barcode and id_unique
    counts_ds_f <- counts_ds[which(id_unique %in% feats), .(n_ds = sum(n_ds)), by = c('id_unique', 'barcode')]
  
    smat    <- Reputils::longToSparse(counts_ds_f)
    cor_mat <- tgs_cor(log2(1 + 7 * as.matrix(t(smat))))
    
    # cluster based on cor mat
    #ki <- sort(kmeans( cor_mat, centers = 8)$cluster)
    #kj <- sort(kmeans(t(cor_mat_fvsg_f), centers = 8)$cluster)
    # dist_mat <- tgs_dist(cor_mat)
    # ki       <- sort(cutree(hclust(dist_mat, 'ward.D'), k = nclust))
    
    cors_df  <- as.data.table(melt(cor_mat, varnames = c('id_unique', 'id_unique2'), value.name = 'cor')) 
   
    # add cluster to long data.frame and reorder factor
    cors_df <- 
      merge(cors_df,
            gstats[, .(id_unique = id_unique, row_clust = module)], all.x = TRUE, by = 'id_unique')[, id_unique := forcats::fct_reorder(id_unique, row_clust)]
    cors_df <- 
      merge(cors_df,
            gstats[, .(id_unique2 = id_unique, col_clust = module)], all.x = TRUE, by = 'id_unique2')[, id_unique2 := forcats::fct_reorder(id_unique2, col_clust)]
  
    # add family annotation columns
    #cors_df <- merge(cors_df, unique(counts[, c('id_unique', 'name')]), all.x = TRUE, by = 'id_unique')
    
    # plot
    steps <- as.numeric(cumsum(table(gstats[which(type == 'te'),]$module)))
    
    # cap cor valid
    cap <- 0.4
    cors_df[which(cor > 0.4), cor := cap]
    cors_df[which(cor < -0.4), cor := -cap]
    
    p_cor_hmap <-
      cors_df %>%
          ggplot(aes(id_unique, id_unique2, fill = cor)) + 
            geom_raster() +
            scale_fill_gradientn(colors = colorspace::diverging_hcl(palette = 'Blue-Red 3', n = 256), limits = c(-cap, cap)) +
            geom_hline(yintercept = steps, lwd = 0.25) +
            geom_vline(xintercept = steps, lwd = 0.25) + 
            xlab('TE features') +
            ylab('TE features') +
            theme(legend.position = 'none',
                  legend.key.width = unit(.3,"cm"),
                  legend.key.height = unit(.15,"cm"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
  }
  
  if (what == 'family')
  {
  # aggregate per family and gene
  counts_ds <- counts_ds[, .(n_ds = sum(n_ds)), by = c('name', 'barcode', 'type', 'class')]
  
  # add class counts
  counts_ds <-
    rbindlist(list(counts_ds, 
                   counts_ds[, .(n_ds = sum(n_ds), name = class), by = c('barcode', 'class')
                    ][which(class %in% c('DNA', 'SINE', 'LINE', 'LTR')), ]), 
              use.names = TRUE, 
              fill = TRUE)  
  
  # filter for feature fams/genes only
  counts_ds_f <- 
    counts_ds[which(name %in% c('DNA', 'SINE', 'LTR', 'LINE', unique(gstats[which(feature & class %in% c('DNA', 'LINE', 'SINE', 'LTR', 'protein_coding')), ]$name))), ]
  
  # vectors of TE families and protein coding genes for later subsetting
  fams  <- c(unique(counts_ds_f[which(type == 'te'), ]$name), 'DNA', 'SINE', 'LTR', 'LINE')
  genes <- unique(counts_ds_f[which(type == 'gene'), ]$name)
  
  # create sparse matrix, filter for min expression and log transform
  smat      <- Reputils::longToSparse(counts_ds_f[, c('name', 'barcode', 'n_ds')])
  smat_f    <- smat[Matrix::rowSums(smat) >= min_expr, ]
  mat_f_trans <- as.matrix(log2(1 + 7 * smat_f))
  
  # calculate pearson correlation matrix
  cor_mat <- tgs_cor(t(mat_f_trans))
  
  # subset for rows to be families and columns to be genes
  cor_mat_fvsg <- cor_mat[rownames(cor_mat) %in% fams, colnames(cor_mat) %in% genes]
  
  # subset genes based on min cor threshold
  above_t_mat    <- cor_mat_fvsg > min_cor
  cor_mat_fvsg_f <- cor_mat_fvsg[rowSums(above_t_mat) > 0, colSums(above_t_mat) > 0]
  
  # cluster based on cor mat
  #ki <- sort(kmeans(cor_mat_fvsg_f, centers = 8)$cluster)
  #kj <- sort(kmeans(t(cor_mat_fvsg_f), centers = 8)$cluster)
  
  ki <- hclust(tgs_dist(cor_mat_fvsg_f), 'ward.D2')
  kj <- hclust(tgs_dist(t(cor_mat_fvsg_f)), 'ward.D2')
  
  # tidy fam gene correlations
  fg_cors_df <- as.data.table(melt(cor_mat_fvsg_f, varnames = c('fam', 'gene'), value.name = 'cor')) 
  
  # add cluster to long data.frame and reorder factor
  fg_cors_df[, ':=' (fam = factor(fam, levels = ki$labels[ki$order], ordered = TRUE), gene = factor(gene, levels = kj$labels[kj$order], ordered = TRUE))]
  # fg_cors_df <- 
    # merge(fg_cors_df,
        # data.table(fam = names(ki), fam_clust = ki), all.x = TRUE, by = 'fam')[, fam := forcats::fct_reorder(fam, fam_clust)]
        
  # label top correlated te or gene per cluster
  labels_gene <- levels(fg_cors_df$gene)
  labels_gene <- labels_gene[nchar(labels_gene) < 8] 
  labels_gene <- labels_gene[seq(1, length(labels_gene), l = min(length(labels_gene), 50))]
  labels_fam <- levels(fg_cors_df$fam)
  labels_fam <- unique(c('DNA', 'SINE', 'LINE', 'LTR', labels_fam[seq(1, length(labels_fam), l = 20)]))
  
  # cap cor
  cap <- 0.4
  fg_cors_df[which(cor > 0.4),  cor := cap]
  fg_cors_df[which(cor < -0.4), cor := -cap]
  
  # plot
  p_cor_hmap <-
    fg_cors_df %>%
        ggplot(aes(gene, fam, fill = cor)) + 
          geom_raster() +
          scale_x_discrete(breaks = as.character(labels_gene)) + 
          scale_y_discrete(breaks = as.character(labels_fam)) +
          scale_fill_gradientn(colors = colorspace::diverging_hcl(palette = 'Blue-Red 3', n = 256), limits = c(-cap, cap)) +
          theme(legend.position = 'none') +
          xlab('') +
          ylab('')
  }
  return(p_cor_hmap)
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

plotEnr <- function(scSet,
                    type = c('gene', 'te'),
                    what = c('transcript', 'cells', 'module'))
{
  what_type <- type
  counts <- suppressWarnings(counts(annoClass(scSet)))[which(type == what_type), ]
  genes  <- genes(scSet)
  
  # Enrichment by transcript length
  if (what == 'transcript')
  {
    gene_size_df <- as.data.table(genes)[, .(size = max(position_start, position_end), n_junc = .N), by = 'name']
    gene_enr_df  <- 
      counts[which(type == 'gene'), 
        ][, .(n_tot = sum(n_all), class = class[1]), by = c('name', 'peak')
          ][which(n_tot > 0), # genes with 0 n but >0 n_all are thrown here
            ] %>%
            tidyr::spread(peak, n_tot, fill = 0, sep = '_') %>%
            as.data.table() %>%
            .[, ':=' (n_tot = peak_FALSE + peak_TRUE, enr = peak_TRUE / (peak_FALSE + peak_TRUE))]
    
    # add gene size and n_junc
    gene_stats <- 
      merge(gene_enr_df, gene_size_df, all.x = TRUE, by = 'name')[, size_bin := cut(size, 
                                                                                breaks = c(0, 1000, 5000, 1e4, 1e6), 
                                                                                include.lowest = TRUE)]
        #][, ':=' (size_median = paste0(range(size), collapse = '-'), enr_median = median(enr)), by = 'size_bin']
    
    # p <-
      # gene_stats %>% 
        # filter(n_tot > 10) %>% 
        # ggplot(aes(x = size_bin, y = enr, fill = class)) + 
          #geom_violin() + 
          # geom_boxplot(outlier.shape = NA) +
          # facet_wrap(~class, nrow = 1) +
          # theme(legend.position = 'none')
    
    p <- gene_stats %>%
      ggplot(aes(peak_TRUE, n_tot, col = class)) +
        geom_point(size = 0.1, alpha = 0.25) +
        scale_y_log10() +
        scale_x_log10() +
        facet_wrap(~size_bin, nrow = 1) +
        theme(legend.position = 'none')
                 
    # p2 <- gene_stats %>% ggplot(aes(size_bin, enr)) + geom_boxplot(outlier.shape = NA)
    # p  <- cowplot::plot_grid(p1, p2, rel_widths = c(4, 1))
  }

  if (what == 'cells')
  {
    ggdata <- 
      counts[, .(n_tot = sum(n_all)), by = c('barcode', 'peak', 'class')
        ][, .(enr = n_tot[peak] / sum(n_tot)), by = c('barcode', 'class')]
  
    p <-
      ggdata %>%
        ggplot(aes(enr, fill = class)) +
          geom_density(alpha = 0.25) +
          #geom_histogram(bins = 75, position = 'identity', alpha = 0.25) +
          theme(legend.position = 'top')
  }
  
  if (what == 'module')
  {
    gstats <- scSet@gstats[which(type == 'te' & feature), ]
    
    # Fishers
    pop_size   <- nrow(gstats)
    module_size <- gstats %>% count(module) %>% dplyr::rename(module_size = n)
    fams_in_pop <- gstats %>% count(name) %>% dplyr::rename(inpop_size = n)
    fams_in_module <- gstats %>% count(name, module) %>% dplyr::rename(inmodule = n)

    dt <-
      left_join(fams_in_module, fams_in_pop) %>%
      left_join(., module_size) %>%
      mutate(outmodule = inpop_size - inmodule) %>% 
      mutate(others_inmodule = module_size - inmodule,
             others_outmodule = pop_size - inpop_size) %>%
      as.data.table()
  
    ggdata <- 
      dt[, pval := fisher.test(matrix(c(inmodule, outmodule, others_inmodule, others_outmodule), ncol = 2), alternative = 'greater')$p.value, by = c('name', 'module')
        ][, qval := p.adjust(pval, 'fdr')]  
    
    # plot
    p <-
      ggdata %>%
      filter(qval < 0.05) %>%
      #mutate(qval = cut(qval, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), include.lowest = TRUE)) %>%
      ggplot(aes(as.character(module), name, fill = -log10(qval))) + 
        geom_raster() +
        geom_text(aes(label = inmodule), size = 2.5) + 
        xlab('') + 
        scale_fill_distiller(palette = 'BuPu', direction = 1) +
        theme(legend.position = 'none',
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(1,"cm"))
  }
  
  return(p)
}               

plotTECov <- function(scSet, 
                      pattern = NULL,
                      center = TRUE,
                      min_expr = 100)
{
  protocol  <- scSet@protocol
  counts    <- Repsc::counts(scSet)
  counts_te <- counts[which(type == 'te'),] 
  cpn_df    <- scSet@cpn[which(type == 'te'), ]
  
  # aggregate counts across cells and TE loci
  counts_aggr <- aggr(counts_te, min_expr = min_expr)
  
  # only retain expressed fams
  cpn_df <- cpn_df[which(name %in% unique(counts_aggr$name)), ]
  
  # calc density
  con_signal <- 
    merge(cpn_df, counts_aggr, all.x = TRUE, by = c('name', 'pos_con', 'type'))[which(is.na(n_tot)), n_tot := 0
        ][, dens := n_tot / bps
          ][which(is.na(dens)), dens := 0]
  
  # exclude NA pos_con
  con_signal <- con_signal[!which(is.na(pos_con)),]
  
  # centering and transformations
  con_signal[, ':=' (max_dens = dens == max(dens), max_umis = n_tot == max(n_tot)), by = 'name']
  con_signal[, dens_scaled := dens / dens[max_dens][1], by = 'name']
  con_signal[, umis_scaled := n_tot / n_tot[max_umis][1], by = 'name']
  #con_signal[, dens_log := cut(dens, breaks = c(0, 0.001 * 2^(1:16)), include.lowest = TRUE) , by = 'name']
  #con_signal[, umis_log := cut(n_tot, breaks = c(0, 2^(1:20)), include.lowest = TRUE) , by = 'name']
  if (center)
  {
    con_signal[, pos_con_cent := as.numeric(pos_con) - as.numeric(pos_con[max_dens][1]), by = 'name']
  }
  
  # add repclass annotation
  con_signal <- merge(con_signal, unique(as.data.table(scSet@tes)[, c('name', 'class')]), all.x = TRUE, by = 'name')
  
  con_signal[which(is.na(peak)), peak := FALSE]
  
  # plot heatmap
  con_signal %>%
    ggplot(aes(pos_con_cent, name, fill = dens_scaled)) +
      geom_raster() +
      xlim(-50, 50) + 
      scale_fill_gradientn(na.value = "pink", colours = colorRampPalette(c('lightgrey', 'orange', 'red', 'darkred', 'purple4'))(256)) +
      #scale_fill_distiller(palette = 'YlGnBu', direction = 1) +
      facet_grid(class~peak, scales = 'free_y') +
      theme(legend.position  = 'top',
            legend.key.width = unit(.3,"cm"),
            legend.key.height = unit(.15,"cm"),
            axis.text        = element_blank(),
            axis.ticks       = element_blank(),
            axis.title       = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank())
}

plotGeneCov <- function(scSet,
                        min_expr = 100,
                        n_genes = 3000)
{
  protocol <- scSet@protocol
  counts   <- Repsc::counts(scSet)
  cpn_df   <- scSet@cpn
  
  counts_gene <- counts[which(type == 'gene'), ]
  
  ggdata <- 
    counts_gene[, .(tot = sum(n)), by = c('name', 'pos_con', 'peak')
      ][, gene_tot := sum(tot), by = 'name'
        ][which(gene_tot >= min_expr), ]
        
  top_genes <- 
    ggdata %>% 
    select(name, gene_tot) %>% 
    distinct %>%
    mutate(expr_group = ntile(gene_tot, n = 3)) %>%
    group_by(expr_group) %>%
    sample_n(round(n_genes / 3)) %>%
    ungroup
    
  # add info for non expressed bins
  ggdata <- 
    merge(cpn_df[which(name %in% top_genes$name), ], ggdata, all.x = TRUE, by = c('name', 'pos_con'))[which(is.na(tot)), ':=' (tot = 0, peak = FALSE)]
    
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
      inner_join(., top_genes, by = 'name') %>%
      ggplot(aes(pos_con_cent, name, fill = log10(tot + 1))) +
        geom_raster() +
        xlim(-100, 100) +
        #scale_fill_distiller(palette = 'Spectral') +
        scale_fill_gradientn(na.value = "pink", colours = colorRampPalette(c('lightgrey', 'orange', 'red', 'darkred', 'purple4'))(256)) +
        theme(legend.position  = 'top',
              axis.text        = element_blank(),
              axis.ticks       = element_blank(),
              axis.title       = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(.15,"cm"),
              panel.background = element_blank()) +
        facet_grid(expr_group~peak, scales = 'free_y')
              
  return(p)
}

plotMiri <- function(scSet)
{
  cstats <- scSet@cstats
  
  if (nrow(cstats < 5e4))
  {
    p <- 
      cstats %>% 
        ggplot(aes(perc_ribo, perc_mito, col = meta)) + 
          geom_point(size = 1, alpha = 0.25)  +
          scale_colour_viridis_d() + 
          theme(legend.position = 'none') +
          facet_wrap(~meta)
  } else {
    p <-
      cstats %>% 
        ggplot(aes(perc_ribo, perc_mito, fill = ..density..)) + 
          #geom_point(size = 1, alpha = 0.25)  +
          geom_hex(bins = 100) +
          #geom_rug(alpha = 0.1) +
          scale_fill_viridis_c() +
          theme(legend.position = 'none') +
          facet_wrap(~meta)
  }
  return(p)
}

plotUnique <- function(scSet,
                       per = c('class', 'family'),
                       label_n = 5)
{
  if (per == 'cells')
  {
    counts_te <- counts(scSet)[which(type == 'te'), ]
    ggdata    <- counts_te[, .(uniq = sum(n), all = sum(n_all)), by = c('barcode', 'meta')][, ':=' (perc_unique = uniq / all * 100)]
    
    p <- ggdata %>% ggplot(aes(perc_unique, group = meta, fill = meta)) + geom_density(alpha = 0.25) + theme(legend.position = 'none')
  }
  
  if (per == 'family')
  {
    counts <- counts(scSet)[which(type == 'te'), ]
    ggdata <- counts[, .(uniq = sum(n), all = sum(n_all)), by = c('class', 'name')][, ':=' (perc_unique = uniq / all * 100, delta = all - uniq)]
    
    p1 <-
      ggdata %>% 
        ggplot(aes(delta, perc_unique, col = class)) + 
          geom_point() + 
          scale_x_log10() + 
          scale_colour_brewer(palette = 'Set1') +
          geom_rect(xmin = 2, xmax = log10(max(ggdata$delta) + 5e5), ymin = 0, ymax = 50, fill = "transparent", col = 'black') +
          theme(legend.position = 'none')
    
    p2 <- 
      ggdata[which(delta >= 100 & perc_unique <= 50), c('class', 'name', 'uniq', 'all')] %>% 
        top_n(20, all - uniq) %>%
        mutate(name = factor(name, levels = name[order(all)], ordered = TRUE)) %>%
        tidyr::gather('mapping', 'n', all, uniq) %>% 
        ggplot(aes(name, n, fill = class, alpha = mapping)) + 
          geom_bar(stat = 'identity') + 
          scale_y_sqrt() +
          scale_fill_brewer(palette = 'Set1') +
          theme(legend.position = 'none')
          
    p <- cowplot::plot_grid(p1, p2, rel_widths = c(1,2))
  }

  return(p)
}

plotLoad <- function(scSet,
                     what = c('slim', 'full'))
{
  cstats    <- scSet@cstats
  
  if (what == 'slim')
  {
    ggdata <- 
      cstats[, .(te_load = TE_uniq / (TE_uniq + size_genic) * 100), by = c('barcode', 'meta')]
    
    # calc te_load x-axis limits
    xlimits <- c(0, quantile(ggdata$te_load, 0.995))
    
    # plot
    p <-
      ggdata %>% 
          ggplot(aes(x = te_load, y = meta, fill = meta)) + 
            ggridges::geom_density_ridges(alpha = 0.95) +
            scale_fill_viridis_d() +
            xlim(xlimits) +
            theme(legend.position = 'none')
  }
  
  if (what == 'full')
  {
    # tidy data.frame
    ggdata <- 
      bind_cols(
                cstats[, c('barcode', 'meta')], 
                cstats[, .SD, .SDcols = c('barcode', 'size_genic', grep('all|uniq', colnames(cstats), value = TRUE))
                        ][, lapply(.SD, function(x) (x + 1) / (size_genic + TE_uniq) * 100), .SDcols = -1:-2]
                ) %>%
          gather('type', 'n', -barcode, -meta) %>% 
          tidyr::separate('type', sep = '_', into = c('class', 'mapping'))
    
    p <-   
      ggdata %>% 
        ggplot(aes(meta, n, fill = class, alpha = mapping)) + 
          geom_violin() + 
          facet_grid(mapping ~ class) + 
          scale_y_log10(breaks = c(0.1,1, 10, 100)) + 
          scale_fill_brewer(palette = 'Set1') +
          theme(legend.position = 'noqne')

    # add total TE/GENE load per cell
    # ggdata <- rbindlist(list(te_load , te_load[, .(class = toupper(what_type), perc = sum(perc)), by = c('barcode', 'load')]), fill = TRUE)
   
    # p <-
      # ggdata %>%
        # ggplot(aes(class, perc, fill = class, col = load, alpha = load)) +
          # geom_violin() +
          # geom_boxplot(width=.1, outlier.colour=NA, position = position_dodge(width = 0.9), show_guide = FALSE) +
          # scale_y_log10() + 
          # scale_color_manual(values = c('black', 'black', 'black')) +
          # scale_fill_brewer(palette = 'Set1') +
          # guides(fill = 'none', colour = guide_legend(override.aes = list(fill = c('black'))))
  }

  return(p)
}

plotDistTSS <- function(scSet)
{
  tss      <- scSet@tss
  peaks    <- scSet@peaks
  bin_size <- scSet@bin_size
  
  hits <- distanceToNearest(peaks, tss)
  dt <- 
    data.table(peak_center = start(mutate(anchor_center(peaks[from(hits)]), width = 1)),
               tss_start  = start(tss[to(hits)]))[, distance := peak_center - tss_start]
  
  
  dist_n <-
   dt[, dist_bin := cut(distance, breaks = c(-1e9, seq(-250, 250, by =50), 1e9), include.lowest = TRUE)
      ][, .(n = .N), by = 'dist_bin'
        ][, perc := n / sum(n) * 100]
  #, labels = c('[0, 25]', '(25, 50]', '(50, 75]', '(75, 100]', '>100')
  p_dists <-
    dist_n %>%
      ggplot(aes(dist_bin, perc)) + 
        geom_bar(stat = 'identity')
        
  return(p_dists)
}

plotCells <- function(scSet)
{  
  # Cell calling QC
  p_csize_rank <- plotCellSize(scSet, what = 'rank')
  p_csize_hist <- plotCellSize(scSet, what = 'hist')
  p_miri       <- plotMiri(scSet)
  p_load       <- plotLoad(scSet, what = 'slim')
  p_summary    <- plotCellSummary(scSet)
  
  p1 <- plot_grid(p_csize_rank,
                     p_csize_hist,
                     p_miri,
                     p_load,
                     labels = 'auto')
                     
  p2 <- plot_grid(p_summary, labels = 'e')
  
  p_comb <- plot_grid(p1, p2, rel_heights = c(3, 1), ncol = 1)
                     # nrow = 2,
                     # rel_widths = c(3, 3, 1))
  return(p_comb)                     
}

plotFeatureComp <- function(scSet,
                            n_feats = 20,
                             type = c('gene', 'te'))
{
  gstats <- scSet@gstats
  tes    <- scSet@tes
  
  # add superfamily to gstats
  gstats <- merge(gstats, as.data.table(tes)[, c('id_unique', 'repfamily')], all.x = TRUE)
  
  ggdata <- gstats[which(type == 'te' & feature), ] %>% count(class, repfamily, name)
  
  p <-
    ggdata %>% 
      ggplot(aes(area = n, fill = class, label = name, group = class, subgroup = repfamily)) + 
        treemapify::geom_treemap() + 
        treemapify::geom_treemap_text(colour = 'black') + 
        theme(legend.position = 'none') + 
        scale_fill_brewer(palette = 'Set1') + 
        treemapify::geom_treemap_subgroup_border(colour = 'black', size = 3) + 
        treemapify::geom_treemap_subgroup_text(place = 'middle', colour = 'white', reflow = TRUE)
  
  
  
  # what_type <- type
  # gstats <- scSet@gstats[which(type == what_type), ]
  
  # if (type == 'te')
  # {
    # feats <- gstats %>% count(name, feature) %>% filter(feature) %>% top_n(25, n) %>% pull(name)

    # p <-
      # gstats[which(name %in% feats), ] %>% 
        # count(name, feature, class) %>% 
        # mutate(n = ifelse(feature, n , -1 * n)) %>% 
        # ggplot(aes(name, n, fill = class)) + 
          # geom_bar(stat = 'identity') + 
          # facet_wrap(~feature, ncol = 2, scales = 'free_x') +
          # scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3")) +
          # theme(legend.position = 'none') + coord_flip()
  # } else {
    # p <-
      # gstats %>% 
      # count(feature, class) %>% 
      # filter(feature) %>%  
      # ggplot(aes(x = feature, y = n, fill = class)) + 
        # geom_bar(stat = 'identity') + 
        # theme(legend.position = 'none')
  # }
  
  return(p)
}

plotGenes <- function(scSet)
{
  scSet <- suppressWarnings(annoClass(scSet))
  
  p_peaks         <- plotPeaks(scSet, type = 'gene', what = 'enrichment')
  p_peaks_summary <- plotPeaks(scSet, type = 'gene', what = 'summary')
  p_dist_tss      <- plotDistTSS(scSet)
  p_gene_cov      <- plotGeneCov(scSet)
  p_umis          <- plotUMIs(scSet, type = 'gene')
  p_enr_cell      <- plotEnr(scSet, type = 'gene', what = 'cells')
  p_enr_tlen      <- plotEnr(scSet, type = 'gene', what = 'transcript')
  p_varmean       <- plotVarMean(scSet, type = 'gene')

  p1 <- cowplot::plot_grid(p_peaks,
                           p_dist_tss,
                           p_peaks_summary,
                           rel_widths = c(2, 1.5, 0.75),
                           labels = c('a', 'b', 'c'),
                           nrow = 1)
  
  p2 <- cowplot::plot_grid(p_umis,
                           p_enr_cell,
                           rel_widths = c(2,1),
                           nrow =1,
                           labels = c('e', 'f'))
                           
  p3 <- cowplot::plot_grid(p1,
                           p2,
                           ncol = 1)
                           
  p2 <- cowplot::plot_grid(p3,
                           p_gene_cov,
                           rel_widths = c(1, 0.5),
                           nrow = 1,
                           labels = c('', 'd'))                        
                           
  p3 <- cowplot::plot_grid(p_enr_tlen,
                           p_varmean,
                           rel_widths = c(1, 0.5),
                           labels = c('g', 'h'))

  p4 <- cowplot::plot_grid(p2,
                           p3,
                           ncol = 1,
                           rel_heights = c(1, 0.5))
  return(p4)                           
}

plotMapping <- function(scSet)
{
  # Unique mapping
  p_unique_cells  <- plotUnique(scSet, per = 'cells')
  p_unique_fams   <- plotUnique(scSet, per = 'family')
  
  p1 <- plot_grid(p_unique_cells, p_unique_fams, rel_widths = c(1, 3), nrow = 1, labels = 'auto')
  
  if (sum(dim(scSet@pstats)) > 0)
  {
    # TE peak calling
    p_te_peaks      <- plotPeaks(scSet, type = 'te', what = 'enrichment')
    p_te_peak_stats <- plotPeaks(scSet, type = 'te', what = 'summary')  
    p_te_peak_fams  <- plotPeaks(scSet, type = 'te', what = 'family')
    
    p2 <- plot_grid(p_te_peaks, p_te_peak_fams, p_te_peak_stats, rel_widths = c(2, 2, 1), labels = c('d', 'e', 'f'))
    
    #p_load_cor      <- plotLoad(scSet, type = 'te', what = 'cor')
    #p_fams          <- plotFams(scSet)
  } else {
    p2 <- plotEmpty('No peaks called', size = 6)
  }
  
  p_umis   <- plotUMIs(scSet, type = 'te')
  p_load   <- plotLoad(scSet, type = 'te', what = 'cells')
  p_te_cov <- plotTECov(scSet)
  
  p3 <- plot_grid(p_umis, p_load, p_te_cov)
} 

plotFeatures <- function(scSet)
{
  
  if (sum(dim(scSet@gstats)) > 0)
  {
    # features
    p_varmean_gen   <- plotVarMean(scSet, type = 'gene')
    p_varmean_te    <- plotVarMean(scSet, type = 'te')
    p_featcor_hmap  <- plotCorMat(scSet, what = 'features')
    p_feat_comp     <- plotFeatureComp(scSet) 
    p_umi_hmap      <- plotUMIMat(scSet)
    p_enr_mod       <- plotEnr(scSet, what = 'module')
    
    # TE-gene comparisons
    p_famcor_hmap     <- plotCorMat(scSet, what = 'family')
    
    p1 <- plot_grid(p_varmean_gen, p_varmean_te, p_feat_comp, labels = 'auto', nrow = 1)
    p2 <- plot_grid(p_featcor_hmap, p_umi_hmap, p_enr_mod, nrow = 1, labels = c('d', 'e', 'f'))
    p3 <- plot_grid(p_famcor_hmap, labels = 'g')
    
    p_comb <- plot_grid(p1, p2, p3, ncol = 1)

  } else {
    p_comb <- plotEmpty('Run selectFeatures function before plotting', size = 6)
  }
  
  return(p_comb)
}

plotUMIs <- function(scSet,
                     type = 'te')
{
  if (sum(type %in% c('gene', 'te')) != 1) { stop ('Please define valid type') }
  what_type <- type
  counts <- Repsc::counts(annoClass(scSet))[which(type == what_type), ]
  
  if (sum(dim(scSet@pstats)) > 0)
  {
    ggdata <- 
      counts[, .(n_uniq = sum(n), n_all = sum(n_all), n_uniq_peak = sum(n[peak])), by = 'class'] %>%
        tidyr::gather('mapping', 'umis', n_uniq, n_all, n_uniq_peak)
  } else {
    ggdata <- 
      counts[, .(n_uniq = sum(n), n_all = sum(n_all)), by = 'class'] %>%
        tidyr::gather('mapping', 'umis', n_uniq, n_all)
  }
  
  p_umis <-
    ggdata %>%
      ggplot(aes(x = class, y = umis, fill = class, alpha = mapping)) +
        geom_bar(stat = 'identity', position = 'dodge') +
        theme(axis.title.y = element_blank(),
              legend.position = 'top',
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(.1,"cm")) + 
        coord_flip() +
        guides(fill = "none", alpha = "legend")
        
  if (type == 'te')
  {
    p_umis <- p_umis + scale_fill_brewer(palette = 'Set1')
  }

  return(p_umis)
}

plotUMIMat <- function(scSet)
{
  set.seed(scSet@seed)
  counts_ds  <- scSet@counts_ds
  cstats     <- scSet@cstats
  gstats     <- scSet@gstats
  
  feats_gene <- gstats[which(type == 'gene' & feature), ][order(module), ]$id_unique
  feats_te   <- gstats[which(type == 'te' & feature), ][order(module), ]$id_unique
  
  # only retain cstats for cells surviving downsampling
  cstats_f <- cstats[which(barcode %in% unique(counts_ds$barcode)), ]
  
  # sample cells if too many
  if (nrow(cstats_f > 1e4))
  {
    n_meta <- length(unique(cstats$meta))
    cstats_f <- cstats_f[, .SD[sample(.N, min(1e4 / n_meta, .N))], by = 'meta']
  }
  
  # filter ds counts for features and sampled cells
  counts_ds_f <- counts_ds[which(id_unique %in% feats_te & barcode %in% cstats_f$barcode), ]
  
  # create mat
  mat = as.matrix(longToSparse(counts_ds_f[, c('id_unique', 'barcode', 'n_ds')]))
  
  # cluster cells
  hc_cells <- hclust(tgs_dist(tgs_cor(log10(1 + 7 * mat))), 'ward.D2')
  
  # reorder factors according to clusters
  counts_ds_f[, ':=' (id_unique = factor(id_unique, levels = feats_te, ordered = TRUE), barcode = factor(barcode, levels = hc_cells$labels[hc_cells$order], ordered = TRUE))]

  # plot
  p <- 
    counts_ds_f %>% 
      ggplot(aes(barcode, id_unique, fill = log10(n_ds))) + 
        geom_raster() +
        theme(legend.position = 'none',
        legend.key.width = unit(.3,"cm"),
        legend.key.height = unit(.15,"cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  
  return(p)
}

plotVarMean <- function(scSet,
                        n_label = 10,
                        type = c('gene', 'te'))
{
  if (sum(dim(scSet@gstats)) > 0)
  {

    if (type == 'gene')
    {
      gstats    <- scSet@gstats[which(type == 'gene'), ]
    }
    if (type == 'te' | type == 'te_fam')
    {
       gstats    <- scSet@gstats[which(type %in% c('te', 'te_fam')), ]
    }
    
    if (is.null(gstats$class))
    {
      gstats$class <- ''
    }
    
    n_feats <- nrow(gstats[which(feature), ])  
    
    #which(class == 'protein_coding')
    gstats[, above := (log2(varmean) - vm_trend) >= sort(log2(varmean) - vm_trend, decreasing = TRUE)[n_label]
      ][which(type == 'te_fam'), fontface := 'bold'
        ][which(type != 'te_fam'), fontface := 'italic']
    
    p_varmean <-
      gstats %>% 
        ggplot(aes(log2(mean_ds), log2(varmean), size = feature, alpha = feature, shape = type, col = class)) + 
          scale_size_manual(values  = c(0.25, 1.5)) +
          scale_alpha_manual(values = c(0.25, 0.75)) +
          scale_shape_manual(values = c(20, 8)) +
          geom_point() +
          ggrepel::geom_text_repel(data = subset(gstats, above), size = 2, alpha = 1, aes(label = id_unique, fontface = fontface)) + 
          ggtitle(paste0('n = ', n_feats)) + 
          theme(plot.title = element_text(margin = margin(b = -25)),
                legend.position = 'none')
       
    if (type != 'gene') 
    {
      p_varmean <- p_varmean + scale_color_brewer(palette = 'Set1')
    }
  } else {
    p_varmean <- plotEmpty(string = 'compute Gstats first')
  }
  return(p_varmean)
}

plotCellSummary <- function(scSet)
{
  cstats <- scSet@cstats
  
  # ggdata <- 
    # cstats[, .(all = length(unique(barcode)), 
               # pass_size = sum(pass_size), 
               # pass_miri_size = sum(is_cell, na.rm = TRUE)), by = 'meta'] %>% 
    # tidyr::gather('pass', 'n', -meta) %>% 
    # dplyr::rename(barcodes = n) %>%
    # group_by(meta) %>%
    # mutate(perc_thrown = (1 - round(barcodes / lag(barcodes), 4)) * 100, 
           # n_thrown    = lag(barcodes) - barcodes) %>% 
    # tidyr::replace_na(list(perc_thrown = 0, n_thrown = 0))      
  
  # p <-  gridExtra::tableGrob(ggdata %>% filter(pass != 'all') %>% arrange(meta), theme = gridExtra::ttheme_default(base_size = 6)) #+facet_wrap(~meta)
  
  p <- 
    gridExtra::tableGrob(
                         cstats[, .('valid barcodes (#)' = .N, 'genic umis (median)' = median(size_genic), 'TE umis (median)' = median(TE_uniq)), by = 'meta'],
                         theme = gridExtra::ttheme_default(base_size = 8)
                         )
  
  return(p)
}

plotPeaks <- function(scSet,
                      type = c('gene', 'te'),
                      what = c('enrichment', 'summary'))
{
  what_type <- type
  pstats    <- annoClass(scSet)@pstats[which(type == what_type), ]

  if (what == 'enrichment')
  {
    p_peaks <-
      pstats[which(peak), ] %>%
        ggplot(aes(n_samp_neighb, fc, col = class)) + 
          geom_point(size = 1, alpha = 0.25) +
          scale_x_log10() +
          theme(legend.position = 'none')
          
    p_peaks <- 
      p_peaks + 
        stat_binhex(data = pstats[which(!peak), ], aes(x=n_samp_neighb, y=fc, alpha = sqrt(..count..)), fill = 'black', inherit.aes = FALSE, bins = 100) +
        scale_alpha_continuous(range = c(0, 1))
        
    if (type == 'te') { p_peaks <- p_peaks + scale_color_brewer(palette = 'Set1') }
  }

  if (what == 'summary')
  {
    # n peaks and fams
    p_peaks <-
      pstats %>% 
        count(class, name, peak) %>% 
        filter(peak) %>% 
        group_by(class) %>% 
        summarize(peaks = sum(n), 
                  fams = length(class)) %>% 
        tidyr::gather('type', 'n', -class) %>% 
        ggplot(aes(x = type, y = n, fill = class)) + 
          geom_bar(stat = 'identity') + 
          theme(legend.position = 'none')
    
    if (type == 'te') { p_peaks <- p_peaks + scale_fill_brewer(palette = 'Set1') }
  }
  
  if (what == 'family')
  {
    p_peaks <-
      pstats[which(type == 'te'), 
        ][, .(peak = sum(n_tot[peak]), all = sum(n_tot)), by = c('class', 'name')
          ][, ':=' (perc_peak = peak / all * 100, delta = all - peak)] %>% 
      ggplot(aes(delta, perc_peak, col = class)) +
        geom_point() +
        scale_x_log10() +
        scale_colour_brewer(palette = 'Set1') +
        theme(legend.position = 'none')
  }

  return(p_peaks)  
}
