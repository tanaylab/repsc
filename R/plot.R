plotCellSize <- function(scSet,
                         what = c('hist', 'rank', 'scatter'))
{
  cstats <- scSet@cstats[which(size_genic >= 100), ]
  
  if (sum(dim(cstats)) == 0)
  {
    cstats <- Repsc::cstats(scSet)    
  }
 
  if (what == 'hist')
  {  
    if (length(unique(cstats$meta)) > 1)
    {
      p <-
        cstats %>%
          ggplot(aes(x = size_genic, y = meta, fill = meta)) + 
            ggridges::geom_density_ridges(alpha = 0.95, quantile_lines = TRUE, quantiles = 2) +
            scale_fill_viridis_d() + 
            theme(legend.position = 'none')
    } else {
      p <-
        cstats %>% 
          ggplot(aes(x = size_genic, fill = meta)) + 
            geom_density() + 
            scale_fill_viridis_d() + 
            theme(legend.position = 'none')
    }
    
    # decide of log scale or not
    dens <- density(cstats$size_genic)
    if (dens$x[which.max(dens$y)] < 5000)
    {
      p <- p + scale_x_log10()
    } else {
      p <- p + coord_cartesian(xlim = c(min(cstats$size_genic), min(max(cstats$size_genic), 7.5e4)))
    }
    
    p <- p + theme(axis.title.y = element_blank())
  }
  
  if (what == 'rank')
  {
    p <-
      cstats %>%
        ggplot(aes(rank, size_genic, col = meta)) + 
          geom_point(size = 0.1) + 
          scale_x_log10() + 
          scale_y_log10() +
          #geom_hline(yintercept = size_t) +
          #scale_size_manual(values = c(8, 0.01)) + 
          scale_colour_viridis_d() + 
          theme(legend.position = 'none') 
  }
  
  if (what == 'scatter')
  {
    p <-
      cstats %>%
        ggplot(aes(size_genic, n_genes, color = meta)) + 
          geom_point(size = 0.1, alpha = 0.25) +
          scale_colour_viridis_d() +
          scale_x_log10() +
          scale_y_log10() +
          #facet_wrap(~meta) +
          theme(legend.position = 'none')
  }
  
  return(p)
}

plotCorMat <- function(scSet, 
                       what = 'features',
                       min_cor = NULL,
                       min_expr = 100)
{
  set.seed(scSet@seed)
  scSet     <- suppressWarnings(annoClass(scSet))
  gstats    <- scSet@gstats
  counts    <- Repsc::counts(scSet)
  counts_ds <- scSet@counts_ds
  tes       <- scSet@tes
 
    
  if (what == 'features')
  {
    gstats_f  <- gstats[which(type == 'te' & feature), ]
    
    if (nrow(gstats_f) > 1e3)
    {
      n_class <- length(unique(gstats_f$class))
      gstats_f <- gstats_f[, .SD[sample(.N, min(1e3 / n_class, .N))], by = 'class']
    }
    
    feats <- gstats_f$id_unique
        
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
    #cors_df <- merge(cors_df, unique(gstats[, c('id_unique', 'class')]), all.x = TRUE, by = 'id_unique')
    
    # plot
    steps <- as.numeric(cumsum(table(gstats_f[which(type == 'te'),]$module))) + 0.5
    
    # cap cor valid
    cap <- 0.2
    cors_df[which(cor > 0.2), cor := cap]
    cors_df[which(cor < -0.2), cor := -cap]
    
    p_cor_hmap <-
      cors_df %>% 
          ggplot(aes(id_unique, id_unique2, fill = cor)) + 
            geom_raster() +
            #scale_fill_viridis_c() +
            scale_fill_gradientn(colors = colorspace::diverging_hcl(palette = 'Berlin', n = 256), limits = c(-cap, cap)) +
            geom_hline(yintercept = steps, lwd = 0.25, col = 'white') +
            geom_vline(xintercept = steps, lwd = 0.25, col = 'white') + 
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
    
    # filter for feature fams/genes only
    counts_ds_f <- 
      counts_ds[which(name %in% gstats[which(feature & (type == 'te_fam' | class == 'protein_coding')), ]$name), ]
    
    # vectors of TE families and protein coding genes for later subsetting
    fams  <- unique(counts_ds_f[which(type == 'te'), ]$name)
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
    if (is.null(min_cor)) { min_cor <- quantile(cor_mat_fvsg, 0.95) }
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
    labels_gene <- labels_gene[nchar(labels_gene) < 10] 
    labels_gene <- labels_gene[seq(1, length(labels_gene), l = min(length(labels_gene), 40))]
    labels_fam <- levels(fg_cors_df$fam)
    labels_fam <- unique(labels_fam[seq(1, length(labels_fam), l = 20)])
    
    # cap cor
    cap <- 0.4
    fg_cors_df[which(cor > 0.4),  cor := cap]
    fg_cors_df[which(cor < -0.4), cor := -cap]
    
    # plot
    
    # color fam labels by class
    name_class_rel <- 
      merge(data.table(name = levels(fg_cors_df$fam)), 
            unique(as.data.table(tes[tes$name %in% levels(fg_cors_df$fam)])[, c('name', 'class')]), 
            all.x = TRUE)[order(match(name, levels(fg_cors_df$fam))),]
    colrs <- c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3")
    
    p_cor_hmap <-
      fg_cors_df %>%
          ggplot(aes(gene, fam, fill = cor)) + 
            geom_raster() +
            scale_x_discrete(breaks = as.character(labels_gene)) + 
            scale_y_discrete(breaks = as.character(labels_fam)) +
            #scale_fill_viridis_c() + 
            scale_fill_gradientn(colors = colorspace::diverging_hcl(palette = 'Berlin', n = 256), limits = c(-cap, cap)) +
            theme(legend.position = 'none',
                  axis.text.y = element_text(size = 8, colour = colrs[name_class_rel[which(name %in% labels_fam), ]$class])) +
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
    # p  <- plot_grid(p1, p2, rel_widths = c(4, 1))
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
    tes    <- scSet@tes
    
    # add superfamily and replace Sines by repfamily
    gstats <- merge(gstats, unique(data.table(name = tes$name, repfamily = tes$repfamily)), all.x = TRUE, by = 'name')
    gstats[which(class == 'SINE'), ]$name <- gstats[which(class == 'SINE'), ]$repfamily
    
    # Fishers
    pop_size    <- nrow(gstats)
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
    
    # select fams to plot
    ggdata_subs <- ggdata %>% filter(qval < 0.05 & inmodule > 2)
    if (length(unique(ggdata_subs$name)) > 20) 
    {
      fams_per_mod <- floor(20 / max(ggdata$module, na.rm = TRUE))
      fams_to_plot <- ggdata %>% filter(qval < 0.05 & inmodule > 2) %>% group_by(module) %>% top_n(fams_per_mod, -qval) %>% ungroup %>% pull(name)
      ggdata_subs  <- ggdata %>% filter(name %in% fams_to_plot)
    }
    
    # plot
    p <-
      ggdata_subs %>%
      ggplot(aes(factor(module, levels = as.character(sort(as.numeric(unique(ggdata_subs$module))))), name, fill = -log10(qval))) + 
        geom_raster() +
        geom_text(aes(label = inmodule, col = -log10(qval)), size = 2.25) + 
        xlab('') + 
        scale_fill_distiller(palette = 'BuPu', direction = 1) +
        scale_color_distiller(palette = 'PuRd', direction = -1) +
        theme(legend.position = 'none',
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(1,"cm"))
  }
  
  return(p)
}               

plotTECov <- function(scSet, 
                      what     = 'heatmap',
                      pattern  = NULL,
                      center   = TRUE,
                      min_expr = 100)
{
  protocol  <- scSet@protocol
  counts    <- Repsc::counts(scSet)
  counts_te <- counts[which(type == 'te'),] 
  cpn_df    <- scSet@cpn
  
  # aggregate counts across cells and TE loci
  counts_aggr <- aggr(counts_te, min_expr = min_expr)
  
  # only retain expressed fams
  cpn_df <- cpn_df[which(name %in% unique(counts_aggr$name)), ]
  
  # calc density
  con_signal <- 
    merge(cpn_df, counts_aggr, all.x = TRUE, by = c('name', 'pos_con'))[which(is.na(n_tot)), n_tot := 0
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
  con_signal[, pos_con_cent := 1:.N - which(max_dens)[1], by = 'name']
  
  # add repclass annotation
  con_signal <- merge(con_signal, unique(as.data.table(scSet@tes)[, c('name', 'class')]), all.x = TRUE, by = 'name')
  
  con_signal[which(is.na(peak)), peak := FALSE]
  
  # cluster
  smat     <- longToSparse(con_signal[, c('name', 'pos_con', 'dens_scaled')])
  mat      <- as.matrix(smat[, Matrix::colSums(smat) != 0])
  cor_mat  <- tgs_cor(t(mat))
  dist_mat <- tgs_dist(mat)#as.dist(sqrt(2 * (1 - cor_mat)))
  hc       <- hclust(dist_mat, 'ward.D2')
  con_signal[, name := factor(name, levels = hc$labels[hc$order], ordered = TRUE)]
  
  if (what == 'heatmap')
  {
    # plot heatmap
    p1 <-
      con_signal %>%
        ggplot(aes(pos_con, name, fill = dens_scaled)) +
          geom_raster() +
          #xlim(-100, 100) + 
          scale_fill_gradientn(na.value = "pink", colours = colorRampPalette(c('lightgrey', 'orange', 'red', 'darkred', 'purple4'))(256)) +
          #scale_fill_distiller(palette = 'YlGnBu', direction = 1) +
          facet_grid(rows = vars(class), scales = 'free', space = 'free') +
          theme(legend.position  = 'none',
                legend.key.width = unit(.3,"cm"),
                legend.key.height = unit(.15,"cm"),
                strip.text = element_blank(),
                axis.text.x      = element_blank(),
                axis.text        = element_blank(),
                axis.ticks       = element_blank(),
                axis.title       = element_blank(),
                panel.spacing    = unit(0.25, 'lines'),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank())

    p2 <-
      con_signal %>%
        ggplot(aes(pos_con_cent, name, fill = dens_scaled)) +
          geom_raster() +
          xlim(-20, 20) + 
          scale_fill_gradientn(na.value = "pink", colours = colorRampPalette(c('lightgrey', 'orange', 'red', 'darkred', 'purple4'))(256)) +
          #scale_fill_distiller(palette = 'YlGnBu', direction = 1) +
          facet_grid(class ~ peak, scales = 'free', space = 'free') +
          theme(legend.position  = 'none',
                legend.key.width = unit(.3,"cm"),
                legend.key.height = unit(.15,"cm"),
                strip.text.x      = element_blank(),
                axis.text.x      = element_blank(),
                axis.text        = element_blank(),
                axis.ticks       = element_blank(),
                axis.title       = element_blank(),
                panel.spacing    = unit(0.25, 'lines'),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank())
                
    p <- plot_grid(p1,p2, nrow = 1, rel_widths = c(1, 1))
  }
  
  if (what == 'histogram')
  {
    # bin by size
    con_signal <- con_signal[, con_size := .N, by = 'name'][, size_bin := cut(con_size, c(0, 50, 100, 200, Inf), include.lowest = TRUE, labels = FALSE)]
    
    if (sum(dim(scSet@gstats)) > 0)
    {
      fams_to_highl <- scSet@gstats[which(type == 'te_fam' & feature), ]$name
      fams_to_highl <- sample(fams_to_highl, min(40, length(fams_to_highl)))
    } else {
      fams_to_highl <-
        con_signal[, .(dens_max = max(dens), class = class[1], size_bin = size_bin[1]), by = 'name'] %>% group_by(size_bin) %>% top_n(10, dens_max) %>% pull(name)
      # fams_to_highl <-
        # con_signal %>% group_by(size_bin) %>% sample_n(15) %>% pull(name)
    }
   
    con_sig_f <- con_signal[which(name %in% fams_to_highl), ]

    # order by dens position
    con_sig_f[, name := factor(name, levels = con_sig_f[, .(dens_max = which.max(dens)), by = 'name'][order(dens_max),]$name, ordered = TRUE)]
    
    
    # create data.frame of TE/genome boundaries to add lines
    borders_df <- 
      merge(con_sig_f[, .(border = sum(!grepl('-', pos_con, fixed = TRUE))), by = 'name'],
            con_sig_f[, c('name', 'size_bin')], all.x = TRUE)
    colrs <- 
      c('DNATRUE' = "#E41A1C", 'LINETRUE' = "#377EB8", 'LTRTRUE' = "#4DAF4A", 'SINETRUE' = "#984EA3", 'DNAFALSE' = alpha("#E41A1C", 0.25), 'LINEFALSE' = alpha("#377EB8", 0.25), 'LTRFALSE' = alpha("#4DAF4A", 0.25), 'SINEFALSE' = alpha("#984EA3", 0.25))    
    
    test = con_sig_f[, pos_con_num := 1:.N, by = 'name'][, test := as.factor(paste0(class, peak))]
    
    # p <-    
      # test %>%
        # ggplot(aes(x = pos_con_num, y = name, height = dens_scaled, fill = as.numeric(test))) +
        # ggridges::geom_density_ridges_gradient(stat = 'identity', scale = 1) +
        # geom_segment(data = borders_df, aes(x = border, xend = border, y = as.numeric(as.factor(name)),
                                      # yend = as.numeric(as.factor(name)) + .9), inherit.aes = FALSE) +
        
        # test %>% ggplot(aes(x = pos_con_num, y = dens_scaled, fill = test)) + geom_bar(stat = 'identity') + facet_wrap(~name, ncol = 1) +      
        # theme(legend.position = 'none',
              # strip.text = element_blank(),
              # panel.spacing.y = unit(-0.05, "lines"),
              # axis.ticks.y = element_blank(),
              # axis.text.y  = element_blank(),
              # axis.title.y = element_blank(),
              # axis.title.x = element_blank(),
              # axis.ticks.x = element_blank(),
              # axis.text.x = element_blank()) +
        # scale_fill_gradientn(colours = colrs[levels(test$test)]) +
        # facet_wrap(~size_bin, scales = 'free', ncol = 1) 
  }

  return(p)
}

plotFAM <- function(scSet, fam = NULL)
{
  fam <- 'MT-int'
  
  pstats <- scSet@pstats[which(name == fam),]
  
  # consensus coverage before and after normalization
  p_umi_cov <- 
    pstats %>% select(pos_con, n_tot, n_samp) %>% tidyr::gather('type', 'n', -pos_con) %>% ggplot(aes(pos_con, n, fill = type)) + geom_bar(stat = 'identity', position = 'dodge')
    
  # loci expressed
  
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
  cstats <- scSet@cstats[which(size_genic >= 100), ]
  
  if (nrow(cstats < 5e4))
  {
    p <- 
      cstats %>% 
        ggplot(aes(perc_ribo, perc_mito, col = meta)) + 
          geom_point(size = 0.25, alpha = 0.25)  +
          scale_colour_viridis_d() + 
          theme(legend.position = 'none') +
          facet_wrap(~meta) +
          theme(panel.margin = unit(0, "lines"))
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
                       what = 'te',
                       highlight = FALSE,
                       highlight_n = 50)
{
  mstats <- scSet@mstats[which(type == what), ]
   
  if (!highlight)
  {  
    p <-
      mstats %>% 
        ggplot(aes(all, perc_unique, col = class)) + 
          geom_point(size = 0.25, alpha = 0.5) + 
          scale_x_log10() + 
          scale_colour_brewer(palette = 'Set1') +
          annotate('rect', xmin = 100, xmax = Inf, ymin = 0, ymax = 50, fill = "lightgrey", col = 'grey', alpha = 0.3) +
          theme(legend.position = 'none')
  } else {
    # filt Rpl/Gm genes
    if (what == 'gene')
    {
      organism      <- scSet@genome@common_name
      ribo_prefixes <- c('Human' = '^RP[LS]', 'Mouse' = '^Rp[ls]')
      ribo_prefix   <- ribo_prefixes[organism]
      mstats        <- mstats[which(!grepl(ribo_prefix, name)), ]
      mstats        <- mstats[which(!grepl('^Gm[1-9]', name)), ]
    }
    
    # filt for worst cases
    if (nrow(mstats[which(delta >= 100 & perc_unique <= 50),]) > 0)
    {
      mstats <- mstats[which(delta >= 100 & perc_unique <= 50),]
    }
    
    # filt for protein-coding if gene
    if (what == 'gene')
    {
      mstats <- mstats[which(class == 'protein_coding'), ]
    }
    
    # order by all n expr, filter, transform
    mstats_t <- 
      mstats %>%
      top_n(highlight_n, fc) %>%
      mutate(name = factor(name, levels = name[order(perc_unique, decreasing = TRUE)], ordered = TRUE)) %>%
      mutate(uniq = perc_unique) %>%
      select(name, uniq, all, class) %>% 
      tidyr::gather('mapping', 'n', -name, -class) %>% 
      mutate(n = ifelse(mapping == 'all', log10(n), n * -1)) %>%
      as.data.table()
      
    # correct log transform for 0 counts
    mstats_t <- mstats_t[which(!is.finite(n)), n := 0]
    
    # color vect by class
    if (what == 'te')
    {
      colrs <- c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3")
    } else {
      colrs <- structure(RColorBrewer::brewer.pal(4, 'Set2'), names = unique(mstats_t$class))
    }
    
    p <-
      mstats_t %>% 
        ggplot(aes(x = name, y = n, alpha = mapping, fill = class)) +
          geom_bar(data = mstats_t[which(mapping == 'uniq'), ], stat = "identity") +
          geom_bar(data = mstats_t[which(mapping == 'all'), ], stat = "identity") +
          scale_y_continuous(breaks = c(seq(-100, 0, 1), seq(1, 10, 1)) , labels = as.character(c(rev(seq(0, 100, by = 1)), 10^c(1:10)))) +
          scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3", 'protein_coding' = 'black')) +
          theme(legend.position = 'none',
                plot.background = element_rect(fill = "#D3D3D34D"),
                #panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) +
                #axis.text.y = element_text(colour = colrs[unique(mstats_t[order(name), c('class', 'name')])$class])) +
          coord_flip()
  }

  return(p)
}

plotLoad <- function(scSet,
                     what = c('slim', 'scatter', 'full'))
{
  cstats    <- scSet@cstats[which(size_genic >= 100), ]
  
  if (what == 'slim')
  {
    ggdata <- 
      cstats[, .(te_load = TE_uniq / (TE_uniq + size_genic) * 100), by = c('barcode', 'meta')]
    
    # calc te_load x-axis limits
    xlimits <- c(0, quantile(ggdata$te_load, 0.995, na.rm = TRUE))
    
    # plot
    p <-
      ggdata %>% 
          ggplot(aes(x = te_load, y = meta, fill = meta)) + 
            ggridges::geom_density_ridges(alpha = 0.95, quantile_lines = TRUE, quantiles = 2) +
            scale_fill_viridis_d() +
            xlim(xlimits) +
            theme(legend.position = 'none',
                  axis.title.y = element_blank())
  }
  
  if (what == 'scatter')
  {
    p <-
      cstats %>% 
        ggplot(aes(size_genic, y = TE_uniq / (TE_uniq + size_genic) * 100, col = meta)) + 
          geom_point(size = 0.1, alpha = 0.25) + 
          scale_colour_viridis_d() + 
          scale_x_log10() + 
          #facet_wrap(~meta) +
          theme(legend.position = 'none')  
  }
  
  if (what == 'full')
  {
    ggdata <- 
      cstats[, grepl('barcode|meta|size_genic|uniq|peak', colnames(cstats)), with = FALSE] %>% 
      tidyr::gather('class', 'n', -c('barcode', 'meta', 'size_genic')) %>% 
      mutate(load = n / (size_genic + n) * 100) %>%
      as.data.table()

    ggdata <- tidyr::separate(ggdata, 'class', into = c('class', 'type'), sep = '_')
  
    p <-   
      ggdata %>% 
        ggplot(aes(x = class, y = load, fill = class, alpha = type)) + 
          geom_split_violin() + 
          #facet_wrap(~meta) +
          scale_y_log10(breaks = c(0.1,1, 10, 100)) + 
          scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3", 'TE' = "orange")) +
          theme(legend.position = 'none')
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

#' Plot cell summary QC
#' @export
plotCells <- function(scSet)
{  
  # trim meta label length
  cstats           <- scSet@cstats
  cstats_mod       <- cstats
  meta_levels      <- unique(cstats$meta)
  
  # decide till what position to cut and substring
  pos <-min(which(sapply(1:100, function(x) length(unique(substring(meta_levels, 1, x)))) == length(meta_levels)))
  cstats_mod$meta  <- substring(cstats$meta, 1, pos) 
  scSet@cstats     <- cstats_mod
  
  # Cell calling QC
  p_csize_rank   <- plotCellSize(scSet, what = 'rank')
  p_csize_hist   <- plotCellSize(scSet, what = 'hist')
  p_umi_vs_gen   <- plotCellSize(scSet, what = 'scatter')
  p_miri         <- plotMiri(scSet)
  p_load         <- plotLoad(scSet, what = 'slim')
  p_load_vs_size <- plotLoad(scSet, what = 'scatter')
  p_n_barcodes   <- 
    ggplot(scSet@cstats[, .(n_cells = .N), by = 'meta'], aes(meta, n_cells, fill = meta)) + geom_bar(stat = 'identity') + scale_fill_viridis_d() + theme(legend.position = 'none', axis.title.x = element_blank())
  p_summary      <- plotCellSummary(scSet)
  
  p1 <- plot_grid(p_csize_rank,
                  p_csize_hist,
                  p_miri,
                  p_load,
                  p_load_vs_size,
                  p_umi_vs_gen,
                  labels = 'auto')
                     
  p2 <- plot_grid(p_n_barcodes, p_summary, labels = c('g', 'e'), rel_widths = c(1, 3))
  
  p_comb <- plot_grid(p1, p2, rel_heights = c(3, 1), ncol = 1)

  return(p_comb)                     
}

plotFeatureComp <- function(scSet,
                            n_feats = 20,
                             type = c('gene', 'te'))
{
  gstats <- scSet@gstats
  tes    <- scSet@tes
  tes_f  <- tes[tes$name %fin% unique(gstats$name)]
  
  # add superfamily to gstats
  gstats <- merge(gstats, unique(as.data.table(tes_f)[, c('name', 'repfamily')]), by = 'name', all.x = TRUE)[which(is.na(repfamily)), repfamily := '?']
  ggdata <- gstats[which(type == 'te' & feature), ] %>% count(class, repfamily, name)
  
  p <-
    ggdata %>%  
      ggplot(aes(area = n, fill = class, label = name, group = class, subgroup = repfamily)) + 
        treemapify::geom_treemap() + 
        treemapify::geom_treemap_text(colour = 'black') + 
        theme(legend.position = 'none') + 
        scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3")) + 
        treemapify::geom_treemap_subgroup_border(colour = 'black', size = 3) + 
        treemapify::geom_treemap_subgroup_text(place = 'middle', colour = 'white', reflow = TRUE)
  
  return(p)
}

# plotGenes <- function(scSet)
# {
  # scSet <- suppressWarnings(annoClass(scSet))
  
  # p_peaks         <- plotPeaks(scSet, type = 'gene', what = 'enrichment')
  # p_peaks_summary <- plotPeaks(scSet, type = 'gene', what = 'summary')
  # p_dist_tss      <- plotDistTSS(scSet)
  # p_gene_cov      <- plotGeneCov(scSet)
  # p_umis          <- plotUMIs(scSet, type = 'gene')
  # p_enr_cell      <- plotEnr(scSet, type = 'gene', what = 'cells')
  # p_enr_tlen      <- plotEnr(scSet, type = 'gene', what = 'transcript')
  # p_varmean       <- plotVarMean(scSet, type = 'gene')

  # p1 <- plot_grid(p_peaks,
                           # p_dist_tss,
                           # p_peaks_summary,
                           # rel_widths = c(2, 1.5, 0.75),
                           # labels = c('a', 'b', 'c'),
                           # nrow = 1)
  
  # p2 <- plot_grid(p_umis,
                           # p_enr_cell,
                           # rel_widths = c(2,1),
                           # nrow =1,
                           # labels = c('e', 'f'))
                           
  # p3 <- plot_grid(p1,
                           # p2,
                           # ncol = 1)
                           
  # p2 <- plot_grid(p3,
                           # p_gene_cov,
                           # rel_widths = c(1, 0.5),
                           # nrow = 1,
                           # labels = c('', 'd'))                        
                           
  # p3 <- plot_grid(p_enr_tlen,
                           # p_varmean,
                           # rel_widths = c(1, 0.5),
                           # labels = c('g', 'h'))

  # p4 <- plot_grid(p2,
                           # p3,
                           # ncol = 1,
                           # rel_heights = c(1, 0.5))
  # return(p4)                           
# }

plotMapping <- function(scSet,
                        type = 'te')
{
  counts <- scSet@counts
  if (sum(counts$n, na.rm = TRUE) != sum(counts$n_all, na.rm = TRUE))
  {
    # UMI loss
    p_umis <- plotUMIs(scSet, type = type, what = 'mapping') # slow, can be improved
    
    # Unique mapping
    p_uniq_scat   <- plotUnique(scSet, what = type)
    p_uniq_highl  <- plotUnique(scSet, what = type, highlight = TRUE)
    
    p1 <- plot_grid(p_umis,
                    p_uniq_scat,
                    labels = 'auto',
                    rel_heights = c(1.5, 1),
                    ncol = 1)
    
    p2 <- plot_grid(p_uniq_highl, labels = 'c')
    
    p_comb <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.5))
  } else {
    p_comb <- plotEmpty('No multi-mapper found in alignment file(s)', size = 6)
  }
  
  return(p_comb)
}

#' Plot feature summary QC
#' @export
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
    p_enr_mod       <- plotEnr(scSet, type = 'te', what = 'module')
    
    # TE-gene comparisons
    p_famcor_hmap     <- plotCorMat(scSet, what = 'family')
    
    p1 <- plot_grid(p_varmean_gen, p_varmean_te, p_feat_comp, labels = 'auto', nrow = 1)
    p2 <- plot_grid(p_featcor_hmap, p_umi_hmap, p_enr_mod, nrow = 1, labels = c('d', 'e', 'f'), rel_widths = c(1.5, 2, 2))
    p3 <- plot_grid(p_famcor_hmap, labels = 'g')
    
    p_comb <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1, 1, 1.5))

  } else {
    p_comb <- plotEmpty('Run selectFeatures function before plotting', size = 6)
  }
  
  return(p_comb)
}

plotUMIs <- function(scSet,
                     type = 'te',
                     what = 'mapping')
{
  if (sum(type %in% c('gene', 'te')) != 1)      { stop ('Type must be "te" or "gene"') }
  if (sum(what %in% c('mapping', 'peaks')) != 1) { stop ('What must be "mapping" or "peaks"') }
  what_type <- type
  
  
  if (what == 'peaks')
  {
    counts <- Repsc::counts(annoClass(scSet))[which(type == what_type), ]
    
    ggdata <- 
      counts[, .(n_uniq = sum(n), n_uniq_peak = sum(n[peak])), by = 'class'] %>%
        tidyr::gather('mapping', 'umis', n_uniq, n_uniq_peak)
  }
  
  if (what == 'mapping')
  {
    mstats <- scSet@mstats[which(type == what_type), ]
    ggdata <- 
      mstats[, .(all = sum(all), uniq = sum(uniq)), by = 'class'] %>%
        mutate(perc = round(uniq / all * 100, 2)) %>%
        mutate(all  = all - uniq) %>%
        tidyr::gather('mapping', 'umis', uniq, all, -perc) %>% 
        mutate(perc = ifelse(mapping == 'all', '', perc)) %>%
        as_tibble()
  }
  
  p_umis <-
    ggdata %>%
      ggplot(aes(x = class, y = umis, fill = class, alpha = mapping)) +
        geom_bar(stat = 'identity') +
        geom_text(aes(label = perc, vjust = 0)) +
        theme(axis.title.y = element_blank(),
              legend.position = 'top',
              legend.key.width = unit(.3,"cm"),
              legend.key.height = unit(.1,"cm")) + 
        guides(fill = "none", alpha = "legend")
        
  if (type == 'te')
  {
    p_umis <- p_umis + scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3"))
  }

  return(p_umis)
}

plotUMIMat <- function(scSet)
{
  set.seed(scSet@seed)
  counts_ds  <- scSet@counts_ds
  cstats     <- scSet@cstats
  gstats     <- scSet@gstats
  
  gstats_f  <- gstats[which(type == 'te' & feature), ]
  if (nrow(gstats_f) > 1e3)
  {
    n_class <- length(unique(gstats_f$class))
    gstats_f <- gstats_f[, .SD[sample(.N, min(1e3 / n_class, .N))], by = 'class']
  }
  feats_te   <- gstats_f[which(type == 'te' & feature), ][order(module), ]$id_unique
  
  # only retain cstats for cells surviving downsampling
  cstats_f <- cstats[which(barcode %in% unique(counts_ds$barcode)), ]
  
  # sample cells if too many per meta
  if (nrow(cstats_f > 1e4))
  {
    n_meta <- length(unique(cstats$meta))
    cstats_f <- cstats_f[, .SD[sample(.N, min(1e4 / n_meta, .N))], by = 'meta']
  }
  
  # filter ds counts for features and sampled cells
  counts_ds_f <- counts_ds[which(id_unique %in% feats_te & barcode %in% cstats_f$barcode), ]
  
  # create mat
  mat = longToSparse(counts_ds_f[, c('id_unique', 'barcode', 'n_ds')])
  
  # transform
  mat_trans                 <- as.matrix(log2(1 + 7 * mat))

  # cluster cells per meta
  hc_cells <- 
    tapply(1:ncol(mat_trans), sapply(strsplit(colnames(mat_trans), '|', fixed = TRUE), function(x) x[2]), function(x)
    {
      hc <- hclust(tgs_dist(tgs_cor(mat_trans[, x])), 'ward.D2')
      hc$labels[hc$order]
    }, simplify = FALSE)
  hc_cells <- as.character(unlist(hc_cells))
  
  # scale
  mat_trans                 <- mat_trans  - apply(mat_trans, 1, median)
  mat_trans[mat_trans  < 0] <- 0
  
  # convert to long and smooth
  counts_ds_f_smoo <- 
    as.data.table(melt(mat_trans, varnames = c('id_unique', 'barcode'), value.name = 'n_ds'))[, rollmean := frollmean(n_ds, n = 15, fill = 0, align = "center")
    ][which(rollmean > 0), ]
    
  # add meta column
  #counts_ds_f_smoo <- merge(counts_ds_f_smoo, cstats[, c('barcode', 'meta')], all.x = TRUE, by = 'barcode')
  
  # reorder factors according to clusters
  counts_ds_f_smoo[, ':=' (id_unique = factor(id_unique, levels = feats_te, ordered = TRUE), barcode = factor(barcode, levels = hc_cells, ordered = TRUE))]

  # plot
  steps      <- as.numeric(cumsum(table(gstats_f[which(type == 'te'),]$module))) + 0.5
  steps_v    <- as.numeric(cumsum(table(sapply(strsplit(hc_cells, '|', fixed = TRUE), function(x) x[2])))) + 0.5
  label_at   <- hc_cells[steps_v - round(table(sapply(strsplit(hc_cells, '|', fixed = TRUE), function(x) x[2])) / 2)]
  label_meta <- sapply(strsplit(label_at, '|', fixed = TRUE), function(x) x[2])
  
  p <- 
    counts_ds_f_smoo %>%  
      ggplot(aes(barcode, id_unique, fill = rollmean)) + 
        geom_raster() +
        scale_fill_gradientn(na.value = "pink", colours = colorRampPalette(c('lightgrey', 'orange', 'red', 'darkred', 'purple4', 'black'))(256)) +
        geom_hline(yintercept = steps, lwd = 0.25) +
        geom_vline(xintercept = steps_v, lwd = 0.25) +
        #scale_x_discrete(breaks = label_at, labels = label_meta) +
        theme(legend.position = 'none',
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
  
  return(p)
}

plotVarMean <- function(scSet,
                        n_label = 15,
                        cap     = Inf,
                        type = c('gene', 'te'))
{
  if (class(scSet) == "scSet")
  {
    gstats <- scSet@gstats
  } else {
    gstats <- scSet
  }
  
  if (sum(dim(gstats)) > 0)
  {

    if (type == 'gene')
    {
      gstats     <- gstats[which(type == 'gene'), ]
      gstats_subs <- gstats[which(class == 'protein_coding' & feature), 
        ][which(nchar(id_unique) < 10), ]
    }
    if (type == 'te' | type == 'te_fam')
    {
      gstats    <- gstats[which(type %fin% c('te', 'te_fam')), ]
      gstats_subs <- gstats[which(type == 'te_fam' & feature), ]
      if (nrow(gstats_subs) == 0)
      {
        gstats_subs <- gstats[which(type == 'te' & feature), ]
      }
    }
    
    if (is.null(gstats$class))
    {
      gstats$class <- ''
    }
    
    # cap varmean
    gstats$varmean <- ifelse(gstats$varmean > cap, cap, gstats$varmean)
    
    n_feats <- nrow(gstats[which(feature), ])  
    
    gstats_subs <- gstats_subs[, above := (log2(varmean) - vm_trend) >= sort(log2(varmean) - vm_trend, decreasing = TRUE)[min(nrow(gstats_subs), n_label)]][which(above), ]
    
    p_varmean <-
      gstats %>% 
        ggplot(aes(log2(mean_ds), log2(varmean), size = feature, alpha = feature, shape = type, col = class)) + 
          scale_size_manual(values  = c(0.1, 1)) +
          scale_alpha_manual(values = c(0.1, 0.5)) +
          scale_shape_manual(values = c(20, 8)) +
          geom_point() +
          ggrepel::geom_text_repel(data = gstats_subs, size = 2, alpha = 1, aes(label = id_unique)) + 
          ggtitle(paste0(n_feats, ' features')) + 
          theme(plot.title = element_text(size = 10, margin = margin(b = -15)),
                legend.position = 'none')
       
    if (type == 'te') 
    {
      p_varmean <- p_varmean + scale_color_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3"))
    } else {
      p_varmean <- p_varmean + scale_color_viridis_d(direction = -1)
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
                         theme = gridExtra::ttheme_default(base_size = 6)
                         )
  
  return(p)
}
#' Plot peak calling QC
#' @export
plotPeaks <- function(scSet)
{
  p_enr_peaks <- plotPeakEnr(scSet, type = 'te', what = 'enrichment')
  p_enr_stats <- plotPeakEnr(scSet, type = 'te', what = 'summary')
  p_enr_fams  <- plotPeakEnr(scSet, type = 'te', what = 'family')
  p_umi_loss  <- plotUMIs(scSet, what = 'peaks')

  # p_load      <- plotLoad(scSet, what = 'full')
  p_te_cov    <- plotTECov(scSet)
  #p_te_cov_h  <- plotTECov(scSet, what = 'histogram')
  
  p1 <- plot_grid(p_enr_peaks, 
                  p_enr_stats,
                  p_umi_loss,
                  nrow = 1,
                  rel_widths = c(2,1,2),
                  labels = 'auto')
                  
  # p2 <- plot_grid(p_te_cov,
                  # p_te_cov_h,
                  # nrow = 1,
                  # rel_widths = c(1, 1))
                  
  p_final <- plot_grid(p1, p_te_cov, rel_heights = c(1, 2), ncol = 1)
                  
  # p3 <- plot_grid(p2,
                  # plot_grid(p_load, labels = 'f'),
                  # ncol = 1,
                  # rel_heights = c(3, 1))
  
  return(p_final)
}

plotPeakEnr <- function(scSet,
                      type = 'te',
                      what = c('enrichment', 'summary', 'family'))
{
  if (sum(type %in% c('gene', 'te')) != 1)      { stop ('Type must be "te" or "gene"') }
  if (sum(what %in% c('enrichment', 'summary', 'family')) != 1) { stop ('What must be "enrichment", "summary", or "family"') }
  
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
    
    if (type == 'te') { p_peaks <- p_peaks + scale_fill_manual(values = c('DNA' = "#E41A1C", 'LINE' = "#377EB8", 'LTR' = "#4DAF4A", 'SINE' = "#984EA3")) }
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