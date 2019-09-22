resolveOverlaps <- function(scSet,
                            min_expr = 50)
{
  message ('Resolving TE/gene exon overlaps')
  tes       <- scSet@tes
  tes_exon  <- data.table(id_unique = tes$id_unique, fam = tes$name, exonic = tes$exonic, host_gene = tes$gene)[which(exonic), ]
  counts    <- scSet@counts
  cstats    <- scSet@cstats
  set.seed(scSet@seed)
  
  # select exonic TEs with decent expression (to do stats)
  loc_expr <- top_n(counts[, .(n = sum(n)), by = 'id_unique'][which(n >= min_expr), ][which(id_unique %fin% tes_exon$id_unique),], 1e4, n)$id_unique
  tes_exon_f <- tes_exon[which(id_unique %fin% loc_expr), ]
  
  # select barco des (with highest expression per locus)
  n_lim  <- round(1e8 / length(loc_expr))
  n_samp <- 10
  counts_te_ex <- counts[which(id_unique %fin% loc_expr), ]
  bcs    <- unique(counts_te_ex[, .(barcode = barcode[order(n, decreasing = TRUE)][1:n_samp]), by = 'id_unique'][which(!is.na(barcode)), ]$barcode)
  while (length(bcs) < n_lim)
  {
    bcs    <- unique(counts_te_ex[, .(barcode = barcode[order(n, decreasing = TRUE)][1:n_samp]), by = 'id_unique'][which(!is.na(barcode)), ]$barcode)
    n_samp <- n_samp + 10
  }
  
  # filter counts for relevant families/host genes
  genefams <- unique(c(tes_exon_f$host_gene, tes_exon_f$fam))
  counts_f <- counts[which(name %fin% genefams & barcode %fin% bcs), ]
    
  # aggregate per fam and gene (all)
  famgene_aggr <- counts_f[, .(n = sum(n)), by = c('barcode', 'name', 'type')]
  gene_aggr    <- famgene_aggr[, ':=' (host_gene = name, n_gene = n)][which(type == 'gene'), c('barcode', 'host_gene', 'n_gene')]
  fam_aggr     <- famgene_aggr[,   ':=' (fam = name, n_fam = n)][which(type == 'te'), c('barcode', 'fam', 'n_fam')]
  
  # aggre loci (only exonic)
  loci_aggr <- 
    counts_f[which(id_unique %fin% loc_expr), .(n_loc = sum(n)), by = c('barcode', 'id_unique')]
  
  # expand
  bc_univ <- 
    merge(tes_exon_f,
          loci_aggr[CJ(barcode = unique(loci_aggr$barcode), 
                       id_unique = unique(loci_aggr$id_unique)), 
                       on = c('barcode', 'id_unique')],
          by = 'id_unique',
          all = FALSE)
 
  # combine locus, host gene fam counts
  locfamgene <- 
    merge(
          merge(bc_univ, 
                fam_aggr, 
                all.x = TRUE, 
                by = c('barcode', 'fam')), 
          gene_aggr, 
          all.x = TRUE, 
          by = c('barcode', 'host_gene'))
  locfamgene[is.na(locfamgene)] <- 0
  

  # subtract locus counts and transform
  locfamgene[, ':=' (n_fam  = log2(1 + 7 * (n_fam - n_loc)), 
                     n_gene = log2(1 + 7 * pmax(0, n_gene - n_loc)), 
                     n_loc  = log2(1 + 7 * n_loc))] # negative gene count possible if TE partially maps outside exon and read maps there
  
  # comp cors
  cors <- locfamgene[, .(cor_fam = cor(n_loc, n_fam), cor_gene = cor(n_loc, n_gene)), by = 'id_unique']
  cors[is.na(cors)] <- 0 # NA in gene cor can happen if gene is all 0
  
  # rescue TEs
  tes_rescue <- cors[which(cor_fam >= 0.4 & cor_gene < 0.4),]$id_unique
  
  # correct genic counts
  test <- 
    merge(counts_f[which(id_unique %fin% tes_rescue), .(n_te = sum(n), n_all_te = sum(n_all)), by = c('barcode', 'id_unique')], tes_exon_f)[, .(n_te = sum(n_te), n_all_te = sum(n_all_te)), by = c('barcode', 'host_gene')][, name := host_gene][, !'host_gene']
  merge(counts, test, all.x = TRUE, by = c('barcode', 'name'))[which(is.na(id_unique)),]
  
  
  # label suspicious genes
  genes_label <- cors[which(cor_fam >= 0.4 & cor_gene >= 0.4), ]$id_unique
  genes_label <- sapply(strsplit(genes_label, '|', fixed = TRUE), function(x) x[3])
  
  # update counts
  counts_f <- counts[which(!id_unique %fin% tes_exon), ]
  
  scSet@counts <- counts_f
  return(scSet)
}


#' Call cell barcodes
#'
#' @param method Character. 'Emptydrops' (the default) or 'size' are currently supported.
#' @param lower Integer. Emptydrops parameter, see
#' @param niters Integer. Emptydrops parameter, see
#' @param min_size Integer. Size parameter, minimal number of molecules per cell.
#' @param fdr Numeric. Emptydrops parameter, see
#'
#' @export
selectCells <- function(scSet,
                        method   = 'size',
                        min_size = 0,
                        max_mito = 1,
                        min_ribo = 0,
                        lower    = 100,
                        niters   = 10000,
                        fdr      = 0.01)
{
  message ('Selecting valid barcodes')
  bcstats <- scSet@cstats
  counts <- scSet@counts
  
  # if (method == 'emptyDrops')
  # {
    # smat <- Reputils::longToSparse(counts_gene[, c('name', 'barcode', 'n')])
    
    # stats <- 
      # DropletUtils::emptyDrops(smat, lower = lower, niters = niters, test.ambient=FALSE,
        # ignore=NULL, alpha=NULL, BPPARAM = BiocParallel::SerialParam())
  
    # barcode_stats <-
      # counts_gene[, .(n_tot = sum(n)), by = 'barcode'
          # ][, is_cell := barcode %in% rownames(stats)[which(stats$FDR < fdr)]]
  
  # }
  
  # call cells
  bcstats <- 
    bcstats[, ':=' (pass_size = size_genic >= min_size,
                   pass_miri = perc_ribo >= 100 * min_ribo & perc_mito <= 100 * max_mito)
             ][, is_cell := (pass_size + pass_miri) == 2]
  
  true_cells <- unique(bcstats[which(is_cell),]$barcode)
  
  # update counts and cstats
  counts <- counts[which(barcode %fin% true_cells), ]
  bcstats <- bcstats[which(barcode %fin% true_cells), ]
  
  # update scSet
  scSet@cells   <- true_cells
  scSet@counts  <- counts
  scSet@cstats  <- bcstats
  
  message ('Retained ', length(true_cells), ' barcodes')

  return(scSet)
}  

#' Adds gene/te stats to scSet
#' @export
selectFeatures <- function(scSet,
                           ds = NULL,
                           class = c('DNA', 'SINE', 'LINE', 'LTR'),
                           peak_only = FALSE,
                           reg = 10,
                           min_expr = 25,
                           min_expr_third = 3,
                           min_norm_var = 0.3,
                           k = NULL)
{
  #params <- scSet@params
  counts     <- Repsc::counts(scSet)
  cstats     <- scSet@cstats
  params_old <- scSet@params
  params_new <- params_old
  params_new[['selectFeatures']]$ds <- ds
  seed      <- scSet@seed
  
  if (length(class) != 4)
  {
    te_classes = class
    counts <- counts[which(class %in% te_classes | type == 'gene'), ]
  }
  
  if (peak_only)
  {
    if (is.null(counts$peak))
    {
      stop ('No peaks found in scSet, run "selectPeaks" function first')
    }
    counts <- counts[which(peak | type == 'gene'),]
  }
  
  if(is.null(params_new[['selectFeatures']]$ds)) 
  {
    ds <- min(median(cstats$size_genic), max(750, round(quantile(cstats$size_genic, 0.05))))
    params_new[['selectFeatures']]$ds <- ds
  }
  
  if (is.null(params_old[['selectFeatures']]$ds))
  {
    downsamp <- TRUE
  } else {
    downsamp <- params_new[['selectFeatures']]$ds != params_old[['selectFeatures']]$ds
  }
  
  if (downsamp)
  {
    # aggregate TEs per id_unique and genes per name and barcode
    counts_aggr <- 
      counts[, .(n = sum(n)), by = c('id_unique', 'name', 'barcode', 'type', 'class')]
    
    # downsample
    counts_ds   <- Reputils::downsamp(counts_aggr, i = 'id_unique', full = FALSE, ds = ds)
    
    # append family ds counts (careful, per cell umis are not equal afterwards!)
    counts_fam_ds <- 
      counts_ds[, .(n_ds = sum(n_ds), n = sum(n)), by = c('name', 'type', 'barcode', 'id_unique', 'class')][which(type == 'te'), ][, type := 'te_fam']
    counts_ds_comb <- rbindlist(list(counts_ds, counts_fam_ds), use.names = TRUE)
    
    gstats <- varmean(counts_ds_comb, reg = reg, min_expr = min_expr)
  } else {
    gstats    <- scSet@gstats
    counts_ds <- scSet@counts_ds
    
    # recompute family ds counts based on counts_ds (for later feature selection)
    counts_fam_ds <- 
      counts_ds[, .(n_ds = sum(n_ds), n = sum(n)), by = c('name', 'type', 'barcode', 'id_unique', 'class')][which(type == 'te'), ][, type := 'te_fam']
    counts_ds_comb <- rbindlist(list(counts_ds, counts_fam_ds), use.names = TRUE)
  }
    
  # select features
  feats <- gstats[, .(feature = tot_ds >= min_expr & log2(varmean) - vm_trend >= min_norm_var)]$feature
  feats <- gstats[feats, ]$id_unique
  # further select features based on 3rd highest expression (only computed on already selected features to save computation)
  feats <- counts_ds_comb[id_unique %fin% feats, .(n_ds_third = sort(n_ds, decreasing = TRUE)[3]), by = 'id_unique'][which(n_ds_third >= min_expr_third),]$id_unique
  
  message (length(feats), ' features selected')
  
  # compute TE modules  
  if (!identical(feats, try(gstats[which(feature), ]$id_unique, silent = TRUE)))
  {
    # aggregate feature counts per barcode and id_unique
    counts_ds_f <- counts_ds[which(id_unique %fin% feats & type == 'te'), .(n_ds = sum(n_ds)), by = c('id_unique', 'barcode', 'class')]
    
    if (nrow(counts_ds_f) > 0)
    {
      smat     <- log2(1 + 7 * as.matrix(Reputils::longToSparse(counts_ds_f[, !'class'])))
      cor_mat  <- tgs_cor(t(smat))
      pca <- prcomp(cor_mat, scale = FALSE, center = FALSE)
     
      #dist_mat <- sqrt(2 * (1 - cor_mat))
      if (is.null(k)) { k <- ceiling(nrow(smat) / 100) }
      set.seed(seed)
      ki  <- clustLouvain(pca$rotation[,1:2], k = k)
      
      # cor_mat  <- tgs_cor(t(smat))
      # dist_mat <- sqrt(2 * (1 - cor_mat))
      # mods <- tapply(1:nrow(cor_mat), as.character(structure(counts_ds_f$class, names = counts_ds_f$id_unique)[rownames(cor_mat)]), simplify = FALSE, function(x)
      # {
        # return(cutree(hclust(as.dist(dist_mat[x, x]), 'ward.D2'), min(length(x), 3)))
      # })
      
      # for (i in 1:length(mods))
      # {
        # mods[[i]] <- structure(paste(names(mods)[i], mods[[i]], sep = '_'), names = names(mods[[i]]))
      # }
      # names(mods) <- NULL
      # ct <- unlist(mods)
            

      # ki <- clustLouvain(smat, k = k)
      # ct <- cutree(hclust(as.dist(dist_mat), 'ward.D2'), 8)
      #ct <- kmeans(dist_mat, 8)$cluster
      #ki <- data.table(id_unique = names(ct), module = ct)[order(module), ][, module :=  cumsum(module != shift(module, n=1L, type = 'lag', fill = FALSE))]
      
      # add module info as column
      gstats <- 
        merge(gstats[, colnames(gstats) != 'module', with = FALSE], 
              ki, all.x = TRUE, by = 'id_unique')
    } else {
      gstats <- gstats[, module := NA]
    }
  }
  
  gstats$feature <- gstats$id_unique %fin% feats  
    
  scSet@gstats    <- gstats
  scSet@ds        <- as.integer(ds)
  scSet@counts_ds <- counts_ds
  
  # add class anno
  scSet <- annoClass(scSet)
    
  return(scSet)
}

#' Selects enriched signal on consensus models
#' @export
selectPeaks <- function(scSet, 
                        neighb_bins = 10,
                        offset = 4,
                        min_expr = 10,
                        min_fc = 1,
                        fdr = 0.01,
                        min_loci = 3,
                        summarize = FALSE)
{
  if (!is.na(scSet@params$addCounts$bin_size))
  {
    set.seed(scSet@seed)
    message( 'Calling peaks')
    bin_size <- scSet@bin_size
    cstats   <- scSet@cstats
    counts   <- Repsc::counts(scSet)
    cpn_df   <- scSet@cpn
    genes    <- scSet@genes
    
    # throw pos_con NA counts
    counts <- counts[which(!is.na(pos_con)), ]
    
    # aggregate counts across cells and TE loci
    counts_aggr <- aggr(counts)
    
    # calc density
    dens_df <- 
      merge(cpn_df, counts_aggr, all.x = TRUE, by = c('name', 'pos_con'))[which(is.na(n_tot)), n_tot := 0
          ][, dens := n_tot / bps
            ][order(name, pos_con), ]
    
    # sample counts with dens as prob
    count_samp <- 
      dens_df[which(dens > 0), 
        ][, .(pos_con = sample(pos_con, prob = dens, replace = TRUE, size = sum(n_tot))), by = 'name'][, .(n_samp = .N), by = c('name', 'pos_con')]
    
    # merge sampled and original counts
    counts_obs_exp <- merge(dens_df, count_samp, all.x = TRUE, by = c('name', 'pos_con'))[which(is.na(n_samp)), n_samp := 0]
    
    # create neighbouring bin vector
    offset_v <- c((-neighb_bins - offset):-offset, offset:(neighb_bins + offset))
    
    # calc neighbouring signal
    n_tot_cols <- paste0('n_tot', offset_v)
    neighb_signal <- 
      counts_obs_exp[, c(n_tot_cols) := c(data.table::shift(n_samp, n = offset_v, type = 'lead')), by = 'name']
     
    # calc neighbouring signal mean
    counts_obs_exp$n_samp_neighb <- 
      rowMeans(neighb_signal[, n_tot_cols, with = FALSE], na.rm = TRUE)
    
    # clean df
    counts_obs_exp <- counts_obs_exp[which(n_samp > 0), !c(n_tot_cols), with = FALSE]
    
    # call enriched regions
    fdr_t <- fdr
    peaks <- 
      counts_obs_exp[, p_val := ppois(q = n_samp, lambda = n_samp_neighb, lower.tail = FALSE, log = FALSE)
        ][, fdr := p.adjust(p_val, 'bonferroni')
          ][, fc := log2(n_samp / n_samp_neighb)
            ][, peak := fdr <= fdr_t & n_tot >= min_expr & fc >=  min_fc
              ][which(is.na(peak)), peak := FALSE] # n_samp and hence peak can be NaN/NA if shift ends in neverland
     
    # add peak column
    counts_new <- 
      merge(scSet@counts[, !'peak'], 
            peaks[which(peak), c('name', 'pos_con', 'peak')], 
            all.x = TRUE, 
            by = c('name', 'pos_con'))[which(is.na(peak)), peak := FALSE]
    
    # update cstats with peak TE load
    cstats_peaks <-
      counts_new[which(type == 'te' & peak), .(n = sum(n)), by = c('barcode', 'class')
        ][, class := paste(class, 'peak', sep = '_')
          ] %>% 
        tidyr::spread(class, n, fill = 0) %>% 
        mutate(TE_peak = DNA_peak + LINE_peak + LTR_peak + SINE_peak) %>% 
        as.data.table()
    cstats_new <- 
      merge(cstats[, !c("DNA_peak", "LINE_peak", "LTR_peak", "SINE_peak", "TE_peak")],       # remove peak columns from previous call
            cstats_peaks, all = TRUE)[, lapply(.SD, function(x) { ifelse(is.na(x), 0, x) })] # replaces any TE NAs with 0
    
    if (summarize)
    {
      counts_new <- counts_new[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'id_unique', 'peak', 'sense', 'type')]
    }
    
    scSet@counts <- counts_new
    scSet@pstats <- peaks
    scSet@cstats <- cstats_new
    
    if (nrow(peaks[which(type == 'gene'),]) > 0)
    {
      scSet        <- translatePeakCoords(scSet)
    }
    
    return(scSet) 
  } else {
    message ('Define bin size when running "addCounts" function to allow peak calling')
  }
}

getTSS <- function(genes)
{
  transcripts <- getTranscripts(genes)
  
  tss <- mutate(anchor_5p(transcripts), width = 1)
  return(tss)
}

getTTS <- function(genes)
{
  transcripts <- getTranscripts(genes)
  
  tts <- mutate(anchor_5p(transcripts), width = 1)
  
  return(tts)
}

getExons <- function(genes,
                     transcripts = NULL)
{
  exons <- genes[genes$type == 'exon']
  
  if (!is.null(transcripts))
  {
    exons <- exons %>% filter(transcript_id %in% transcripts)
  }
  
  return(exons)
}

getTranscripts <- function(genes)
{
  transcripts <- genes[genes$type == 'transcript']
  return(transcripts)
}

getIntrons <- function(genes)
{
  # transc <- getTranscripts(genes)
  # exons  <- getExons(genes)
  
  # if (!is.null(transcripts))
  # {
    # transc <- transc %>% filter(transcript_id %in% transcripts)
    # exons  <- exons %>% filter(transcript_id %in% transcripts)
  # } 
  
  # introns     <- setdiff_ranges_directed(transc, exons)
  
  genes_dt   <- as.data.table(genes)
  gene_intvs <- as_granges(genes_dt[, .(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1]), by = 'id'])
  
  introns <- GenomicRanges::setdiff(gene_intvs, genes, ignore.strand = FALSE)
  
  return(introns)
}