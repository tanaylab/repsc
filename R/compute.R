compAlignments <- function(tes, 
                           BSgenome    = NULL,
                           max_gap     = 0.95, 
                           n_jobs      = 500,
                           outdir      = NULL,
                           settings_df = NULL,
                           overwrite   = FALSE)
{
  if (is.null(tes$name))    { stop ('No name column found in TE data.frame') }
  if (is.null(outdir))      { stop ('Please specificy directory to save MSA files to') }
  
  # add sequence to TE intervals
  tes <- addSeq(BSgenome, tes)
  
  # calc mafft settings based on estimated alignment complexity
  if (is.null(settings_df))
  {
    message ('Calculating settings for multiple sequence alignments')
    settings_df <- 
      data.table(bps = width(tes$seqs), 
                 name = tes$name)[, .(score = max(bps) / 2 * .N), by = 'name'] %>%
        mutate(perc_r     = percent_rank(score),
               group      = cut(perc_r, breaks = c(0, 0.5, 0.75, 0.95, 1.01), include.lowest = TRUE, labels = FALSE),
               threads    = c(1, 1, 1, 5)[group],
               fft        = c(TRUE, TRUE, TRUE, FALSE)[group],
               parttree   = c(FALSE, FALSE, TRUE, TRUE)[group],
               local      = c(TRUE, FALSE, FALSE, FALSE)[group],
               maxiterate = c(1000, 1000, 0, 0)[group],
               retree     = c(2, 2, 1, 0)[group]) %>%
        arrange(-perc_r)
  } 
   
  # commands for mafft alignment and conversion
  cmds <- list()
  for (fam in settings_df$name)
  {
    s <- settings_df %>% filter(name == fam)
    
    cmds[[fam]] <- 
      glue("Repexpr:::mapConMSA(Repexpr:::mafft(tes[which(tes$name == '{fam}'), ]$seqs, 
                                ids = tes[tes$name == '{fam}']$id_unique, 
                                overwrite = {overwrite}, 
                                retree = {s$retree}, 
                                maxiterate = {s$maxiterate}, 
                                parttree = {s$parttree}, 
                                save = TRUE, 
                                local = {s$local}, 
                                outdir = '{outdir}', 
                                auto = FALSE, 
                                long = TRUE, 
                                fft = {s$fft}, 
                                family = '{fam}', 
                                threads = {s$threads}))")
  }

  # create new empty env and fill with relevant
  empty <- new.env(parent = emptyenv())
  empty$tes <- tes
  
  # distribute, compute and combine
  res_gclust <- 
    gcluster.run3(command_list = cmds,  
                  packages = c("data.table", "plyranges", "base", "stats"), 
                  job_names = names(cmds),
                  max.jobs = n_jobs, 
                  envir = empty, 
                  io_saturation = FALSE)
  
  res <- rbindlist(lapply(res_gclust, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL)

  # export results
  write_fst(res, paste0(outdir, BSgenome@provider_version, '_conv_df.fst'))
  write.cvs(settings_df, paste0(outdir, BSgenome@provider_version, 'mafft_settings.csv'), quote = FALSE)
  
  return(res)
}

#' Counts number of reads per interval and cell barcode
compCounts <- function(reads, 
                       tes,
                       genes,
                       bin_size = 25,
                       sense_only = FALSE,
                       conversion = NULL,
                       protocol = 'fiveprime')
{
	if (class(reads) != 'GRanges') { reads <- as_granges(reads) }
  if (class(tes)   != 'GRanges') { tes   <- as_granges(tes) }
  if (class(genes) != 'GRanges') { genes <- as_granges(genes) }
  
  # combine te and gene intervals
  intervals <- 
    c(granges(tes) %>% 
        mutate(id_unique      = tes$id_unique, 
               name           = tes$name,
               position_start = tes$position_start,
               position_end   = tes$position_end,
               type = 'te'), 
      granges(genes) %>% 
        mutate(id_unique      = genes$id_unique, 
               name          = genes$gene_name,
               position_start = genes$position_start,
               position_end   = genes$position_end,
               type = 'gene'))
        
  GenomeInfoDb::seqlevelsStyle(intervals) <- GenomeInfoDb::seqlevelsStyle(reads)
  
  # exclude intervals not overlapping reads
  intervals_filt <- subsetByOverlaps(intervals, reads)
  
  # join reads by intervals
  hits   <- matchReads(intervals_filt, reads, sense_only = sense_only)
  
  message ('Mapping read position on consensus') # check 0/1 based accuracy!
  if (is.null(conversion)) 
  {
    # using rmsk alignment (crude)
    hits[, ':=' (dist_start = pos, dist_end = width - pos)]
    hits[which(dist_start <= dist_end), pos_con := as.numeric(position_start + dist_start)]
    hits[which(dist_start >  dist_end), pos_con := as.numeric(position_end   - dist_end)]
    hits[which(pos_con > position_end | pos_con < position_start), pos_con := NA] 
  } else {
    # with conversion df
    hits_gene <- hits[which(type == 'gene'), ][, pos_con := pos]
    hits_te   <- hits[which(type == 'te'), ]
    hits_te   <- merge(hits_te, conversion, all.x = TRUE, by = c('id_unique', 'pos', 'name'))
    hits      <- rbind(hits_gene, hits_te)
  }
  
  # add reads downstream TE if 3' protocol
  if (protocol == 'threeprime')
  {
    message ("Extending TE intervals at 3'")
    tes_3p  <- flank3prime(tes, genes, extend = 500)
    
    hits_3p <- setnames(matchReads(tes_3p, reads), 'pos', 'pos_con')[, pos_con := pos_con * -1]
        
    #cpn_3p  <- cpn(tes_3p)[, pos_con := pos_con * -1][which(!is.na(pos_con)), ]
    
    # combine
    hits   <- bind_rows(hits, hits_3p)
    #cpn_df <- bind_rows(cpn_3p, cpn_df)[, con_size := sum(unique(con_size)), by = 'name']
  }
  
  if (bin_size > 1)
  {
    message ("Binning positions")
    hits <- hits[, pos_con := cut(pos_con, breaks = seq(-1e6, 1e6, by = bin_size), include.lowest = TRUE)]
  }
    
	message ('Computing single-cell counts')
  # count how often interval occurs for given barcode (alternative might be to count unique and multi seperately)
  counts_long <- hits[,.(n = sum(NH_weight)), by = c('barcode', 'id_unique', 'NH_flag', 'pos_con', 'name', 'type', 'sense')]
  counts_wide <- dcast.data.table(counts_long, barcode + id_unique + name + type + pos_con + sense ~ NH_flag, value.var = 'n', fill = 0, sep = '_')
  res <- counts_wide[, .(id_unique, barcode, n = unique, n_all = unique + multi, pos_con, name, type, sense)]

	return(res)
}

callCells <- function(scSet,
                      method = 'emptyDrops',
                      lower  = 100,
                      niters = 10000,
                      fdr    = 0.01)
{
  counts      <- scSet@counts
  counts_gene <- counts[which(type == 'gene'), ]
  
  smat <- dbutils::tidyToSparse(counts_gene[, c('name', 'barcode', 'n')])
  
  stats <- 
    DropletUtils::emptyDrops(smat, lower = lower, niters = niters, test.ambient=FALSE,
      ignore=NULL, alpha=NULL, BPPARAM = BiocParallel::SerialParam())

  scSet@counts <- counts[which(barcode %in% rownames(stats)[which(stats$FDR < fdr)]), ]
  
  return(scSet)
}                      

callPeaks <- function(scSet, 
                      background_correction = TRUE,
                      neighb_bins = 5,
                      offset = 2,
                      min_umis = 10,
                      min_loci = 1,
                      enrichment = 1.5,
                      summarize = TRUE)
{
  bin_size <- scSet@bin_size
  counts   <- scSet@counts
  cpn_df   <- scSet@cpn
  #tes    <- scSet@tes
  
  # aggregate counts across cells and TE loci
  counts_aggr <- aggr(counts)
  
  # calc density
  dens_df <- 
    merge(counts_aggr, cpn_df, all = TRUE, by = c('name', 'pos_con'))[which(is.na(bps)), bps := bin_size
      ][which(is.na(n_tot)), n_tot := 0
        ][, dens := n_tot / bps
          ][which(!is.na(pos_con)),
            ][order(name, pos_con), ]
  
  # create neighbouring bin vector
  offset_v <- c((-neighb_bins - offset):-offset, offset:(neighb_bins + offset))
  
  # calc neighbouring signal
  n_tot_cols <- paste0('n_tot', offset_v)
  bps_cols   <- paste0('bps', offset_v)
  neighb_signal <- 
    dens_df[, c(n_tot_cols, bps_cols) := c(shift(n_tot, n = offset_v, type = 'lead'), shift(bps, n = offset_v, type = 'lead')), by = 'name']
   
  # calc neighbouring dens
  dens_df$dens_neighb <- 
    rowSums(neighb_signal[, n_tot_cols, with = FALSE], na.rm = TRUE) /
    rowSums(neighb_signal[, bps_cols, with = FALSE], na.rm = TRUE)
	
  # clean df
  dens_df <- dens_df[, !c(n_tot_cols, bps_cols), with = FALSE]
  
  # calc fc
  dens_df[, fc := dens / dens_neighb]
  
  # plot
  dens_df %>% mutate(fc = dens / dens_neighb) %>%  ggplot(aes(dens, fc, col = type)) + geom_point(size = 0.2, alpha = 0.25) + scale_x_log10() + ylim(0, 100)
  
  # add peak columns
  counts_peaks <- 
    merge(counts, 
          dens_df[, peak := n_tot >= 10 & fc >= 10][which(peak), c('name', 'pos_con', 'peak')], 
          all.x = TRUE, 
          by = c('name', 'pos_con'))[which(is.na(peak)), peak := FALSE]
  
  if (summarize)
  {
    res <- counts_peaks[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'id_unique', 'peak', 'sense', 'type')]
  }
  return(res)
}

compTELoad <- function(counts)
{
  gene_te_tot <- counts[, .(tot = sum(n)), by = c('barcode', 'type')]
  
  res <- 
    gene_te_tot %>%
    spread(type, tot) %>%
    mutate(te_load = te / (te + gene) * 100) %>%
    as.data.table()

  return(res)
}

normCounts <- function(scSet,
                       N = 20000,
                       separate = FALSE)
{
  counts <- scSet@counts
  seed   <- scSet@seed
  set.seed(seed)
  
  if (separate) 
  {
    counts_frac <- counts[, n_prob := n / sum(n), by = c('barcode', 'type')]
  } else {
    # divide n by total number of reads per cell
    counts_frac <- counts[, n_prob := n / sum(n), by = 'barcode']
  }
  
  # get integer counts based on weighted sampling
  counts_sampled <- 
    counts_frac[, .(id_unique = sample(id_unique, prob = n_prob, replace = TRUE, size = N)), by = 'barcode'
      ][, .(n_norm = .N), by = c('id_unique', 'barcode')]
  
  # add to input counts
  # if (!long)
  # {
    # res <- dbutils::tidyToSparse(counts_sampled)
    
  # } else {
    scSet@counts <- merge(counts, counts_sampled, all.x = TRUE, by = c('id_unique', 'barcode'))[which(is.na(n_norm)), n_norm := 0]
  # }
  
  return(scSet)
}                      