compAlignments <- function(tes, 
                           BSgenome,
                           max_gap = 0.95, 
                           n_jobs = 500,
                           outdir = NULL,     
                           overwrite = FALSE)
{
  if (is.null(tes$name)) { stop ('No name column found in TE data.frame') }
  if (is.null(outdir))      { stop ('Please specificy directory to save MSA files to') }
  
  # add sequence to TE intervals
  tes <- addSeq(BSgenome, tes)
  
  # calc mafft settings based on estimated alignment complexity
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

  write_fst(res, paste0(outdir, BSgenome@provider_version, '_conv_df.fst'))
  
  # merge conversion df and te intervals
  #tes_dt      <- as.data.table(tes %>% select(-seqs))
  #tes_aligned <- merge(tes_dt, res_comb, by = 'id_unique', all.x = TRUE)
  
  # make GRanges
  #tes_aligned_gr <- as_granges(tes_aligned)
  
  # add te id as list name
  #names(tes_aligned_gr$conversion) <- tes_aligned_gr$id_unique
  
  return(res)
}

#' Counts number of reads per interval and cell barcode
compCounts <- function(reads, 
                       tes,
                       genes,
                       bin_size = 20,
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
  hits   <- matchReads(intervals_filt, reads, sense_only = sens_only)
  
  message ('Mapping read position on consensus') # check 0/1 based accuracy!
  if (!is.null(conversion)) 
  {
    # with conversion df
    pos_conv_df <- rbindlist(intervals_filt$conversion, idcol = 'id_unique')
    hits        <- merge(hits, pos_conv_df, all.x = TRUE, by = c('id_unique', 'pos'))
  } else {
    # using rmsk alignment (crude)
    hits[, ':=' (dist_start = pos, dist_end = width - pos)]
    hits[which(dist_start <= dist_end), pos_con := as.numeric(position_start + dist_start)]
    hits[which(dist_start >  dist_end), pos_con := as.numeric(position_end   - dist_end)]
    hits[which(pos_con > position_end | pos_con < position_start), pos_con := NA] 
  }
    
  # comp copy-number of consensus pos
  #cpn_df <- cpn(intervals)
  
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
    
	message ('Creating single-cell count matrix')
  # count how often interval occurs for given barcode (alternative might be to count unique and multi seperately)
  counts_long <- hits[,.(n = sum(NH_weight)), by = c('barcode', 'id_unique', 'NH_flag', 'pos_con', 'name', 'type', 'sense')]
  counts_wide <- dcast.data.table(counts_long, barcode + id_unique + name + type + pos_con + sense ~ NH_flag, value.var = 'n', fill = 0, sep = '_')
  res <- counts_wide[, .(id_unique, barcode, n = unique, n_all = unique + multi, pos_con, name, type, sense)]

  # add intervals to counts
  #intv_counts <- counts_wide[as.data.table(intervals_filt), nomatch = 0, on = 'id_unique']
  
  # add cpn and consize to counts
  #res <- merge(intv_counts, cpn_df, all.x = TRUE)

	return(res)
}

callCells <- function(counts,
                      method = 'emptyDrops',
                      lower  = 100,
                      niters = 10000,
                      fdr    = 0.01)
{
  smat <- dbutils::tidyToSparse(counts[which(type == 'gene'), c('name', 'barcode', 'n')])
  
  stats <- 
    DropletUtils::emptyDrops(smat, lower = lower, niters = niters, test.ambient=FALSE,
      ignore=NULL, alpha=NULL, BPPARAM=BiocParallel::SerialParam())

  counts_filt <- counts[which(barcode %in% rownames(stats)[which(stats$FDR < fdr)]), ]
  
  return(counts_filt)
}                      

callPeaks <- function(counts, 
                      tes,
                      background_correction = TRUE,
                      bin_size = 20,
                      smoothed = FALSE, 
                      min_umis = 10,
                      min_loci = 1,
                      enrichment = 1.5,
                      reg = 1)
{
  # get TE copynumber
  cpn_df <- cpn(tes)
  
  # aggregate counts across cells and TE loci
  counts_aggr <- aggr(counts)
  
  dens_df <- 
    merge(counts_aggr, cpn_df, all = TRUE)[which(is.na(bps)), bps := 20
      ][which(is.na(n_tot)), n_tot := 0
        ][, dens := n_tot / bps]
  
  
  
  
  
  peaks = dens_df %>% 
		group_by(name) %>%
			mutate(a = lag(dens, n = 5, default = NA),
				   b = lag(dens, n = 4, default = NA), 
				   c = lag(dens, n = 3, default = NA),
				   d = lag(dens, n = 2, default = NA),
				   e = lead(dens, n = 2, default = NA),
				   f = lead(dens, n = 3, default = NA),
				   g = lead(dens, n = 4, default = NA),
				   h = lead(dens, n = 5, default = NA)) %>% 
		    ungroup %>%
			select(name, bin_id, a,b,c,d,e,f,g,h) %>% 
			gather(offset, 'norm_umis', a,b,c,d,e,f,g,h) %>%
			group_by(name, bin_id) %>% 
				summarise(neighbouring_umis_mean = mean(norm_umis, na.rm = TRUE)) %>% 
			ungroup %>%
			left_join(bin_stats_norm, .) %>%
			mutate(fc = norm_umis / neighbouring_umis_mean)
			
	peak_sig = peaks %>% rowwise %>% mutate(p.value = ppois(q = norm_umis, lambda = neighbouring_umis_mean, lower.tail = FALSE, log = FALSE)) %>%
		ungroup %>% mutate(q.value = p.adjust(p.value, 'fdr'))
										  
	tss = peak_sig %>% mutate(tss = umis >= 20 & fc > 3 & q.value < 0.01)
  
  
  
  
  if (!is.null(counts$peak)) { stop ('Counts already contain peaks') }
  if (exclude_exons) { counts <- counts[which(!exonic), ] }
  
  dens_fg <- aggr(counts, sorted = TRUE)
  dens_bg <- data.table(name = tes$name, dens_bg = tes$dens_bg)[, .(dens_bg = mean(dens_bg, na.rm = TRUE)), by = 'name']
  
  dens_df <- merge(dens_fg, dens_bg, by = 'name')
  
  # calc stats
  stats <- # calcs rollstats and accounts for left/right boundaries
    dens_df[which(!is.na(pos_con)), 
            c('dens_bin_left', 'dens_bin_center', 'dens_bin_right', 'n_bin') := list(frollmean(dens, n = bin_size, align = 'left', algo = 'exact'), 
                                                                                     frollmean(dens, n = bin_size, align = 'center', algo = 'exact'),
                                                                                     frollmean(dens, n = bin_size, align = 'right', algo = 'exact'),
                                                                                     zoo::rollsum(n_tot, k = bin_size, fill = 'extend')), 
            by = 'name'][, dens_bin := ifelse(is.na(dens_bin_center), na.omit(dens_bin_left, dens_bin_right), dens_bin_center)]
  stats <- stats[, fc := dens_bin / min(dens_bin), by = 'name'] # + min(dens_bin)
  
  peaks <- stats[which(fc >= enrichment & n_bin >= min_umis & n_loci >= min_loci), .(peak = TRUE), by = c('name', 'pos_con')]
  
  res <- merge(counts, peaks, all.x = TRUE)[which(is.na(peak)), peak := FALSE]
  
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

normCounts <- function(counts,
                       N = 20000,
                       separate = FALSE,
                       long = TRUE)
{
  coutns <- as.data.table(counts)
  
  if (separate) 
  {
    # divide n by total number of reads per cell
    counts_frac <- counts[, n_norm := n / sum(n), by = 'barcode']
  } else {
    counts_frac <- counts[, n_norm := n / sum(n), by = c('barcode', 'type')]
  }
  
  # get integer counts based on weighted sampling
  counts_sampled <- 
    counts_frac[, .(id_unique = sample(id_unique, prob = n_norm, replace = TRUE, size = N)), by = 'barcode'
      ][, .(n_norm = .N), by = c('id_unique', 'barcode')]
  
  # add to input counts
  if (!long)
  {
    smat <- dbutils::tidyToSparse(counts_sampled)
    return(smat)
  } else {
    return(counts_sampled)
  }
}                      