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
compCounts <- function(bams, 
                       tes,
                       tes_3p,
                       genes,
                       hierarchy = 'locus',
                       bin_size = NA,
                       sense_only = TRUE,
                       conversion = NULL,
                       protocol = 'fiveprime',
                       use_gcluster = FALSE)
{
  if (class(tes)   != 'GRanges') { tes   <- as_granges(tes) }
  if (class(genes) != 'GRanges') { genes <- as_granges(genes) }
  
  # slim metadata columns
  tes   <- tes %>% select(id_unique, name, position_start, position_end)
  genes <- genes %>% select(id_unique, name, position_start, position_end)
  
  # import reads from bam file
  reads <- Reputils::importBAM(bams, anchor = protocol, use_gcluster = use_gcluster, tag = 'NH')
  reads <- reads[which(!is.na(barcode)), ]
  
  # remove duplicated read entries (can happen if neglected mate (2nd) maps to different positions). NH doesn't need to be adjusted
  reads <- unique(reads)
  
  # add dummy meta column if none provided
  if (is.null(reads$meta))
  {
    reads <- reads[, ':=' (meta = 'sample', barcode = paste(barcode, 'sample', sep = '|'))]
  }
  
  # combine TE intervals downstream if threeprime protocol is used
  if (protocol == 'threeprime')
  {
    tes <- c(tes_3p, tes)
  }
  
  # combine te and gene intervals if no hdf5 is provided (will also quantify genic counts)
  if (is.null(bams$hdf5))
  {
    intervals <- 
      c(mutate(tes, type = 'te'), 
        mutate(genes, type = 'gene'))
  } else {
    intervals <- mutate(tes, type = 'te')
  }
        
  # convert to data.table
  intervals <- as.data.table(intervals)
  
  # join reads by intervals
  hits   <- matchReads(intervals, reads, sense_only = sense_only)
  
  # split genic and te counts
  hits_gene <- hits[which(type == 'gene'), ][, ':=' (pos_con = 1L, id_unique = name)] # removes information about position on gene
  hits_te   <- hits[which(type == 'te'), ]
  
  message ('Mapping read position on consensus') # check 0/1 based accuracy!
  if (is.na(conversion)) 
  {
    # using rmsk alignment (crude)
    hits_te[, ':=' (dist_start = start_read - start, dist_end = end - start_read)]
    hits_te[which(dist_start <= dist_end & strand == '+'),  pos_con := position_start + dist_start]
    hits_te[which(dist_start > dist_end & strand == '+'),  pos_con := position_end   - dist_end]
    hits_te[which(dist_start <= dist_end & strand == '-'),  pos_con := position_start - dist_start]
    hits_te[which(dist_start > dist_end & strand == '-'),  pos_con := position_end   + dist_end]
    hits_te[which(strand == '+' & (pos_con > position_end | pos_con < position_start)), pos_con := NA
         ][which(strand == '-' & (pos_con > position_start | pos_con < position_end)), pos_con := NA]
  } else {
    # with conversion df

    hits_te   <- merge(hits_te, conversion, all.x = TRUE, by = c('id_unique', 'pos', 'name'))
  }
  
  # combine te and genic
  hits      <- rbindlist(list(hits_gene, hits_te[, !c('dist_start', 'dist_end')]))
  
  # modify downstream extended counts by multiplying * -1
  if (protocol == 'threeprime')
  {
    hits <- hits[which(downstream), pos_con := as.integer(pos_con * -1)]
  }
  
  if (!is.na(bin_size))
  {
    message ("Binning positions")
    hits <- hits[, pos_con := cut(pos_con, breaks = seq(-1e6, 1e6, by = bin_size), include.lowest = TRUE)] # use labels = FALSE might accelerate speed
  } else {
    hits[, pos_con := 1L]
  }
  
  if (protocol == 'threeprime' & !is.na(bin_size))
  {
    # reorder factor levels so that downstream extension is at the end
    hits <- hits[, pos_con := forcats::fct_relevel(pos_con, rev(grep('-', levels(pos_con), fixed = TRUE, value = TRUE)), after = Inf)]
  }
  
	message ('Computing single-cell counts')
  # count how often interval occurs for given barcode 
  if (hierarchy == 'locus')
  {
    res <- hits[, ':=' (uniq = NH == 1, NH_weight = 1 / NH)][,.(n = sum(uniq), n_all = sum(NH_weight)), by = c('barcode', 'id_unique', 'pos_con', 'name', 'type', 'meta')]
  } else {
    res <- hits[, ':=' (uniq = NH == 1, NH_weight = 1 / NH)][,.(n = sum(uniq), n_all = sum(NH_weight)), by = c('barcode', 'pos_con', 'name', 'type', 'meta')]
    res <- res[, id_unique := name]
  }

  # round n_all
  res <- res[, n_all := round(n_all, 3)]
  
  # remove pos_con info from genes and tes if not defined
  if (is.na(bin_size))
  {
    res <- res[, pos_con := NA]
  } else {
    res <- res[which(type == 'gene'), pos_con := NA]
  }
  
	return(res)
}

cstats <- function(scSet) # slow, can be improved!
{
  if (sum(dim(scSet@cstats)) == 0)
  {
    message ('Computing cell stats')
    ribo_prefixes <- c('Human' = '^RP[LS]', 'Mouse' = '^Rp[ls]')
    mito_prefixes <- c('Human' = '^MT-', 'Mouse' = '^mt-')   
    genes         <- scSet@genes
    organism      <- scSet@genome@common_name
    counts        <- scSet@counts
    counts_gene   <- counts[which(type == 'gene'), ]
    counts_te     <- counts[which(type == 'te'), ]
    ribo_prefix   <- ribo_prefixes[organism]
    mito_prefix   <- mito_prefixes[organism]
    
    mt_genes <- 
      genes %>%
        filter(grepl(mito_prefix, name)) %>%
        as_tibble %>%
        pull(name) %>%
        unique
        
    ribo_genes <-
      genes %>%
        filter(grepl(ribo_prefix, name)) %>%
        as_tibble %>%
        pull(name) %>%
        unique
     
    # tag mito/ribo genes
    counts_gene <- 
      counts_gene[, .(gene_tot = sum(n), gene_tot_all = sum(n_all)), by = c('barcode', 'name', 'meta')][, ':=' (mito = name %fin% mt_genes, ribo = name %fin% ribo_genes)]
    
    # summarize percentage  
    summary_genes <-
      counts_gene[, .(size_genic = sum(gene_tot), 
                      size_genic_all = sum(gene_tot_all),
                      n_genes    = length(unique(name)),
                      n_ribo = sum(gene_tot[ribo]), 
                      n_mito = sum(gene_tot[mito])), by = c('barcode', 'meta')
        ][, ':=' (perc_mito = round(n_mito / size_genic * 100, 2),
                  perc_ribo = round(n_ribo / size_genic * 100, 2))]
    
    # summarize umis per class and mapping (all, uniq)
    summary_tes <-
      counts_te[, .(all = sum(n_all), uniq = sum(n)), by = c('barcode', 'class')] %>%
        tidyr::gather('mapping', 'n', -1:-2) %>%
        mutate(class = paste(class, mapping, sep = '_')) %>% 
        select(-mapping) %>% 
        tidyr::spread(class, n, fill = 0) %>% 
        as.data.table()

    # calc total TE umis
    patterns  <- c('all', 'uniq')
    TE_counts <- sapply(patterns, function(xx) rowSums(summary_tes[,grep(xx, names(summary_tes)), drop=FALSE, with = FALSE]))
    colnames(TE_counts) <- paste0('TE_', patterns)
    summary_tes <- cbind(summary_tes, TE_counts)
    
    # cbind genic and te stats
    bcstats <- 
      merge(summary_genes, 
            summary_tes, 
            by = c('barcode'), 
            all.x = TRUE)[, lapply(.SD, function(x) { ifelse(is.na(x), 0, x) }) # replaces any TE NAs with 0 
        ][, !c('n_ribo', 'n_mito')]
    
    # rank barcodes
    bcstats <- bcstats[order(meta, -size_genic), rank := 1:.N, by = 'meta']
  } else {
    bcstats <- scSet@cstats
  }  
  
  return(bcstats) 
}

mstats <- function(scSet,
                   reg = 100)
{
  message ('Computing mapping stats')
  counts <- scSet@counts
  
  mstats <- counts[, .(uniq = sum(n), all = sum(n_all)), by = c('name', 'type', 'class')
    ][, ':=' (perc_unique = uniq / all * 100, delta = all - uniq, fc = (all + reg) / (uniq + reg))]
  
  return(mstats)
}
  
compClassLoad <- function(scSet = NULL)
{
  counts <- suppressWarnings(counts(annoClass(scSet)))

  # get counts per class, peak, barcode
  class_tot <- 
    counts[, .(tot_all = sum(n_all), tot_uniq = sum(n), tot_uniq_p = sum(n[peak])), by = c('barcode', 'class')] %>%
    tidyr::complete(., barcode, class, fill = list(tot_all = 0, tot_uniq = 0, tot_uniq_p = 0)) %>%
    as.data.table()
  
  # calc class load
  class_load <- 
    class_tot[, ':=' (load_all    = tot_all / sum(tot_all) * 100, 
                        load_uniq   = tot_uniq / sum(tot_uniq) * 100,
                        load_uniq_p = tot_uniq_p / sum(tot_uniq_p) * 100),
                        by = 'barcode']

  # add type info
  res <-
    merge(class_load, 
          unique(counts[, c('class', 'type')]),
          all.x = TRUE,
          by = 'class')
        
  return(res)
}

#' Normalize read counts
#' @export
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