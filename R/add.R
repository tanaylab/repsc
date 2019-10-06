addAlignments <- function(scSet,
                          max_gap     = 0.95,
                          outdir      = NULL,
                          overwrite   = FALSE,
                          settings_df = NULL,
                          n_chunks    = FALSE,
                          msa_dir     = NULL)
{
  if (is.null(msa_dir)) { stop ('Please specify directory with sufficient space to store MSA output!') }
  scSet@msa_dir <- msa_dir
  tes           <- scSet@tes
  genome        <- scSet@genome
  
  compAlignments(tes, 
                 BSgenome    = genome, 
                 max_gap     = max_gap, 
                 n_jobs      = n_chunks,
                 outdir      = msa_dir,     
                 overwrite   = overwrite,
                 settings_df = NULL)
}

#' Add TE/Gene counts to scSet
#' @export
addCounts <- function(scSet,
                      bams            = NULL,
                      msa_dir         = NULL,
                      te_exon_overlap = 'filter',
                      hierarchy       = 'locus',
                      bin_size        = NA,
                      min_cell_size   = 100,
                      use_gcluster    = FALSE)
{
  bams           <- as.data.table(bams)
  bin_size       <- as.integer(bin_size)
  scSet@params[['addCounts']]$bin_size        <- bin_size
  scSet@params[['addCounts']]$te_exon_overlap <- te_exon_overlap
  scSet@input    <- bams
  scSet@bin_size <- bin_size
  tes            <- scSet@tes
  tes_3p         <- scSet@tes_3p
  genes          <- scSet@genes
  protocol       <- scSet@protocol
  
  if (te_exon_overlap == 'filter')
  {
    tes       <- tes %>% filter(is.na(gene))
    scSet@tes <- tes

    if (length(tes_3p > 0))
    {
      tes_3p       <- tes_3p %>% filter(!id_unique %in% tes$id_unique)
      scSet@tes_3p <- tes_3p
    }
  }
  
  # check input
  checkInput(bams)
  
  # convert factor to character
  bams$paths <- as.character(bams$paths)
  
  # overwrite downstream results
  scSet@counts    <- data.frame()
  scSet@counts_ds <- data.frame()
  scSet@cstats    <- data.frame()
  scSet@pstats    <- data.frame()
  scSet@gstats    <- data.frame()
  scSet@plots     <- list()
  scSet@cells     <- character()
  
  if (!is.null(msa_dir))
  {
    message ('Importing conversion data.frame')
    conv_df <- fst::read_fst(paste0(msa_dir, '/', scSet@genome@provider_version, '_conv_df.fst'), as.data.table = TRUE)
    message ('Done')
  } else {
    conv_df <- NA
  }
  
  if (is.null(bams$chunk))
  {
    bams$chunk <- 1
  }
    
  #' Count number of reads per interval and cell barcode
  if (use_gcluster)
  {
    cmds <- list()
    for (i in unique(bams$chunk))
    {
      cmds[[as.character(i)]] <- 
        glue("Repsc:::compCounts(bams = bams[bams$chunk == {i}, ],
                                 tes,
                                 tes_3p,
                                 genes,
                                 hierarchy    = '{hierarchy}',
                                 bin_size     = {bin_size},
                                 conversion   = {conv_df},
                                 protocol     = '{protocol}',
                                 use_gcluster = FALSE)")
    }
    # create new empty env and fill with relevant
    empty        <- new.env(parent = emptyenv())
    empty$bams   <- bams
    empty$tes    <- tes
    empty$tes_3p <- tes_3p
    empty$genes  <- genes
    
    # distribute, compute and combine
    res <- 
      Reputils:::gcluster.run3(command_list = cmds,  
                    packages = c("data.table", "plyranges", "base", "stats"), 
                    max.jobs = 50, 
                    envir = empty, 
                    io_saturation = FALSE)
    
    message ('Combining datasets')
    res <- rbindlist(lapply(res, function(x) x$retv))
    message ('Done')
  } else {
    res <- 
      compCounts(bams, 
                 tes,
                 tes_3p,
                 genes,
                 hierarchy   = hierarchy,
                 bin_size    = bin_size,
                 sense_only  = TRUE,
                 conversion  = conv_df,
                 protocol    = protocol)
  }
  # summarize results in case of dumb chunking
  if (hierarchy == 'locus')
  {
    if (!is.na(bin_size))
    {
      counts <- res[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'name', 'id_unique', 'pos_con', 'type', 'meta')]
    } else {
      counts <- res[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'name', 'id_unique', 'type', 'meta')]
      counts <- counts[, pos_con := NA]
    }
  } else {
    if (!is.na(bin_size))
    {
      counts <- res[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'name', 'pos_con', 'type', 'meta')]
    } else {
      counts <- res[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'name', 'type', 'meta')]
      counts <- counts[, pos_con := NA]
    }
    counts <- counts[, id_unique := name]
  }
  
  # import 10x gene counts if defined
  if (!is.null(bams$hdf5))
  {
    message ('Importing 10x gene counts from hdf5 file(s)')
    # make hdf5 files unique
    bams <- unique(bams[, intersect(c('hdf5', 'meta'), colnames(bams)), with = FALSE])
    
    # import counts as list and name according to meta
    counts_10x_l <- lapply(bams$hdf5, import10x)
    names(counts_10x_l) <- bams$meta
    
    # unlist and format
    counts_10x_dt <-
      rbindlist(counts_10x_l, idcol = 'meta')[, ':=' (pos_con = NA, id_unique = name, type = 'gene', n_all = NA, barcode = paste(barcode, meta, sep = '|'))]
  
    # combine with TE counts
    counts <- rbindlist(list(counts, counts_10x_dt), use.names = TRUE)
  
    # mark valid barcodes based on 10x
    valid_cells <- unique(counts_10x_dt$barcode)
    counts      <- counts[, is_cell := barcode %in% valid_cells][which(is_cell), ]
    scSet@cells <- valid_cells  
  }
  
  # initially throw small cells based on unique genic reads/UMIs
  bcs_keep <- counts[which(type == 'gene'), .(csize = sum(n)), by = 'barcode'][which(csize >= min_cell_size), ]$barcode
  counts_f <- counts[which(barcode %fin% bcs_keep), ]
  
  scSet@counts   <- counts_f
  
  # if (te_exon_overlap == 'resolve')
  # {
    # scSet <- resolveOverlaps(scSet, min_expr = 100)
  # }
    
  # add class annotation
  scSet <- annoClass(scSet)
  
  scSet@cpn      <- cpn(scSet)
  scSet@cstats   <- cstats(scSet) # needs class annotation!
  scSet@mstats   <- mstats(scSet) # needs class annotation!
    
  return(scSet)
}