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

#' Adds single-cell TE and genic counts 
#'
#' @param bams data.frame. Input data.frame specifying paths to mandatory bam files and reading specifics (e.g. paired status, mate, anchoring) and optional hdf5 files, as well as optional metadata. See example below for further details.
#' @param bin_size Integer. Basepair bin size to group reads/UMIs on TE consensus/gene models. Default is 25 bps.
#' @param use_gcluster Boolean. Use gcluster.run function to distribute jobs across HPC nodes. Requires misha package. If FALSE (the default), users can distribute jobs using the future infrastructure.
#' @export
addCounts <- function(scSet,
                      bams     = NULL,
                      msa_dir  = NULL,
                      bin_size = 25,
                      use_gcluster = FALSE)
{
  bin_size       <- as.integer(bin_size)
  tes            <- scSet@tes
  tes_3p         <- scSet@tes_3p
  genes          <- scSet@genes
  protocol       <- scSet@protocol
  
  # test if input files exist
  if (sum(file.exists(bams$paths)) != nrow(bams)) { stop (bams$paths[file.exists(bams_paths)], 'does not exist') }
  if (!is.null(bams$hdf5))
  {
    if (sum(file.exists(bams$hdf5)) != nrow(bams)) { stop (bams$paths[file.exists(bams_paths)], 'does not exist') }
  }
  
  # overwrite downstream results
  scSet@cstats  <- data.frame()
  scSet@pstats  <- data.frame()
  scSet@gstats  <- data.frame()
  
  if (!is.null(msa_dir))
  {
    message ('Importing conversion data.frame')
    conv_df <- fst::read_fst(paste0(msa_dir, '/', scSet@genome@provider_version, '_conv_df.fst'), as.data.table = TRUE)
    message ('Done')
  } else {
    conv_df <- NA
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
                                 bin_size     = {bin_size},
                                 conversion   = {conv_df},
                                 protocol     = '{protocol}',
                                 use_gcluster = FALSE)")
    }
    # create new empty env and fill with relevant
    empty <- new.env(parent = emptyenv())
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
    
    res <- rbindlist(lapply(res, function(x) x$retv))
  } else {
    res <- 
      compCounts(bams, 
                 tes,
                 tes_3p,
                 genes,
                 bin_size   = bin_size,
                 sense_only = TRUE,
                 conversion = conv_df,
                 protocol   = protocol)
  }
  # summarize results in case of dumb chunking
  counts <- res[, .(n = sum(n), n_all = sum(n_all)), by = c('barcode', 'name', 'id_unique', 'pos_con', 'type', 'meta')]
  
  # import 10x gene counts if defined
  if (!is.null(bams$hdf5))
  {
    message ('Importing 10x gene counts from hdf5 file(s)')
    # make hdf5 files unique
    bams <- unique(bams[, colnames(bams) %in% c('hdf5', 'meta')])
    
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
  
  scSet@counts   <- counts
  scSet@cpn      <- cpn(scSet)
  
  # add class annotation
  scSet <- annoClass(scSet)
  
  scSet@cstats   <- cstats(scSet) # needs class annotation!
    
  return(scSet)
}