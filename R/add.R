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

addCounts <- function(scSet,
                      reads    = NULL,
                      chunk    = FALSE,
                      msa_dir  = NULL,
                      bin_size = 25)
{
  bin_size       <- as.integer(bin_size)
  scSet@bin_size <- bin_size
  tes            <- scSet@tes
  genes          <- scSet@genes
  protocol       <- scSet@protocol
  
  if (!is.null(msa_dir))
  {
    message ('Importing conversion data.frame')
    conv_df <- fst::read_fst(paste0(msa_dir, '/', scSet@genome@provider_version, '_conv_df.fst'), as.data.table = TRUE)
    message ('Done')
  } else {
    conv_df <- NULL
  }
  
  if (!chunk)
  {
    # import reads
    reads <- importBAM(reads, anchor = protocol)
    
    #' Count number of reads per interval and cell barcode
    counts <- 
      compCounts(reads, 
                 tes,
                 genes,
                 bin_size   = bin_size,
                 sense_only = FALSE,
                 conversion = conv_df,
                 protocol   = protocol)
  }
  
  # add cpn df
  scSet@cpn      <- cpn(scSet, bin_size = bin_size)
  scSet@counts   <- counts
  
  return(scSet)
}