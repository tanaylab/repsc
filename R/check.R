#' Check validity of input
#' @export
checkInput <- function(df)
{
  if (sum(class(df) == 'data.frame') == 0)
  {
    stop ('Input must be a data.frame')
  }
  
  cols_inval <- setdiff(colnames(df), c('paths', 'hdf5', 'paired', 'mate', 'barcode', 'meta', 'chunk'))
  if (length(cols_inval) > 0)
  {
    warning ('Unrecognized columns in input data.frame (will be ignored): ', cols_inval)
  }
  
  if (is.null(df$paths))
  {
    stop ('No "paths" column to alignment files provided')
  }

  if (is.null(df$paired))
  {
    stop ('No "paired" column defining library status found')
  }
  
  if (is.null(df$barcode))
  {
    stop ('No "barcode" column defining how to import barcode found, e.g. "CB" for 10x')
  }
  
  if (!is.null(df$hdf5) & !requireNamespace('hdf5r', quietly = TRUE))
  {
    stop ('hdf5 files provided but package "hdf5r" not found, please install')
  }
  
  if (is.null(df$mate) & sum(df$paired) > 0)
  {
    stop ('No "mate" column defining what mate to import from alignments found')
  }
  
  if (sum(df$paired) > 0 & sum(is.na(df$mate)) > 0)
  {
    stop ('For paired alignments, "mate" column defining what mate to import required')
  }
  
  if (is.null(df$chunk))
  {
    warning ('No "chunk" column defining how to parallelize input processing found')
  }
  
  if (is.null(df$meta))
  {
    warning ('No "meta" column defining sample relationship/annotation found')
  }

  # test if input files exist
  bams_exist  <- file.exists(df$paths)
  if (sum(bams_exist) != nrow(df)) 
  { 
    stop (bams$paths[!bams_exist], 'do(es) not exist') 
  }
  if (!is.null(df$hdf5))
  {
    hdf5_exists <- file.exists(df$hdf5) 
    if (sum(hdf5_exists) != nrow(df))
    {
      stop (df$hdf5[!hdf5_exists], 'do(es) not exist') 
    }
  }
  
  # modify for printing
  df$paths <- basename(df$paths)
  if (!is.null(df$hdf5)) 
  {
    df$hdf5 <- basename(df$hdf5)
  }
  print (head(df))
  message ('Input is valid')
}

checkCounts <- function(counts)
{
  




}