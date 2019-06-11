#' Extends TE bins
#'
#' @param bins Data.frame. TE bin intervals and annotation.
#' @param downstream Integer. Basepairs to extend downstream.
#' @param binsize Integer. Tile size of extended intervals.
#'
addBinExtension <- function(bins, 
                            downstream = 1000, 
                            binsize = 20, 
                            min_dist = 100,
                            tss = GRanges(), 
                            polyA = GRanges())
{
  intv_col <- makeGRangesFromDataFrame(.collapseBins(bins), keep.extra.columns = TRUE)
  
  # throw nested repeat intervals
  intv_col_filt <- intv_col[countOverlaps(intv_col) == 1]
  
  if (length(tss) > 0 | length(polyA) > 0)
  {
    # combine blacklist intervals
    rep_tss_polyA <- c(makeGRangesFromDataFrame(intv_col), makeGRangesFromDataFrame(tss), makeGRangesFromDataFrame(polyA)) # makeGRangesFromDataFrame to drop mcol
    
    # get downstream blacklist feature
    intv_col_filt$down <- precede(intv_col_filt, rep_tss_polyA, ignore.strand = FALSE) 
    intv_col_filt$dist <- NA
    
    # get distance to downstream blacklist feature
    intv_col_filt[!is.na(intv_col_filt$down)]$dist <- 
      distance(intv_col_filt[!is.na(intv_col_filt$down)], 
               rep_tss_polyA[na.omit(intv_col_filt$down)])
    
    # set extend threshold by downstream bps value or next blacklist feature
    intv_col_filt$extend <- ifelse(is.na(intv_col_filt$dist) | (intv_col_filt$dist - min_dist) > downstream, downstream, intv_col_filt$dist)
  
    # get flanking regions until downstream or blacklist feature
    intv_ext <- flank(intv_col_filt, width = intv_col_filt$extend, start = FALSE)
    
    # filter flanking regions of width < 2 (basically no flanking regions due to neighbouring TSS, polyA or TE)
    intv_ext <- intv_ext[width(intv_ext) > 1]
    
    # filter flanking regions overlapping tss/polyA (precede ignores overlapping features)
    intv_ext <- intv_ext[!1:length(intv_ext) %in% findOverlaps(intv_ext, c(makeGRangesFromDataFrame(tss), makeGRangesFromDataFrame(polyA)))@from]
    #
  } else {
    intv_ext <- flank(intv_col_filt, width = downstream, start = FALSE)
  }
    
  # bin intervals
  message('Binning extended regions')
  intv_ext_bins_list   <- tile(intv_ext, width = binsize)
  intv_ext_bins        <- unlist(intv_ext_bins_list)
  intv_ext_bins$row_id <- rep(intv_ext$row_id,  S4Vectors::elementNROWS(intv_ext_bins_list)) 
  intv_ext_bins$repname <- rep(intv_ext$repname,  S4Vectors::elementNROWS(intv_ext_bins_list)) 
  
  # add bin_id
  intv_ext_bins_plus         <- as.data.table(intv_ext_bins[strand(intv_ext_bins) == '+'])
  intv_ext_bins_plus$bin_id  <- intv_ext_bins_plus[, .(bin_id = -1:-.N), by = 'row_id']$bin_id
  intv_ext_bins_minus        <- as.data.table(intv_ext_bins[strand(intv_ext_bins) == '-'])
  intv_ext_bins_minus$bin_id <- intv_ext_bins_minus[, .(bin_id = -.N:-1), by = 'row_id']$bin_id
  
  # combine and format
  res <-
    bind_rows(intv_ext_bins_plus, intv_ext_bins_minus) %>%
    as_tibble %>%
    dplyr::rename(chrom = seqnames, 
                  bin_size = width) %>% 
    bind_rows(bins, .)
  
  return(res)
}

addBinStats <- function(bins)
{
	bins <- as_tibble(bins)
  
  # add median bin size and copy-number stats
  bins <-
    bins %>%
    group_by(repname, bin_id) %>%
    mutate(bin_cpn      = n(),
           bin_size_med = median(bin_size)) %>% 
    group_by(repname) %>%
    mutate(fam_cpn = median(bin_cpn)) %>%
    ungroup
}

addNormCounts <- function(bins, reg = 0)
{
  if (class(bins) == 'GRanges')
  {
    bins_gr <- bins
    bins <- 
      bins %>%
      as_tibble %>%
      mutate(umis_norm = (umis + reg) / bin_size_med / bin_cpn * fam_cpn) %>%
      group_by(repname) %>%
        mutate(umis_tot      = sum(umis),
               umis_norm_tot = sum(umis_norm)) %>%
      ungroup
      
    bins_gr$umis_norm     <- bins$umis_norm
    bins_gr$umis_tot      <- bins$umis_tot
    bins_gr$umis_norm_tot <- bins$umis_norm_tot
    
    return(bins_gr)
  } else {
    bins <- 
      bins %>%
      mutate(umis_norm = (umis + reg) / bin_size_med / bin_cpn * fam_cpn) %>%
      group_by(repname) %>%
        mutate(umis_tot      = sum(umis),
               umis_norm_tot = sum(umis_norm)) %>%
      ungroup
    return(bins)
  }
}





