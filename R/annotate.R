#' Add TE and gene class annotation to counts
#' @export
annoClass <- function(scSet = NULL,
                      n_levels = 3)
{
  f <- function(df, tes, genes)
  {
    if (is.null(df$class) & sum(dim(df)) > 0)
    {
      fam_class  <- unique(as.data.table(tes)[, c('name', 'class')])
      gene_class <- unique(as.data.table(genes)[, c('name', 'class')])
    
      df_anno <- merge(df, rbind(fam_class, gene_class), all.x = TRUE, by = 'name')
      
      # refactor levels to reduce complexity
      top_classes  <- unique(df_anno[which(type == 'gene'), c('name', 'class')]) %>% count(class) %>% top_n(n_levels, n) %>% pull(class)
      res <- df_anno[, class := as.character(forcats::fct_other(class, keep = c('DNA', 'LTR', 'SINE', 'LINE', top_classes)))]
      return(res)
    } else {
      return(df)
    }
  }
  
  tes    <- tes(scSet)
  genes  <- genes(scSet)
  counts <- scSet@counts # don't do counts(scSet)!
  gstats <- scSet@gstats
  pstats <- scSet@pstats

  scSet@counts <- f(counts, tes, genes)
  scSet@gstats <- f(gstats, tes, genes)
  scSet@pstats <- f(pstats, tes, genes)
  return(scSet)
}                      