#' @include scSet.R
NULL

setMethod("export", signature('scSet'), function(scSet, 
                                                 outdir    = NULL, 
                                                 by_meta   = FALSE,
                                                 output    = 'sparse',
                                                 intervals = FALSE,
                                                 force     = FALSE, 
                                                 workspace = FALSE, 
                                                 dpi = 300) 
{ 
  if(is.null(outdir) & !force) 
  { 
    stop('No outdir provided. Either specify the outdir or use force = TRUE to save reports to the current working directory') 
  }
  oldwd <- getwd()
  
  counts <- Repsc::counts(scSet)
  gstats <- scSet@gstats
  pstats <- scSet@pstats
  cstats <- scSet@cstats
  
  # create output directory
  datum <- gsub('-', '', Sys.Date())
  zeit <- gsub(' ', '', gsub(':', '', format(Sys.time(), "%X")))

  suffix <- paste0('Repsc_output', '_', datum, '_', zeit, '/')
  outdir <- ifelse(is.null(outdir), suffix, paste0(outdir, '/', suffix))
  dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  setwd(outdir)
  message(paste0('Exporting scSet to ', outdir))
  
  # export curated TE/gene intervals
  if (intervals)
  {
    tes   <- as.data.table(tes(scSet))
    genes <- as.data.table(genes(scSet))
    
    fwrite(tes,   paste0(outdir, '/tes_curated.txt'), sep = '\t')
    fwrite(genes, paste0(outdir, '/genes_curated.txt'), sep = '\t')
  }
    
  # export counts as sparse matrices
  if (output == 'sparse')
  {
    if (by_meta)
    {
      # factorizes loci, genes, and fams so output is comparable between metas
      counts_te   <- counts[which(type == 'te'), ][, ':=' (id_unique = as.factor(id_unique), name = as.factor(name))]
      counts_gene <- counts[which(type == 'gene'), ][, name := as.factor(name)]
      
      for (meta_id in unique(counts$meta))
      {
        # factorizes barcodes with genic levels, so TE/genic output is comparable
        counts_gene_meta <- counts_gene[which(meta == meta_id), ][, barcode := as.factor(barcode)]
        counts_te_meta   <- counts_te[which(meta == meta_id), ][, barcode := factor(barcode, levels = levels(counts_gene_meta$barcode))]
      
        # export
        .exportSparse(counts_te_meta, counts_gene_meta, outdir = outdir, prefix = paste0(meta_id, '_'))
      }
    } else {
      # factorizes barcodes so cells are identical between genes, loci, and fams
      counts <- counts[, barcode := as.factor(paste(barcode, meta, sep = '|'))]
      
      counts_te   <- counts[which(type == 'te'), ]
      counts_gene <- counts[which(type == 'gene'), ]
      
      # export
      .exportSparse(counts_te, counts_gene, outdir = outdir)
    }
  }
  
  # export counts in long format
  if (output == 'long')
  {
    fwrite(counts, paste0(outdir, '/counts.tsv'), sep = '\t')
  }
  
  # export stats
  fwrite(cstats, paste0(outdir, '/cell_stats.tsv'), sep = '\t')
  fwrite(pstats, paste0(outdir, '/peak_stats.tsv'), sep = '\t')
  fwrite(gstats, paste0(outdir, '/feature_stats.tsv'), sep = '\t')

  # export workspace
  if(workspace) { save(scSet, file = 'workspace.RData') }
  
  # export history
  savehistory(file = 'command_history.txt')
  
  setwd(oldwd) 
})     

.exportSparse <- function(counts_te,
                          counts_gene,
                          outdir = NULL,
                          prefix = '')
{
  message ('Exporting sparse matrices')
  outdir_meta <- paste0(outdir, prefix)
  
  # format to matrix market format for export
  smat_te_loc <- longToSparse(counts_te[,   c('id_unique', 'barcode', 'n')])
  smat_te_fam <- longToSparse(counts_te[,   c('name', 'barcode', 'n_all')])
  smat_gene   <- longToSparse(counts_gene[, c('name', 'barcode', 'n')])
  
  barcodes      <- levels(counts_gene$barcode)
  if (prefix != '') 
  {
    barcodes_orig <- sapply(strsplit(barcodes, '|', fixed = TRUE), function(x) glue::glue_collapse(x[-max(2, length(x))], sep = '|')) # removes last meta barcode suffix
  } else {
    barcodes_orig <- barcodes
  }
  genes         <- rownames(smat_gene)
  te_loci       <- rownames(smat_te_loc)
  te_fams       <- rownames(smat_te_fam)
   
  Reputils::writeMM2(smat_te_loc, paste0(outdir_meta, 'te_loci_counts.mtx'))
  Reputils::writeMM2(smat_te_fam, paste0(outdir_meta, 'te_fam_counts.mtx'))
  Reputils::writeMM2(smat_gene,   paste0(outdir_meta, 'gene_counts.mtx'))
  
  fwrite(list(barcodes_orig), paste0(outdir_meta, 'barcodes.tsv'), sep = '\t')
  fwrite(list(genes),         paste0(outdir,      'genes.tsv'),    sep = '\t')
  fwrite(list(te_loci),       paste0(outdir,      'te_loci.tsv'),  sep = '\t')
  fwrite(list(te_fams),       paste0(outdir,      'te_fams.tsv'),  sep = '\t')
    
  message ('Done')
}

.exportPlots <- function(scSet, dpi = 300)
{
 # export target plots
 i <- 1
 for (p in guideSet@plots$targets)
 {
    fn <- paste0('targets', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
 
 # export guide plots
 i <- 1
 for (p in guideSet@plots$guides)
 {
    fn <- paste0('guides', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
 
 # export guide plots
 i <- 1
 for (p in guideSet@plots$combinations)
 {
    fn <- paste0('combinations', '_', i, '.png')
    ggsave(fn, p, device = 'png', dpi = dpi)
    i <- i + 1
 }
}