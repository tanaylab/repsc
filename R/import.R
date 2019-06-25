importBAM <- function(bam_df, 
                      bai_pos = 'suffix', 
                      paired = NULL, 
                      mate = NA,
                      anchor = NULL,
                      filt = TRUE,
                      barcode = 'CB')
{
  bam_files <- as.character(bam_df$paths)
  paired    <- bam_df$paired
  mate      <- as.character(bam_df$mate)
  barcode   <- as.character(bam_df$barcode)
  
  if (sum(dirname(bam_files) == '.') > 0) { stop ('Please provide full path to files') }
  if (is.null(paired))                    { paired <- unique(Rsamtools::testPairedEndBam(bam_files)) }
   
  # get matching indicies
  if (bai_pos == 'suffix')
  {
    bam_indices <- gsub('.bai$', '', bam_files)
  } else {
    bam_indices <- gsub('.bam$', '', bam_files) 
  }
  
  # actual reading function
  readBamF <- function(bam_file,
                       bam_index,
                       paired,
                       barcode, 
                       mate)
  {
    message(glue("Reading {bam_file}"))
    if (barcode == 'CB')
    {
      res <- GenomicAlignments::readGAlignments(bam_file, 
                                                index = bam_index,
                                                use.names = FALSE,
                                                param = Rsamtools::ScanBamParam(tag = c("NH", 'CB'),
                                                                                flag=Rsamtools::scanBamFlag(isPaired = paired, 
                                                                                                 isProperPair     = paired, 
                                                                                                 isFirstMateRead  = mate == 'first',
                                                                                                 isSecondMateRead = mate == 'second')))
      if (class(res) == 'list') { res <- unlist(res) }
      res <- granges(res, use.mcols = TRUE, use.names = TRUE)      
    } else {
      res <- GenomicAlignments::readGAlignments(bam_file, 
                                                index = bam_index,
                                                use.names = FALSE,
                                                param = Rsamtools::ScanBamParam(tag = c("NH"),
                                                                                flag=Rsamtools::scanBamFlag(isPaired = paired, 
                                                                                                 isProperPair     = paired, 
                                                                                                 isFirstMateRead  = mate == 'first',
                                                                                                 isSecondMateRead = mate == 'second')))    
      if (class(res) == 'list') { res <- unlist(res) }  
      res <- granges(res, use.mcols = TRUE, use.names = TRUE)      
      res$barcode <- barcode
    }                                                                                            
    return(res)
  }
  
  # import in parallel
  message ('Importing BAM file(s)')
  results <- list()
  for (i in 1:length(bam_files))
  {
    results[[i]] <- future({ readBamF(bam_files[i], bam_indices[i], barcode = barcode[i], paired = paired[i], mate = mate[i]) })
  }
  reads <- unlist(as(lapply(results, value), 'GRangesList'))
  message ('Done')

  # add NH weighted column and rename CB to barcode and filter bad barcodes
  reads <- 
      reads %>%
      mutate(NH_weight = 1 / NH, 
             NH_flag = ifelse(NH > 1, 'multi', 'unique'))
             
  if (!is.null(reads$CB))
  {
    reads$barcode <- reads$CB
    reads <- reads %>% select(-CB)
  }
  
  if (filt)
  {
    reads <-
      reads %>%
      filter(!is.na(barcode)) #%>%   
      #filter(!is.na(UB))
  }
  
  # make UCSC seqlevelsStyle and filter chroms
  GenomeInfoDb::seqlevelsStyle(reads) <- 'UCSC'       
  reads <- GenomeInfoDb::keepStandardChromosomes(reads, pruning.mode = 'coarse')
  
  # anchor coordinates
  if (!is.null(anchor))
  {
    if (anchor == 'fiveprime')
    {
      reads <- mutate(anchor_5p(reads), width = 1)
    }
    if (anchor == 'threeprime')
    {
      reads <- mutate(anchor_3p(reads), width = 1)
    }
  }
  
  return(reads)
} 

importDFAM <- function(path)          
{
  dfam_hits <- fread(path)
  colnames(dfam_hits) <- gsub('-', '_', colnames(dfam_hits))
  colnames(dfam_hits) <- gsub('#', '', colnames(dfam_hits))
  
  res <- 
    dfam_hits %>% 
      mutate(start = ifelse(strand == '+', ali_st, ali_en),
             end = ifelse(strand == '+', ali_en, ali_st),
             seqnames = seq_name, 
             name = family_name,
             id_unique = paste(1:nrow(.), name, sep = '|')) %>%
      select(-ali_st, -ali_en, -env_st, -env_en, -seq_name, -family_name, -family_acc) %>%    
      as_granges()
    
  return(res)
}

importMSA <- function(fnames,
                      long = TRUE,
                      n_jobs = 500)
{
  # f_sizes <- 
    # tibble(fname = fnames, size = file.info(fnames)$size) %>%
    # arrange(size) %>%
    # mutate(chunk = cut(cumsum(size), breaks = seq(1, max(cumsum(size)) + 1e9, by = 1e9), labels = FALSE, include.lowest = TRUE))
  # test <- future.apply::future_by(f_sizes[1:5,], f_sizes$chunk[1:5], function(x) readDNAStringSet(x$fname))
  
  if (length(fnames) > 1)
  {
    # create commands for distribution
    cmds <- list()
    for (fname in fnames)
    {
      if (long)
      {
        cmds[[basename(fname)]] <- glue("Repexpr:::msaToLong(readDNAStringSet('{fname}'))")
      } else {
        cmds[[basename(fname)]] <- glue("readDNAStringSet('{fname}')")
      }
    }
  
    # create new empty env and fill with relevant
    empty <- new.env(parent = emptyenv())
    
    # distribute, compute and combine
    res_gclust <- 
      gcluster.run3(command_list = cmds,  
                    packages = c("Biostrings", "base", "stats"), 
                    job_names = names(cmds),
                    max.jobs = n_jobs, 
                    envir = empty, 
                    io_saturation = FALSE)
    
    res <- rbindlist(lapply(res_gclust, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL)
  } else {
    if (long)
    {
      res <- msaToLong(readDNAStringSet(fnames))
    } else {
      res <- readDNAStringSet(fnames)
    }
  }
  
  return(res)
}
  
importRMSK <- function(path, 
                       assembly = 'main',
                       main_classes = TRUE,
                       proper_alignments = TRUE)          
{
  rmsk_hits <- 
    fread(path, skip = 3, 
          fill = TRUE,
          header = FALSE,
          data.table = TRUE,
          col.names = c('sw_score', 'perc_div', 'perc_del', 'perc_ins', 'seqnames', 'start', 'end', 'left_chrom', 'strand', 'name', 'repclass_family', 
								        'position_start', 'position_end', 'left_rep', 'id_unique'))
  
  rmsk_hits[, c("repclass", "repfamily") := tstrsplit(repclass_family, "/", fixed=TRUE)]
  rmsk_hits[, strand := gsub('C', '-', strand, fixed = TRUE)]
  
  # keep only major repclasses
  if (main_classes)
  {
    rmsk_hits <- rmsk_hits[which(repclass %in% c('DNA', 'SINE', 'LTR', 'LINE')), ]
  }

  # adjust position on alignment based on strand
  rmsk_hits <- rmsk_hits[which(strand == '-'), position_start := (left_rep)][, !'left_rep']
  rmsk_hits$position_start <- as.integer(rmsk_hits$position_start)
  
  if (proper_alignments)
  {
    rmsk_hits <- rmsk_hits[which(position_start < position_end), ]
  }
   
  res <-
    rmsk_hits %>%
    mutate(id_unique = paste(name, 1:nrow(.), sep = '|')) %>%
    select(-left_chrom, -repclass_family) %>%    
    as_granges()
    
  if (assembly == 'main')
  {
    res <- GenomeInfoDb::keepStandardChromosomes(res, pruning.mode = 'coarse')
  }
    
  return(res)
}