bamToCounts <- function(bam_files,
                        genes, 
                        tes, 
                        type = '10x', 
                        protocol = 'fiveprime')
{
  if (protocol == 'fiveprime')
  {
    mate   <- 'first'
    paired <- TRUE  
    genes  <- curateGenes(genes, extend = 100)
    anchor <- 'fiveprime'
  }
  
  if (protocol == 'threeprime')
  {
    mate   <- NA
    paired <- FALSE
    genes  <- curateGenes(genes, extend = 0)
    anchor <- 'threeprime'
  }
  
  # import reads
  reads <- importBAM(bam_files, bai_pos = 'suffix', paired = paired, mate = mate, anchor = anchor)

  # count reads in features
  gene_counts <- compCounts(reads, intervals = genes, tidy = TRUE)
  te_counts   <- compCounts(reads, intervals = tes, tidy = TRUE)





}