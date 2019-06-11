msaToLong <- function(msa)
{
  n_alignments <- as.numeric(length(msa))
  n_positions  <- as.numeric(width(msa)[1])
  msa_size     <- n_alignments * n_positions
  
  if (msa_size < 1.5e7)
  {
    mat  <- as.matrix(msa)
    colnames(mat) <- 1:ncol(mat)
    
    res <- 
      data.table(locus = rep(rownames(mat), ncol(mat)), 
                 pos   = as.numeric(rep(colnames(mat), each = nrow(mat))), 
                 nuc = as.vector(mat))[which(nuc != '-'), ]
  } else {  
    chunk_size <- max(1e4, round(n_alignments / 20))
    results    <- list()
    chunks     <- split(1:length(msa), ceiling(seq_along(1:length(msa)) / chunk_size))
    for (i in 1:length(chunks)) 
    {
      chunk <- msa[chunks[[i]]]
      results[[i]] <- future::future({
                                      rbindlist(lapply(chunk, function(locus) {data.table(pos=1:length(locus), nuc=as.vector(locus))[which(nuc!='-'),]}), idcol = 'locus')
                                     })
    }

    res <- rbindlist(lapply(results, future::value))
  }
  
  return(res)
}

