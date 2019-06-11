getTSS <- function(genes)
{
  transcripts <- getTranscripts(genes)
  
  tss <- mutate(anchor_5p(transcripts), width = 1)
}

getTTS <- function(genes)
{
  transcripts <- getTranscripts(genes)
  
  tts <- mutate(anchor_5p(transcripts), width = 1)
  
  return(tts)
}

getExons <- function(genes,
                     transcripts = NULL)
{
  exons <- genes[genes$type == 'exon']
  
  if (!is.null(transcripts))
  {
    exons <- exons %>% filter(transcript_id %in% transcripts)
  }
  
  return(exons)
}

getTranscripts <- function(genes)
{
  transcripts <- genes[genes$type == 'transcript']
  return(transcripts)
}

getIntrons <- function(genes)
{
  # transc <- getTranscripts(genes)
  # exons  <- getExons(genes)
  
  # if (!is.null(transcripts))
  # {
    # transc <- transc %>% filter(transcript_id %in% transcripts)
    # exons  <- exons %>% filter(transcript_id %in% transcripts)
  # } 
  
  # introns     <- setdiff_ranges_directed(transc, exons)
  
  genes_dt   <- as.data.table(genes)
  gene_intvs <- as_granges(genes_dt[, .(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1]), by = 'id'])
  
  introns <- GenomicRanges::setdiff(gene_intvs, genes, ignore.strand = FALSE)
  
  return(introns)
}