addDist <- function(gr1, 
                    gr2,
                    n_chunks = 1)
{
  gr1$dist <- NA
  
  if (n_chunks == 1)
  {
    hits <- distanceToNearest(gr1, gr2)
    gr1[from(hits)]$dist <- mcols(hits)$distance 
  }
  
  if (n_chunks > 1)
  {
    if (is.null(names(gr1)))
    {
      names(gr1) <- 1:length(gr1)
    }
    
    # Calc distance to nearest in parallel
    results <- list()
    chunks  <- chunk(gr1, n_chunks)
    for (i in 1:length(chunks))
    {
      print(i)
      gr1_chunk    <- gr1[chunks[[i]]]
      #results[[i]] <- future({ mcols(distanceToNearest(gr1_chunk, gr2))$distance })
      results[[i]] <- future({ hits <- distanceToNearest(gr1_chunk, gr2); structure(mcols(hits)$distance, names = names(gr1_chunk[from(hits)])) })
      #results[[i]] <- future({ distanceToNearest(gr1_chunk, gr2) })
    }
    
    # add dist as metadata column
    dists                  <- unlist(sapply(results, value))
    gr1[names(dists)]$dist <- as.numeric(dists)
  }

  return(gr1)
}

addSeq <- function(BSgenome, 
                   intervals)
{
# add genomic seq
  message('Adding sequences')
  
  n_threads <- detectCores()
  
  if (n_threads > 5)
  {
    intvs_chunked <- 
      intervals %>% 
        mutate(cum_sum = cumsum(width),
               chunks = cut(cum_sum, seq(1, max(cum_sum), l = n_threads + 1), labels = FALSE, include.lowest = TRUE)) %>%
        as_tibble
    
    seqs <- plyr::dlply(intvs_chunked, "chunks", 
                        .fun = function(x) getSeq(BSgenome, as_granges(x)), 
                        .parallel = TRUE)
     
    intervals$seqs <- unlist(DNAStringSetList(seqs))
  } else {
    intervals$seqs <- getSeq(BSgenome, intervals)
  }
  
  return(intervals)
}

aggr <- function(counts, 
                 min_expr = 0,
                 sorted = FALSE)
{
  if (is.null(counts$peak))
  {
    counts[, peak := FALSE]
  }
  
  res <- 
    counts[, .(n_tot = sum(n_all, na.rm = TRUE), n_tot_uniq = sum(n)), by = c('name', 'pos_con', 'type', 'peak')]
  
  if (min_expr > 0)
  {
    expressed <- res[, .(n_tot = sum(n_tot)), by = 'name'][which(n_tot >= min_expr), ]$name
    res <- res[which(name %in% expressed), ]
  }
  
  return(res)
}

#' Genes should be curated!
calcBackground <- function(BSgenome = Hsapiens,
                           mappable = NULL,
                           tes, 
                           genes, 
                           reads,
                           min_cov = 1,
                           bin_size = 100)
{
  tes_slim  <- tes %>% select(-matches('conversion'))
  steps     <- c(0, 250 * 2^(0:25))
  
  # get intervals
  BSgenome <- GRanges(seqinfo(BSgenome))
  
  # get expressed gene ids
  genes_expr <- # get expressed genes
    genes %>%
    countOverlapsWeighted(., reads) %>% # scores multimappers
    as_tibble %>% 
    group_by(id) %>%
    mutate(covg = sum(score)) %>%
    filter(covg >= min_cov) %>%
    ungroup %>%
    as_granges()
  introns  <- getIntrons(genes_expr)
  
  # add te distance to expressed genes and mark introns
  tes_slim          <- addDist(tes_slim, genes_expr)
  tes_slim$dist_bin <- cut(tes_slim$dist, breaks = steps, include.lowest = TRUE)
  tes_slim <- tes_slim %>% mutate(intron = 1:length(.) %in% from(findOverlaps(tes_slim, introns)))
  
  if (!is.null(mappable))
  {
    BSgenome <- GenomicRanges::setdiff(BSgenome, mappable)
  }
  
  # retrieve and bin intergenic regions
  intg      <- GenomicRanges::setdiff(BSgenome, bind_ranges(granges(tes_slim), granges(genes_expr)), ignore.strand = TRUE)
  intg_bins <- unlist(tile(intg, width = bin_size))
  
  # mark introns
  intg_bins <- intg_bins %>% mutate(intron = 1:length(.) %in% from(findOverlaps(., introns)))
  
  # Calc distance to nearest exon in parallel
  intg_bins          <- addDist(intg_bins, genes_expr, n_chunks = 5)
  intg_bins$dist_bin <- cut(intg_bins$dist, breaks = steps, include.lowest = TRUE)
 
  # add read counts 
  intg_bins <- countOverlapsWeighted(intg_bins, reads) # scores multimappers
 
  # calc density
  dist_dens_conv_df <- 
    as.data.table(intg_bins)[, .(kbs = sum(width), cov = sum(score)), by = c('dist_bin', 'intron')
      ][, dens := cov / kbs
       ][, dens := ifelse(is.na(dist_bin), NA, dens)]
  
  # plotting
  # dist_dens_conv_df[order(dist_bin),] %>% 
    # filter(!is.na(dist_bin)) %>%
    # as_tibble %>% ggplot(aes(dist_bin, dens, col = intron, size = kbs)) + geom_point(size = 0.1) + geom_point() + scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # add bg density to TEs
  tes_dt           <- merge(as.data.table(tes_slim), dist_dens_conv_df, all.x = TRUE)
  dist_dens_conv_v <- structure(tes_dt$dens, names = tes_dt$id_unique)
  
  # add conversion list
  tes$dens_bg <- dist_dens_conv_v[tes$id_unique]
  
  return(tes)
}   

clustLouvain <- function(mat, k) {
  message ('Running Louvain community detection')
  
  # create data structure to be used distances between all neighbors nodes
  knn <- FNN::get.knn(as.matrix(mat), k = k)


  # Create dataframe where each record (row) represents an edge between two neighbors and has a weight = 1/(distance+1)
  #i.e: short distance --> high weight
  knn2 <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))


  #creat graph from data frame (this graph is directional weighted)
  nw <- igraph::graph_from_data_frame(knn2, directed = FALSE, vertices = NULL)
  

  nw <- igraph::simplify(nw)


  lc <- igraph::cluster_louvain(nw)


  clusters <- igraph::membership(lc)
  
  res <- data.table(id_unique = rownames(mat), module = as.integer(clusters))
  
  message (max(clusters), ' communities found')
  return (res)
}

clustSeqs <- function(seqs, 
                      kmer_length = 3,
                      cut_height = 500)
{
  kmer_mat <- oligonucleotideFrequency(seqs, width = kmer_length)
 
  d <- tgs_dist(kmer_mat)
  
  hc <- fastcluster::hclust(d, 'ward.D2')
  
  ct <- cutree(hc, h = cut_height)
  
  return(ct) 
}

countOverlapsWeighted <- function(gr1, gr2)
{
  gr1$score = 0
  
  hits <- findOverlaps(gr1, gr2)
   
  score_mat <- data.table(from = from(hits), score = gr2[to(hits)]$NH_weight)[, .(score = sum(score)), by = 'from']
  
  gr1[score_mat$from]$score <- score_mat$score
  
  return(gr1)
}

data_summary <- function(x) {
   m <- median(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

translatePeakCoords <- function(scSet)
{
  message ('Translating peak coordinates')
  peak_bins  <- unique(scSet@counts[which(type == 'gene' & peak), c('id_unique', 'name', 'pos_con')])
  genes      <- scSet@genes
  genes_expr <- genes[genes$name %in% peak_bins$name]
  
  # create 1 bp res exon intervals
  g_expr_dt <- as.data.table(unlist(tile(genes_expr, width = 1)))
  
  # annotate 1 bp intevals with id unique and binned pos con of gene
  g_expr_anno_dt <-
    g_expr_dt[, name := rep(genes_expr$name, width(genes_expr))                                   # add id unique to 1 bp intervals
      ][which(strand == '+'), pos_con := 1:.N, by = 'name'                                             # add pos con for + strand
        ][which(strand == '-'), pos_con := .N:1, by = 'name'                                           # add pos con for - strand
          ][, pos_con := cut(pos_con, breaks = seq(-1e6, 1e6, by = scSet@bin_size), include.lowest = TRUE)] # bin pos_con
  
  # join peak_bins with bp genomic cooods
  peak_coords_1bp <- merge(peak_bins, g_expr_anno_dt, by = c('name', 'pos_con'), all.x = TRUE)
  
  # reduce bp res peak intervals into genomic coordinate windows
  peak_coords <- reduce(as_granges(peak_coords_1bp))
  
  # add gene anno and only keep unique intervals (in case interval overlaps multiple gene features)
  res <- unique(join_overlap_left_directed(peak_coords, genes_expr))
  
  scSet@peaks <- res
  return(scSet)
}

flank3prime <- function(tes, 
                        genes, # should be curated
                        extend = 500)
{
  message ("Extending TE intervals at 3'")
  tes    <- tes %>% select(id_unique, name)
  tes_dt <- as.data.table(tes)
  genes  <- genes %>% select(id_unique)
  intv_blackl <- c(tes, genes)
  
  # get downstream feature index
  ds_feat <- precede(tes, intv_blackl)
  na_index <- which(is.na(ds_feat))
  ds_feat[na_index] <- 1
  
  # get distance to nearest downstream and calc extension size
  dists <- distance(tes, intv_blackl[ds_feat])
  dists[na_index] <- Inf
  dists[dists > extend] <- extend
  
  # extend downstream till nearest feature and modify position start and end intervals
  tes_flank <- 
    tes_dt[, ':=' (extend = as.integer(dists))
      ][which(strand == '+'), ':=' (start = as.integer(end + 1), end = as.integer(end + extend), position_start = 1L, position_end = extend),
        ][which(strand == '-'), ':=' (end = as.integer(start - 1), start = as.integer(start - extend), position_start = extend, position_end = 1L)
          ][which(start <= end), ]
  
  res <- as_granges(tes_flank[, !'width'][, downstream := TRUE])
  return(res)
}

mapConMSA <- function(msa, 
                      max_gap = 0.95)
{
  #if (!is.character(names(msa)) | sum(names(msa) == '') > 0) { stop ('Provide valid names in MSA') }
  suppressWarnings(if (class(msa) == "DNAStringSet") { msa <- msaToLong(msa) })
  
  al_len  <- max(msa$pos)
  n_loci  <- length(unique(msa$locus))
  name <- suppressWarnings(unlist(strsplit(msa$locus[1], '|', fixed = TRUE))[which(is.na(as.numeric(unlist(strsplit(msa$locus[1], '|', fixed = TRUE)))))])
  
  # create conversion df
  to_omit <- msa[, .(gap_perc = 1 - .N / n_loci), by = 'pos'][order(pos),]$gap_perc > max_gap
  conv_df <- 
    data.table(pos_con_orig = 1:al_len,
               omit = to_omit)[, pos_con_ungapped := .N:1 * -1, by = 'omit' # orig 1:.N
               ][which(omit), pos_con_ungapped := NA
                ][, c('pos_con_orig', 'pos_con_ungapped')]
                
  # format msa df
  msa <- msa[, 1:2]
  colnames(msa) <- c('id_unique', 'pos_con_orig')
    
  # add pos
  msa[, pos := 1:.N, by = 'id_unique']
 
  # convert con pos
  msa_long_conv <- merge(msa, conv_df, all.x = TRUE)[, pos_con := pos_con_ungapped][, !c('pos_con_orig', 'pos_con_ungapped')]
  
  # throw NA on con
  msa_long_conv_filt <- msa_long_conv[which(!is.na(pos_con)),]
  
  # add name column
  res <- msa_long_conv_filt[, name := rep(name, times = length(pos))]
  
  # create nested df
  #conv_nested <- msa_long_conv[, list(conversion = list(.SD)), by = 'id_unique']
    
  return(res)
}

# join reads by intervals
matchReads <- function(intervals, 
                       reads,
                       #resolve_multi = TRUE,
                       sense_only = FALSE)  # not supported right now
{
  # modify read columns 
  reads[, ':=' (start_read = start)]
  
    # remove duplicated read entries (can happen if neglected mate (2nd) maps to different positions). NH doesn't need to be adjusted
  reads <- unique(reads)
  
  # faster than findOverlaps and much faster than join_overlap_inner_directed!
  foverlaps.DTs <- function(DT.x, DT.y)
  {
    setkey(DT.y, seqnames, strand, start, end)
    res <- foverlaps(DT.x, DT.y, type = "any", which = TRUE, nomatch = 0L)
    return(res)
  }
  
  # if (resolve_multi)
  # {
    # if (sum(reads[which(NH > 1), ][, .(nomatch = .N != max(NH)), by = 'qname']$nomatch) > 0)
    # {
      # stop ('NH tag does not match number of alignments. Did you split the BAM by chromosome?')
    # }
  
  # }
  
  message ('Finding read vs TE/Gene overlaps')
  hits = foverlaps.DTs(reads, intervals)
    
  # join based on overlaps
  hits <- bind_cols(intervals[hits$yid, !c('seqnames', 'width')],
                    reads[hits$xid, colnames(reads) %in% c('start_read', 'barcode', 'NH', 'meta', 'start_read', 'qname'), with = FALSE])
                   
  # deduplicate reads that map to same locus (might affect pos_con mapping)
  hits_dedup <- hits[, dupl := duplicated(paste0(id_unique, qname))][, n_dupl := sum(dupl), by = 'qname'][, NH := NH - n_dupl][which(!dupl), ]
  
  
  # bla = hits_dedup[which(NH > 1), ][, .(n_loci = .N, loc_comb  = paste0(id_unique, collapse = '-'), id_unique = id_unique), by = c('barcode', 'qname')][which(n_loci > 1), ][, .(id_unique = unique(id_unique), n = length(unique(qname))), by = c('barcode', 'loc_comb')]
  
  # bla2 = hits_dedup[which(NH == 1), .(n_uniq = .N), by = c('barcode', 'id_unique')]
  
  # test <- merge(bla, bla2, by = c('barcode', 'id_unique'), all.x = TRUE)
  # test[is.na(test)] <- 0
  
  # test[, .(cor = cor(n, n_uniq, method = 'spearman')), by = c('loc_comb', 'id_unique')][which(is.na(cor)), cor := 0][order(loc_comb, cor),][which(cor > 0.1),]
  
  message ('Done')
  return(hits_dedup)
}

# computes genomic cpn and number of bps per bin for TE and gene intervals
cpn <- function(scSet,
                conv_df = NA)
{
  protocol  <- scSet@protocol
  bin_size  <- scSet@params$addCounts$bin_size
  tes       <- scSet@tes
  tes_flank <- scSet@tes_3p
  genes     <- scSet@genes
  
  if (!is.na(bin_size))
  {
    message ('Computing genomic CPN of TE bins')
    if (is.na(conv_df))
    {
      # add downstream flanking regions to TE intvs if 3' protocol was used
      if (protocol == 'threeprime')
      {
        tes          <- # combine original and downstream TE intervals
          c(
            select(tes, id_unique, name, position_start, position_end), 
            tes_flank
            ) 
      } else {
        tes$downstream <- NA
      }
      
      # combine TE and gene intervals
      intvs_dt <- 
        c(
          select(genes, id_unique, name, position_start, position_end), 
          tes
          ) %>% as.data.table
      
      # filter non-aligned intvs, and adjust alignment start, end for GRanges coverage function (start must be smaller than end)
      intvs_adj <- 
        intvs_dt[which(!is.na(position_start)), 
          ][, position_start_old := position_start
           ][which(strand == '-'), ':=' (position_start = position_end, position_end = position_start_old)
            ][, !'position_start_old']
            
      # bp genomic cpn per family/gene
      fam_bp_cpn <- 
        as.data.table(coverage(as_granges(intvs_adj[which(is.na(downstream)), .(seqnames = name, start = position_start, end = position_end, strand)])))[, .(pos_con = 1:.N, cpn = value), by = 'group_name'][, name := group_name][, !'group_name']
      
      # bp genomic cpn per family (flanking region, -1 to indicate downstream) 
      if (protocol == 'threeprime')
      {
        fam_bp_cpn_flank <- 
          as.data.table(coverage(as_granges(intvs_adj[which(!is.na(downstream)), .(seqnames = name, start = position_start, end = position_end, strand)])))[, .(pos_con = 1:.N, cpn = value), by = 'group_name'][, name := group_name][, !'group_name'][, pos_con := -1 * pos_con]
          
        fam_bp_cpn <- rbind(fam_bp_cpn, fam_bp_cpn_flank)
      }
            
      # bin and summarize counts
      fam_bin_bps <- 
        fam_bp_cpn[, pos_con := cut(pos_con, breaks = seq(-1e6, 1e6, by = bin_size), include.lowest = TRUE)
          ][, .(bps = sum(cpn)), by = c('name', 'pos_con')] 
          
      # throw bins with 0 bp coverage (not present in curated TE universe)
      fam_bin_bps <- fam_bin_bps[which(bps > 0), ]
    }
        
    if (!is.na(conv_df))
    {
      # calc cpn on chunks because of size (parallel)
      cpn_list <- 
        foreach (i = 1:(length(chunks))) %dopar%
        {
          # print(i)
          tes_chunk <- tes[chunks[[i]]]
          
          if (is.null(tes_subs$conversion)) # if conversion df missing, assumes 3' extension or rmsk alignments
          {
            chunk <- as.data.table(tes_subs)[, .(pos_con = position_start:position_end), by = c('name', 'id_unique')]
          } else {
            chunk <- rbindlist(tes_subs$conversion, idcol = 'name')
          }
          
          res_chunk <- chunk[, .(cpn = .N), by = c('name', 'pos_con')] 
          
          return(res_chunk)
        }
      
      # combine and sum cpn per family (chunking split families sometimes) and add consensus model size
      cpn_df <- rbindlist(cpn_list)[, .(cpn = sum(cpn)), by = c('name', 'pos_con')]
      
      # add NA cpn counts (NA counts are wrong for rsmk alignment based estimates)
      fam_tot_bps <- data.table(name = tes$name, bps = width(tes))[, .(bps = sum(bps)), by = 'name']
      fam_con_bps <- cpn_df[, .(bps_con = sum(cpn)), by = 'name']
      cpn_na <- merge(fam_tot_bps, fam_con_bps, by = 'name')[, .(pos_con =  NA, cpn = bps - bps_con), by = 'name']
      
      # combine na and con cpn dfs
      cpn_all <- rbind(cpn_df, cpn_na)
        
      # add con model size
      res <- cpn_all[, con_size := max(abs(pos_con), na.rm = TRUE), by = 'name']
    }
    
    # relevel bins so negative is at end
    res <- .relevelBins(fam_bin_bps)
  } else {
    res <- data.frame()
  }
  
  return(res)
}

pairwiseAlign <- function(df)
{
  sub_m      <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
  alignments <- Biostrings::pairwiseAlignment(df$seq, df$con_seq, type = 'global-local', gapExtension = 1, gapOpening = 4, substitutionMatrix = sub_m)
  
  seq_len   <- suppressWarnings(if (class(df) == 'GRanges') { width(df) } else { df$width })
  dels      <- purrr::map(stringr::str_locate_all(Biostrings::pattern(alignments), '-'), ~ .x[, 1])
  ins       <- purrr::map(stringr::str_locate_all(Biostrings::subject(alignments), '-'), ~ .x[, 1])
  align_st  <- Biostrings::start(Biostrings::subject(alignments))
  align_len <- Biostrings::nchar(alignments)
   
  align_df <- tibble::tibble(id_unique = df$id_unique, 
                             seq_len = seq_len,
                             dels = dels, 
                             ins = ins, 
                             align_len = align_len,
                             align_st = align_st)
    
  # throw misaligned repeats (no end to end alignment)
  align_df$aligned <- align_df$align_len - sapply(align_df$dels, length)
  
  align_filt_df <- align_df[align_df$aligned == align_df$seq_len, ] 
  
  res <- align_filt_df[, !colnames(align_filt_df) %in% c('seq_len', 'aligned')]
     
  return(res)
}
#######################

alignTEs <- function(tes, tmp_dir = NULL, n_jobs = 400)
{
  if (is.null(tmp_dir) & n_jobs > 1) { stop('Must provide path to temporary directory') }
  

  
  tes <- as.data.table(tes) 
  
  # distribute alignment
  if (n_jobs > 1)
  {
    # add chunk ids
    chunk_ids <- unique(stringi::stri_rand_strings(n = n_jobs, length = 12))
    
    tes <- 
      tes[, cum_sum := cumsum(as.numeric(nchar(con_seq)) + width + seq(1, n_jobs, l = nrow(tes)))
       ][, chunks   := cut(cum_sum, 
                           breaks = seq(1, max(cum_sum), l = n_jobs + 1), 
                           include.lowest = TRUE, 
                           labels = chunk_ids)
        ][, !'cum_sum']
    
    # save df to disk
    fwrite(tes, glue("{tmp_dir}/chunk_all.txt"))
        
    # generate commands    
    cmds <- list()
    for (chunk in unique(tes$chunks))
    {
      cmds[[chunk]] <- 
        glue("df <- data.table::fread(cmd = paste('grep', '{chunk}', '{tmp_dir}/chunk_all.txt'), 
                                            col.names = c({glue::glue_collapse(glue::single_quote(colnames(tes)), sep = ',')})) 
                    Repexpr:::pairwiseAlign(df)")
    }
    
    # create new empty env and fill with relevant
    empty <- new.env(parent = emptyenv())
    
    # distribute, compute and combine
    res    <- gcluster.run3(command_list = cmds, max.jobs = n_jobs, envir = empty, io_saturation = FALSE)
    res_df <- data.table::rbindlist(lapply(res, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL)
    
    file.remove(glue("{tmp_dir}/chunk_all.txt"))
  } else {
    res_df <- pairwiseAlign(tes)
  }
    
  # merge dfs
  tes_slim <- suppressWarnings(tes[, !c('seq', 'con_seq', 'chunks')])
  
  tes_anno <- 
    merge(tes_slim, res_df, all.x = TRUE, by = 'id_unique')[order(seqnames, start), ] %>%
    as_granges()
  
  return(tes_anno)
}

dedupBAM <- function(bam_files, local = TRUE, paired = FALSE, index = TRUE)
{
  old_wd <- getwd()
  bam_dir <- dirname(bam_files)
  setwd(bam_dir)
  
  commands <- list()
  for (bam_file in bam_files)
  {
      out_file <- gsub('.bam$', '_deduplicated.bam', bam_file)
      bam_name <- basename(bam_file)
      #paired <- Rsamtools::testPairedEndBam(bam_file, index = gsub('.bam$', '.bai', bam_file))
   
    if (paired == TRUE)
    {
      commands[[bam_file]] <- 
        glue::glue('system("umi_tools dedup --output-stats={bam_name} --multimapping-detection-method=NH --per-cell --extract-umi-method=tag --buffer-whole-contig --cell-tag=CR --umi-tag=UR --paired -I {bam_file} -S {out_file}")')
    } else {
      commands[[bam_file]] <- 
        glue::glue('system("umi_tools dedup --output-stats={bam_name} --multimapping-detection-method=NH --per-cell --extract-umi-method=tag --buffer-whole-contig --cell-tag=CR --umi-tag=UR -I {bam_file} -S {out_file}")')
    }
  }
  
  gpatterns::gcluster.run2(jobs_title='UMI-tools', command_list = commands, io_saturation=TRUE, max.jobs = 8)
  
  if (index)
  {
    indexBAM(dir(bam_dir, pattern = 'deduplicated', full.names = TRUE))  
  }
  setwd(old_wd)
}

indexBAM <- function(bam_files, local = TRUE, max_jobs = 20)
{
  commands <- list()
  for (bam_file in bam_files)
  {
    bam_index <- gsub('.bam$', '.bai', bam_file)
    commands[[bam_file]] <- glue::glue('system("samtools index {bam_file} {bam_index}")')
  }
  gpatterns::gcluster.run2(jobs_title='indexBAM', command_list = commands, io_saturation=TRUE, max.jobs = max_jobs)
}

mapConPos = function(tes, nested = FALSE, bin_size = 20, n_threads = 1)
{
	align_ind <- which(!is.na(tes$align_len))
  tes_filt <- as.data.table(tes)[align_ind, ]    #order(id_unique), 
 
  # data.frame with locus id and pos (basepair resolution)
  pos_all <- data.table(name   = rep(tes_filt$name, times = tes_filt$align_len),
                        id_unique = rep(tes_filt$id_unique, times = tes_filt$align_len),
                        pos       = unlist(purrr::map(tes_filt$align_len, ~ 1:.x)),
                        align_st  = rep(tes_filt$align_st, times = tes_filt$align_len))
  # pos of deletions
  dels    <- tes_filt[, .(pos = unlist(dels), seq_dels = TRUE), by = 'id_unique']
  # pos of insertions
  ins     <- tes_filt[, .(pos = unlist(ins), seq_ins = TRUE), by = 'id_unique']
  
  # combine into single df
  comb1   <- merge(pos_all, ins, all.x = TRUE)
  comb2   <- merge(comb1, dels, all.x = TRUE)
  comb2[is.na(comb2)] <- FALSE
  
  con_pos_conv_df <-
    comb2[, .(pos      = pos, 
              pos_con  = align_st[1]:(align_st[1] + .N - 1), 
              seq_dels = seq_dels,
              name  = name), 
              by = c('id_unique', 'seq_ins')
          ][, pos_con := ifelse(seq_ins, NA, pos_con)
           ][which(!seq_dels), list(name, id_unique, pos, pos_con)]
           
  # bin consensus positions
  con_pos_conv_df$pos_con_bin <- cut(con_pos_conv_df$pos_con, seq(1, max(con_pos_conv_df$pos_con, na.rm = TRUE), by = bin_size - 1), include.lowest = TRUE)
            
  if (nested)
  {
    
    res <- merge(tes, 
                 con_pos_conv_df[, list(pos_conv = list(.SD)), by = 'id_unique'][, !'name'], 
                 all.x = TRUE) %>%
           select(-dels, -ins, -align_len, -align_st)
  
  } else {
   
   res <- con_pos_conv_df
  
  }

  return(res)
}

# relevels pos con column fact so that negative levels are at the end
.relevelBins <- function(df)
{
  if (sum(grepl('-', levels(df$pos_con))) > 0)
  {
    res <- df[, pos_con := forcats::fct_relevel(pos_con, rev(grep('-', levels(pos_con), fixed = TRUE, value = TRUE)), after = Inf)]
  } else {
    res <- df
  }
  return(res)
}

splitBams <- function(bam_files, local = TRUE)
{
  for (bam_file in bam_files)
  {
    bam_dir  <- paste0(dirname(bam_file), '/')
    bam_name <- basename(bam_file)
    
    chrom_levels <- 
      system(glue::glue("samtools idxstats {bam_file}"), intern = TRUE) %>% 
      strsplit(., '\t') %>% 
      do.call(rbind, .) %>% 
      as_tibble %>%
      filter(V3 > 0) %>% 
      filter(V1 %in% c(1:22, 'X', 'Y', 'MT')) %>%
      pull(V1)
      
    if (sum(grepl('chr', chrom_levels)) == 0) { chrom_levels <- paste0('chr', chrom_levels) }

    if (!local) 
    {
      cmds <- list()
      for (chrom in chrom_levels) 
      {
        outbam <- paste0(bam_dir, gsub('.bam$', paste0('_', chrom, '.bam'), bam_name))
        cmds[[chrom]] <- glue::glue('system("samtools view -bh {bam_file} {gsub("chr", "", chrom)} > {outbam}")')
      }
      
      gpatterns::gcluster.run2(jobs_title='BamByChrom', command_list = cmds, io_saturation = TRUE, max.jobs = 30)
    }
  }
}