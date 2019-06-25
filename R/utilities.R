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

aggr <- function(counts, sorted = FALSE)
{
  res <- 
    counts[, .(n_tot = sum(n_all)), by = c('name', 'pos_con', 'type')]
  
  #create df for all con pos
  # full_con <- unnest(counts[, .(pos_con = list(1:(con_size[1])), repclass = repclass[1]), by = 'name'])
  
  #aggregate counts (why unique)
  # dens_df <- counts[, .(n_tot = sum(n_all), cpn = cpn[1], dens = sum(n_all) / cpn[1], n_loci = length(unique(id_unique))), by = c('name', 'pos_con')]
  
  #join
  # dens_full <- merge(full_con, dens_df, all = TRUE)[which(is.na(dens)), ':=' (dens = 0, n_tot = 0, n_loci = 0)]
  
  # if (sorted)
  # {
    # dens_full <- dens_full[order(name, pos_con), ]
  # }
  
  # return(dens_full)
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

chunk <- function(vect, n_chunks = NULL, chunk_size = 1e4)
{
  if (!is.null(n_chunks)) 
  {
    chunk_size <- floor(length(vect) / n_chunks)
  }
  
  chunks <- split(1:length(vect), ceiling(seq_along(1:length(vect)) / chunk_size))
  return(chunks)
}

countOverlapsWeighted <- function(gr1, gr2)
{
  gr1$score = 0
  
  hits <- findOverlaps(gr1, gr2)
   
  score_mat <- data.table(from = from(hits), score = gr2[to(hits)]$NH_weight)[, .(score = sum(score)), by = 'from']
  
  gr1[score_mat$from]$score <- score_mat$score
  
  return(gr1)
}

flank3prime <- function(tes, 
                        genes, # should be curated
                        extend = 500)
{
  # message ("Extending TE intervals at 3'")
  intv_blackl <- c(granges(tes, use.mcols = FALSE), unique(granges(genes, use.mcols = FALSE)))
  
  # get downstream feature index
  ds_feat <- precede(tes, intv_blackl)
  na_index <- which(is.na(ds_feat))
  ds_feat[na_index] <- 1
  
  # get distance to nearest downstream and calc extension size
  dists <- distance(tes, intv_blackl[ds_feat])
  dists[na_index] <- Inf
  dists[dists > extend] <- extend
  
  # extend 3' till nearest feature
  tes_flank <- tes %>% flank_downstream(width = dists) %>% filter(width != 0)
  
  # modify position start and end intervals
  res <- tes_flank %>% mutate(position_start = 0, position_end = width)
  
  return(tes_flank)
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

mafft = function(sequences, 
                 ids = NULL,
                 subgroup = FALSE,
                 family = 'te', 
                 save = FALSE, 
                 outdir = NULL,
                 threads = -1, # uses all
                 epsilon = 0.123, 
                 memsave = FALSE,
                 fft = TRUE,
                 overwrite = FALSE,
                 op = 1.53, # default
                 parttree = FALSE,
                 retree = 2,
                 maxiterate = 0,
                 long = FALSE,
                 auto = TRUE,
                 local = FALSE)
{
  # if not outdir saves to Repexpr package
  if (is.null(outdir)) { outdir <- system.file('extdata', 'msa', package = 'Repexpr') }
  
  # create flag vector
  memsave  <- ifelse(memsave, '--memsave', '--nomemsave')
  local    <- ifelse(local, '--localpair', '')
  fft      <- ifelse(fft, '--fft', '--nofft')
  parttree <- ifelse(parttree, '--parttree', '')
  
  flags <- glue("--thread {threads} --ep {epsilon} --op {op} {memsave} {local} {parttree} --maxiterate {maxiterate} {fft}")
  
  names(sequences) <- ids
  
  input_file <- tempfile()
  writeXStringSet(DNAStringSet(sequences), input_file, format = 'fasta')
  	
  outfile <- paste0(outdir, '/', family, '_msa.fasta')

  if (!file.exists(outfile) | overwrite)
  {
    cmd = glue('mafft {flags} {input_file} > {outfile}')
    message(cmd)
    
    system(cmd)
  } else {
    message('file exists already')
  }
  
  file.remove(input_file)
  
  msa <- importMSA(outfile, long = long)
  
  if (!save)
  {
    file.remove(outfile)
  }
  
  return(msa)
}

# join reads by intervals
matchReads <- function(intervals, 
                       reads,
                       sense_only = FALSE)
{
 # modify read columns 
  reads <- reads %>% mutate(start_read = start, strand_read = strand)
  
  hits    <- join_overlap_inner(intervals, reads) %>% mutate(sense = strand == strand_read) 
  hits_dt <- as.data.table(hits)
  
  res <- 
    hits_dt[which(strand == '+'), pos := start_read - start + 1
      ][which(strand == '-'), pos := end - start_read + 1]
  
  return(res)
}

cpn <- function(scSet, 
                conv_df  = NULL,
                n_chunks = 5,
                bin_size = 25)
{
  message ('Calculating genomic TE copy-number')
  bin_size <- scSet@bin_size
  
  if (is.null(conv_df))
  {
    tes <- scSet@tes
 
    # calc chunking steps
    chunks <- chunk(tes, n_chunks = n_chunks)

    tes_dt <- as.data.table(tes)[, al_length := position_end - position_start]
    tes_dt <- tes_dt[, al_length := ifelse(is.na(al_length), width, al_length)] # width of interval for NA alignments
    
    # create cpn universe template df (currently doesn't support NA!)
    cpn_univ_df <- 
      unique(merge(tes_dt[, .(pos_con = 1:max(position_end, na.rm = TRUE)), by = 'name'],
                   data.table(pos_con = 1:1e6, con_size_bin = cut(1:1e6, seq(-1e6, 1e6, by = bin_size), include.lowest = TRUE)), 
                   all.x = TRUE)[, .(name = name, pos_con = con_size_bin)])
    
    cpn_list <- 
      foreach (i = 1:(length(chunks))) %dopar%
      {
        # print(i)
        tes_chunk <- tes_dt[chunks[[i]], ]
        
        # create df with all consensus positions of the genome
        pos_con_all <-
          data.table(name = rep.int(tes_chunk$name, tes_chunk$al_length),  
                     pos_con = sequence(tes_chunk$al_length) + rep.int(tes_chunk$position_start-1, tes_chunk$al_length))
        
        # count how often consensus position occurs per name
        cpn_df <- pos_con_all[, .(cpn = .N), by = c('name', 'pos_con')]
      }  
    # combine results and summarize    
    cpn_df <- rbindlist(cpn_list)[, .(cpn = sum(cpn)), by = c('name', 'pos_con')]
    
    # bin and summarize counts
    cpn_df <- 
      cpn_df[, pos_con := cut(pos_con, breaks = seq(-1e6, 1e6, by = bin_size), include.lowest = TRUE)
        ][, .(bps = sum(cpn)), by = c('name', 'pos_con')]
        
    res <- merge(distinct(cpn_univ_df), cpn_df, all.x = TRUE)[which(is.na(bps)), bps := 0]
  }
  
  if (!is.null(conv_df))
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
    
  return(res)
}

consensusCov <- function(tes, n_jobs = 75, bin_size = 20)
{
  # generate commands    
  tes <- tes %>% as_tibble %>% mutate(chunks = ntile(seqnames, n_jobs))
  
  cmds <- list()
  for (chunk in unique(tes$chunks))
  {
    cmds[[chunk]] <- glue("df <- Repexpr:::cpn(tes %>% filter(chunks == {chunk}), bin_size = {bin_size})")
  }
  
  # create new empty env and fill with relevant
  empty <- new.env(parent = emptyenv())
  empty$tes <- tes
  
  # distribute, compute and combine
  res <- 
    gcluster.run3(command_list = cmds,  
                  packages = c("data.table", "plyranges", "Repexpr", "base", "stats"), 
                  max.jobs = n_jobs, 
                  envir = empty, 
                  io_saturation = FALSE)
  
  # combine and summarize
  res <- rbindlist(lapply(res, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL)[, .(cpn = sum(cpn)), by = c('name', 'pos_con_bin')]

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

gcluster.run3 <-
function (..., command_list = NULL, opt.flags = "", envir = NULL, max.jobs = 400,
    debug = FALSE, R = paste0(R.home(component = "bin"), "/R"),
    packages = NULL, jobs_title = NULL, job_names = NULL, collapse_results = FALSE,
    queue = NULL, memory = NULL, threads = NULL, io_saturation = NULL,
    queue_flag = "-q @{queue}", memory_flag = "-l mem_free=@{memory}G",
    threads_flag = "-pe threads @{threads}", io_saturation_flag = "-l io_saturation=@{io_saturation}",
    script = '/home/davidbr/tgdata/src/gpatterns/exec/sgjob.sh')
{
    qq <- GetoptLong::qq
    
    if (!is.null(command_list)) {
        commands <- purrr::map(command_list, function(x) parse(text = x))
    }
    else {
        commands <- as.list(substitute(list(...))[-1L])
    }
    if (!is.null(queue)) {
        opt.flags <- paste(opt.flags, qq(queue_flag))
    }
    if (!is.null(memory)) {
        opt.flags <- paste(opt.flags, qq(memory_flag))
    }
    if (!is.null(threads)) {
        opt.flags <- paste(opt.flags, qq(threads_flag))
    }
    if (!is.null(io_saturation)) {
        opt.flags <- paste(opt.flags, qq(io_saturation_flag))
    }
    
    if (length(commands) < 1)
        stop("Usage: gcluster.run2(..., command_list = NULL, opt.flags = \"\" max.jobs = 400, debug = FALSE)",
            call. = F)
    if (!length(system("which qsub", ignore.stderr = T, intern = T)))
        stop("gcluster.run2 must run on a host that supports Sun Grid Engine (qsub)",
            call. = F)
    .gcheckroot()
    tmp.dirname <- ""
    submitted.jobs <- c()
    tryCatch({
        tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"),
            "/tmp", sep = ""))
        if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
            stop(sprintf("Failed to create a directory %s", tmp.dirname),
                call. = F)
        cat("Preparing for distribution...\n")
        save(.GLIBDIR, file = paste(tmp.dirname, "libdir", sep = "/"))
       
        if (is.null(envir)) {  
          vars <- ls(all.names = TRUE, envir = parent.frame())
        
          envir <- parent.frame()
          while (!identical(envir, .GlobalEnv)) {
              envir <- parent.env(envir)
              if (!isNamespace(envir)) {
                  vars <- union(vars, ls(all.names = TRUE, envir = envir))
              }
          }
          suppressWarnings(save(list = vars, file = paste(tmp.dirname,
              "envir", sep = "/"), envir = parent.frame()), compress = FALSE)
        } else {
          vars <- ls(all.names = TRUE, envir = envir)
          message(vars)
          suppressWarnings(save(list = vars, file = paste(tmp.dirname,
              "envir", sep = "/"), envir = envir, compress = FALSE))
        }
        .GSGECMD <- commands
        save(.GSGECMD, file = paste(tmp.dirname, "commands",
            sep = "/"))
        opts <- options()
        save(opts, file = paste(tmp.dirname, "opts", sep = "/"))
        if (!is.null(packages)) {
            .GPACKAGES <- as.list(packages)
        }
        else {
            .GPACKAGES <- as.list(.packages())
        }
        save(.GPACKAGES, file = paste(tmp.dirname, "packages",
            sep = "/"))
        cat("Running the commands...\n")
        completed.jobs <- c()
        progress <- -1
        repeat {
            num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
            if (length(submitted.jobs) < length(commands) &&
                num.running.jobs < max.jobs) {
                istart <- length(submitted.jobs) + 1
                iend <- min(length(commands), istart + (max.jobs -
                  num.running.jobs) - 1)
                for (i in istart:iend) {
                  out.file <- sprintf("%s/%d.out", tmp.dirname,
                    i)
                  err.file <- sprintf("%s/%d.err", tmp.dirname,
                    i)
                  if (!is.null(job_names)) {
                    job.name <- job_names[i]
                  }
                  else if (!is.null(jobs_title)) {
                    job.name <- sprintf("%s_%s", jobs_title,
                      i)
                  }
                  else {
                    job.name <- sprintf("sgjob_%s", i)
                  }
                  command <- sprintf("qsub -terse -cwd -S /bin/bash -N %s -o %s -e %s -V %s %s %d '%s' '%s'",
                    job.name, out.file, err.file, opt.flags,
                    script, i, tmp.dirname, R)
                  jobid <- system(command, intern = TRUE)
                  if (length(jobid) != 1)
                    stop("Failed to run qsub", call. = FALSE)
                  if (debug)
                    cat(sprintf("\tSubmitted job %d (id: %s)\n",
                      i, jobid))
                  submitted.jobs <- c(submitted.jobs, jobid)
                }
            }
            Sys.sleep(3)
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            old.completed.jobs <- completed.jobs
            completed.jobs <- setdiff(submitted.jobs, running.jobs)
            if (debug) {
                delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
                if (length(delta.jobs) > 0) {
                  for (jobid in delta.jobs) cat(sprintf("\tJob %d (id: %s) completed\n",
                    match(jobid, submitted.jobs), jobid))
                }
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- length(completed.jobs)
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("\t%d job(s) still in progress\n",
                    length(commands) - progress))
                }
            }
            else {
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- as.integer(100 * length(completed.jobs)/length(commands))
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("%d%%...", progress))
                }
                else cat(".")
            }
        }
        if (!debug && progress != -1 && progress != 100)
            cat("100%\n")
    }, interrupt = function(interrupt) {
        cat("\n")
        stop("Command interrupted!", call. = FALSE)
    }, finally = {
        if (length(submitted.jobs) > 0) {
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            answer <- c()
            for (i in 1:length(commands)) {
                res <- list()
                res$exit.status <- NA
                res$retv <- NA
                res$stdout <- NA
                res$stderr <- NA
                if (submitted.jobs[i] %in% running.jobs)
                  res$exit.status <- "interrupted"
                else {
                  fname <- sprintf("%s/%d.retv", tmp.dirname,
                    i)
                  if (file.exists(fname)) {
                    fst::read_fst(fname)
                    res$exit.status <- "success"
                    res$retv <- retv
                  }
                  else res$exit.status <- "failure"
                }
                out.file <- sprintf("%s/%d.out", tmp.dirname,
                  i)
                if (file.exists(out.file)) {
                  f <- file(out.file, "rc")
                  res$stdout <- readChar(f, 1e+06)
                  close(f)
                }
                err.file <- sprintf("%s/%d.err", tmp.dirname,
                  i)
                if (file.exists(err.file)) {
                  f <- file(err.file, "rc")
                  res$stderr <- readChar(f, 1e+06)
                  close(f)
                }
                answer[[i]] <- res
            }
            for (job in running.jobs) system(sprintf("qdel %s",
                job), ignore.stderr = T, intern = T)
            unlink(tmp.dirname, recursive = TRUE)
            if (collapse_results) {
                canswer <- tryCatch(data.table::rbindlist(lapply(answer, function(x) x$retv), use.names = FALSE, fill = FALSE, idcol = NULL),
                  error = function(e) {
                    message("returning original output due to an error. collapse your reults manually (are all the parts data frames?)")
                    return(NULL)
                  })
                if (!is.null(canswer)) {
                  return(canswer)
                }
            }
            return(answer)
        }
        unlink(tmp.dirname, recursive = TRUE)
    })
}
