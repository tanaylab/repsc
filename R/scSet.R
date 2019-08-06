#' An S4 class to store genomic annotations and results of Repguide functions.
#'
#' @slot genome BSgenome object.
#' @slot tes GRanges object
#' @slot cis GRanges object
#' @slot refdir Character
#' @slot tempdir Character
#' @slot targets GRanges object
#' @slot blacklist GRanges object
#' @slot whitelist GRanges object
#' @slot plots List
#' @slot guide_length Integer
#' @slot families Character
#' @slot PAM Character
#' @slot kmers GRanges
#' @slot alignments DNAStringSetList
#' @slot consensus DNAStringSet
#' @slot combinations Data.frame
#' @slot calls List
#' @slot n_cores Integer
#' @slot seed Integer
#' @exportClass scSet
scSet <- 
  setClass('scSet', 
           slots = c(
                                 genome     = 'BSgenome',
                                 protocol   = 'character', 
                                 tes        = 'GRanges',
                                 tes_3p      = 'GRanges',
                                 mdata      = 'data.frame',
                                 genes      = 'GRanges',
                                 tss        = 'GRanges',
                                 cpn        = 'data.frame',
                                 cells      = 'character',
                                 peaks      = 'GRanges',
                                 cstats     = 'data.frame',
                                 gstats     = 'data.frame',
                                 pstats     = 'data.frame',
                                 msa_dir    = 'character',
                                 ds         = 'integer',
                                 counts     = 'data.frame',
                                 counts_ds  = 'data.frame',
                                 bin_size   = 'integer',
                                 plots      = 'list', 
                                 alignments = 'DNAStringSetList',
                                 consensus  = 'DNAStringSet',
                                 calls      = 'list',
                                 n_cores    = 'numeric',
                                 seed       = 'numeric'
                                )
          )

setMethod(
  "initialize",
  signature = "scSet",
  function(
           .Object, 
           genome,
           protocol,
           tes,
           genes,
           reads,
           n_cores,
           seed
           ) 
    {        
      if (class(genome) != 'BSgenome') { stop ('Provide a BSgenome object as genome') }
      if (is.null(tes)) { stop ('Provide transposable element annotation of genome') }
      
      # Detect number of cores
      if (is.null(n_cores)) { n_cores <- parallel::detectCores() - 5 }
      
      # Import 
      genes     <- rtracklayer::import(genes)
      GenomeInfoDb::seqlevelsStyle(genes) <- 'UCSC'
      tes       <- Reputils::importRMSK(tes)
      tss       <- Reputils::getTSS(genes) %>% select(gene_name)
      
      # Curate
      genes_cur <- Reputils::curateGenes(genes)        # only retains exons, resolves naming ambiguity
      tes_cur   <- Reputils::curateTEs(tes,            # tries to resolve nested TE intervals and adds annotation
                                       genes_cur) 
      if (protocol == 'threeprime')
      {
        tes_3p  <- flank3prime(tes_cur, genes_cur, extend = 500) 
      } else {
        tes_3p  <- GRanges()
      }
      
      
      if (is.null(tes$name) | is.null(tes$id_unique)) { stop ('TE annotation requires name and id_unique metadata columns') }
      
      if (is.element('doMC', installed.packages())[1]) 
      { 
        doMC::registerDoMC(n_cores) 
      } else {
        warning ('Package doMC not found. Install package to support multicore usage')
        n_cores <- 1
      }
      
      .Object@genome      <- genome
      .Object@protocol    <- protocol
      .Object@tes         <- tes_cur
      .Object@tes_3p       <- tes_3p
      .Object@genes       <- genes_cur
      .Object@tss         <- tss
      .Object@ds          <- 1000L
      .Object@plots       <- list()
      .Object@n_cores     <- n_cores  
      .Object@seed        <- seed
            
      timestamp(prefix = 'Created new scSet at ', suffix = '', quiet = FALSE)
      return(.Object)
    }
)

setMethod("show", "scSet", function(object) 
{ 
  genome_id  <- object@genome@pkgname
  #n_barcodes <- object@counts$barcode
  #n_cells    <- object@n_cells
  n_tes      <- length(object@tes)
  n_fams     <- length(unique(object@tes$name))
  n_plots    <- length(object@plots)
  n_cores    <- object@n_cores
  
  #message(glue::glue('scSet with {n_cells} cells'))
  message(glue::glue('scSet object of {genome_id}'))
  message(glue::glue('with {n_tes} loci from {n_fams} families'))
  message(glue::glue('with {n_plots} QC plots'))
  message(glue::glue('registered {n_cores} cores'))
})

#' Exports results from scSet
#' 
#' @param scSet scSet containing the results
#' @param outdir String. Creates a new directory with timestamp prefix in \code{outdir} and exports results. If \code{NULL} and \code{force = TRUE}, the folder is created in the current working directory.
#' @param full Logical. If \code{TRUE}, additionally exports stats for all guideRNAs and combinations (instead of only selected guides).
#' @param force Logical. If \code{TRUE} and \code{outdir = NULL}, writes output to new folder in current working directory.
#' @param workspace Logical. If \code{FALSE} (the default), suppresses additional export of \code{scSet} as .RData file.
#' @param dpi Integer. Resolution of exported images. 
#' @examples
#' \dontrun{
#' gs <- createscSet(Hsapiens, tes = te_annotation_df)
#' gs <- addTargets(gs, targets = 'LTR13')
#' gs <- addGuides(gs, guide_length = 16, n_mismatches = 0, gc_content = c(0.25, 0.9), n_clust = 12)
#' gs <- plotGuides(gs)
#' export(gs, outdir = NULL, force = TRUE) # Creates new folder in current working directory and exports results  
#' }
#' @export
setGeneric('export', function(scSet, ...) standardGeneric('export'), signature = 'scSet') 

#' Create new scSet object
#'
#' @param genome BSgenome object (required). Target genome assembly stored as BSgenome object.
#' @param alt_chromosomes Logical. If \code{FALSE} (the default), restricts genome annotation to the main chromosome assembly.
#' @param tes Path to repeatmasker output file or GRanges object with 'repname' metacolumn (required).
#' @param cis Path to bed file with cis regulatory feature coordinates or GRanges object (optional).
#' @param blacklist Path to bed file with blacklisted regions or GRanges object (optional). Guides binding to \code{blacklist} regions are blacklisted.
#' @param whitelist Path to bed file with whitelisted regions or GRanges object (optional). Guide off-target binding to \code{whitelist} regions are scored neutrally.
#' @param temp Path to directory where temporary files are stored. Needs to be read- and writeable with sufficient storage space for large file sizes. If \code{NULL} (the default), is set to the return value of \code{tempdir()}.
#' @param n_cores Integer. Number of cores to use for downstream functions. If \code{NULL} (the default), detects the number of cores automatically. The [doMC](https://cran.r-project.org/web/packages/doMC/index.html) packge must be installed to register the cores.
#' @param refdir Path to search for bowtie index files. Will create new indeces in \code{refdir} if no corresponding files are found (i.e. do not match BSgenome prefix). If empty (the default), searches in the bowtie_indeces directory of the Repguide installation path.
#' @param seed Integer. Seed for the random number generator. 19 by default.
#' @return scSet object.
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38) 
#'
#' Path to directory containing BSgenome.Hsapiens.UCSC.hg38 bowtie indeces (e.g. BSgenome.Hsapiens.UCSC.hg38.1.ebwt, ...)
#' indexdir <-  system.file(package = 'Repguide', 'bowtie_indeces') 
#'
#' # Path to TE annotation file
#' te_anno <- system.file(package = 'Repguide', 'extdata', 'hg38_ucsc_rmsk_ltr.txt.gz')
#' 
#' gs <- createscSet(genome = BSgenome.Hsapiens.UCSC.hg38, tes = te_anno, refdir = indexdir)
#' gs  
#' }
#' @seealso [BSgenome::available.genomes()], [bowtie manual](http://bowtie-bio.sourceforge.net/manual.shtml), and [repeatmasker](http://www.repeatmasker.org/)
#' @export
createScSet <- function(genome, 
                        protocol = NULL,
                        tes = NULL, 
                        genes = NULL, 
                        reads = NULL,
                        n_cores = NULL,
                        seed = 19)
{
  new('scSet', genome, protocol, tes, genes, reads, n_cores, seed)
}                         