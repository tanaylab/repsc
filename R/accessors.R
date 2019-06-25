#' @include scSet.R
NULL

#' Returns transposable element families from scSet
#'
#' Retrieves all available families in the provided scSet that match defined characteristics.
#'
#' @param scSet scSet object to query.
#' @param pattern Regular expression. Returns only families that match \code{pattern}.
#' @param fixed Logical. If \code{TRUE}, \code{pattern} is matched as is.
#' @return Character vector.
#' @examples
#' \dontrun{
#' gs <- createscSet(genome = Hsapiens, tes = te_annotation_df)
#' families <- names(gs, pattern = 'LTR12') # returns all families containing 'LTR12' in their name.
#'
#' # Using a pattern not found in scSet will suggest similar hits, e.g.:
#' names(gs, pattern = 'LTR5HS')
#'}
#' @seealso [createscSet()]
#' @export
repnames <- function(scSet, 
                     pattern = NULL,
                     fixed = FALSE)
{
  if (length(pattern) > 1) { stop ('Pattern must be of length 1') }
  families <- unique(scSet@tes$name)
  if (!is.null(pattern))
  {
    families <- grep(pattern, families, fixed = fixed, value = TRUE)
  }
  
  if (length(families) == 0)
  {
    families <- unique(scSet@tes$name)
    fam_dists <- structure(stringdist::stringdist(pattern, unique(scSet@tes$name)), names = families)
    suggestions <- paste(names(head(sort(fam_dists), 5)), collapse = ', ')
    stop (glue::glue("{pattern} not found in scSet. Did you mean any of {suggestions}?"))
  }
  
  return(families)
}

#' Returns gene intervals from scSet
#' @export
setGeneric("genes", function(scSet) standardGeneric("genes"))
setMethod("genes", signature("scSet"), function(scSet) {
  out <- scSet@genes
  return(out)
})

#' Returns TE intervals from scSet
#' @export
setGeneric("tes", function(scSet) standardGeneric("tes"))
setMethod("tes", signature("scSet"), function(scSet) {
  out <- scSet@tes
  return(out)
})