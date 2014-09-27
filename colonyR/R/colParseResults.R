#' @title Parse results from a Colony run
#'
#' @description Parse the monitoring information and final best configuration
#'   for a Colony run
#' @param file basename of the output file name given to Colony
#' @param dir directory where the run took place
#'
#' @return A list with two elements: \code{run} with the monitoring
#'   information and \code{pedigree} with the final pedigree
#'
#' @export
#'

colParseResults = function(file, dir = "") {
  f = file.path(dir, file)
  # load the analysis monitoring
  fr = paste0(f, ".MidResult")
  r = read.table(fr, header = T, comment.char = "")
  # load the family data
  fc = paste0(f, ".BestConfig")
  c = read.table(fc, header = T, comment.char = "")
  # output
  o = list()
  o[["run"]] = r
  o[["pedigree"]] = c
  return(o)
}
