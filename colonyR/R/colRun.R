#' @title Run Colony
#'
#' @description Run a Colony pedigree reconstruction analysis. Some description
#'   of the parameters are taken from the ColonyUserGuide.pdf and colony2.dat
#'   provided with Colony, as downloaded on 2014-09-27.
#'
#' @details From the pdf manual accompanying Colony: "Update allele frequency:
#'   Allele frequencies are required in calculating the likelihood of a
#'   configuration. These frequencies can be provided by the user (see below)
#'   or are calculated by Colony using the genotypes in OFS, CMS (optional) and
#'   CFS (optional). In the latter case, you can ask Colony to update allele
#'   frequency estimates by taking into account of the inferred sibship and
#'   parentage relationships during the process of searching for the maximum
#'   likelihood configuration. However, updating allele frequencies could
#'   increase computational time substantially, and may not improve
#'   relationship inference much if the genetic structure of your sample is not
#'   strong (i.e. family sizes small and evenly distributed, most candidates
#'   are not assigned parentage). I suggest not updating allele frequencies
#'   except when family sizes (unknown) are suspected to be large (relative to
#'   sample size) and highly variable."
#'
#' @param ids vector giving the offspring ids
#' @param genotypes data frame with the genotypes. Each row is an individual,
#'   in the same order as for ids. Each pair of columns represents one locus.
#' @param marker_names a vector of strings with the marker names
#' @param random_seed used by Colony
#' @param update_allele_frequency boolean, if \code{TRUE} Colony will update
#'   allele frequencies by taking into account of reconstructed pedigrees.
#' @param inbreeding boolean, absence of presence of inbreeding in Colony
#'   analysis
#' @param length_of_runs Give a value of 1, 2, 3, 4 to indicate short, medium,
#'   long, very long run.
#' @param full_likelihood_precision 1/2/3=low/medium/high Precision for
#'   Fulllikelihood 
#' @param monitor_interval interval between records (number of iterations)
#' @param save_file boolean, save the results to a file?
#'
#' @export
#'


colRun = function(ids, genotypes,
  marker_names,
  random_seed = 1234,
  update_allele_frequency = FALSE,
  inbreeding = FALSE,
  length_of_runs = 1,
  full_likelihood_precision = 3,
  monitor_interval = 10000,
  save_file = FALSE) {

  library(uuid)
  temp_dir = UUIDgenerate(use.time = F)
  
  # create a working directory
  dir.create(temp_dir)
  
  # prepare the data
  id_table = colPrepData(ids = ids,
    genotypes = genotypes,
    marker_names = marker_names,
    dataset_name = "colonyFromR",
    output_file_name = file.path(temp_dir, "colonyFromR"),
    random_seed = random_seed,
    update_allele_frequency = update_allele_frequency,
    inbreeding = inbreeding,
    number_of_runs = 1,
    length_of_runs = length_of_runs,
    full_likelihood_precision = full_likelihood_precision,
    monitor_interval = monitor_interval,
    file = file.path(temp_dir, "colony.dat"))
  
  # run Colony
  system(paste0("colony IFN:", file.path(temp_dir, "colony.dat")), intern = T)
  
  # parse the data
  d = colParseResults(file = "colonyFromR", dir = temp_dir)
  
  # delete the working directory
  files = list.files(temp_dir)
  file.remove(file.path(temp_dir, files))
  file.remove(temp_dir)
  
  # prepare the results
  d[["ids"]] = id_table
  d[["temp_dir"]] = temp_dir
  d[["input_genotypes"]] = genotypes
  d[["input_markers_names"]] = marker_names
  d[["input_ids"]] = ids
  
  # save the results if needed
  if (save_file) {
    file_out = paste0("colony_run_", temp_dir, ".rda")
    colony_run_results = d
    save(colony_run_results, file = file_out)
    d[["file"]] = file_out
  }

  # return
  return(d)
}
