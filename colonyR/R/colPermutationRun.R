#' @title Run Colony on permuted genotypes
#'
#' @description Run Colony on multiple permuted genotypes. See
#'   \code{\link{colRun}} for the details about common parameters.
#'
#' @param n_perm number of permutations to run
#' @param n_cores number of cores to use
#' @param save_file boolean, if \code{TRUE} one result file is saved for each 
#'   permutation (usefu to follow the progress and in case of crash)
#' @param useLocalColonyBinary boolean, use \code{./colony} instead of
#'   \code{colony}?
#'
#' @return A list with the output from \code{\link{colRun}}
#'
#' @export
#'

colPermutationRun = function(ids, genotypes,
  marker_names,
  random_seed = 1234,
  update_allele_frequency = FALSE,
  inbreeding = FALSE,
  length_of_runs = 1,
  full_likelihood_precision = 3,
  monitor_interval = 10000,
  n_perm = 10,
  n_cores = 1,
  save_file = FALSE,
  useLocalColonyBinary = FALSE) {

  # create the permuted genotypes
  permuted_genotypes = list()
  for (i in 1:n_perm) {
    g = genotypes
    stopifnot(ncol(g) %% 2 == 0)
    for (j in 1:(ncol(g) / 2)) {
      alleles = unlist(g[, c(2*j-1, 2*j)])
      stopifnot(length(alleles) == 2*nrow(g))
      alleles = sample(alleles)
      g[, (2*j-1)] = alleles[1:nrow(g)]
      g[, (2*j)] = alleles[(nrow(g) + 1) : (2*nrow(g))]
    }
    permuted_genotypes[[i]] = g
  }

  # analysis function
  analysis_function = function(g) {
    r = colRun(ids = ids,
      genotypes = g,
      marker_names = marker_names,
      random_seed = random_seed,
      update_allele_frequency = update_allele_frequency,
      inbreeding = inbreeding,
      length_of_runs = length_of_runs,
      full_likelihood_precision = full_likelihood_precision,
      monitor_interval = monitor_interval,
      save_file = save_file,
      useLocalColonyBinary = useLocalColonyBinary)
    return(r)
  }

  # run
  if (n_cores == 1) {
    results = lapply(permuted_genotypes,
      analysis_function)
  } else {
    stopifnot(n_cores > 1)
    library(snow)
    cluster = makeSOCKcluster(rep("localhost", n_cores))
    to.export = list("ids", "marker_names", "random_seed",
      "update_allele_frequency", "inbreeding", "length_of_runs",
      "full_likelihood_precision", "monitor_interval",
      "save_file")
    clusterExport(cluster, to.export, envir = environment())
    clusterEvalQ(cluster, {library(colonyR)})
    clusterSetupRNG(cluster)
    results = parLapply(cluster, permuted_genotypes, analysis_function)
    stopCluster(cluster)
  }

  # return
  return(results)
}
