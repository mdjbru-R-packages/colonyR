#' @title Prepare data for analysis with Colony
#'
#' @description Prepare the input file for a Colony run. Some description of
#'   the parameters are taken from the ColonyUserGuide.pdf and colony2.dat
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
#' @param dataset_name for Colony output
#' @param output_file_name prefix used by Colony for the output files
#' @param random_seed used by Colony
#' @param update_allele_frequency boolean, if \code{TRUE} Colony will update
#'   allele frequencies by taking into account of reconstructed pedigrees.
#' @param inbreeding boolean, absence of presence of inbreeding in Colony
#'   analysis
#' @param number_of_runs integer
#' @param length_of_runs Give a value of 1, 2, 3, 4 to indicate short, medium,
#'   long, very long run.
#' @param full_likelihood_precision 1/2/3=low/medium/high Precision for
#'   Fulllikelihood 
#' @param monitor_interval interval between records (number of iterations)
#' @param file file name for the input file produced by this function
#'
#' @return Nothing but writes an input file for Colony
#'
#' @examples
#' # (this requires a library with trout data from our group)
#' library(PKDtroutR)
#' d = subset(trout, trout$field_trip_number == "september" &
#'       trout$full_origin == "vainupea")
#' microsats = c( "Ssosl438_1", "Ssosl438_2", "Ssosl311_1", "Ssosl311_2",
#'   "Str15inraP_1", "Str15inraP_2", "LG.14_1_1", "LG.14_1_2", "Str543inraP_1",
#'   "Str543inraP_2", "Ssa197_1", "Ssa197_2", "LG.15_1_1", "LG.15_1_2",
#'   "Strutta.58_1", "Strutta.58_2", "Str60inra_1", "Str60inra_2", "Str73inra_1",
#'   "Str73inra_2", "Ssosl417_1", "Ssosl417_2", "Str85inraP_1", "Str85inraP_2",
#'   "LG.10_2_1", "LG.10_2_2", "Bs131_1", "Bs131_2", "Ssa407_1", "Ssa407_2" )
#' marker_names = strsplit(microsats[2 * 1:(length(microsats) / 2)], "_")
#' marker_names = unlist(lapply(marker_names, function(x) x[1]))
#' d[, c("fish_global_id", microsats)]
#' ids = d$fish_global_id
#' genotypes = d[, microsats]
#' colPrepData(ids, genotypes, marker_names, monitor_interval = 1000,
#'             length_of_runs = 1)
#' #run
#' a = system("colony IFN:colony.input.toto", intern = T)
#' # load results
#' r = read.table("colonyFromR.MidResult", header = T, comment.char = "")
#' plot(r$NumIterate, r$BtLogL, type = "l", col = "cornflowerblue")
#' lines(r$NumIterate, r$CrLogL, col = "brown")
#'
#' @export
#'

colPrepData = function(ids, genotypes,
  marker_names,
  dataset_name = "colonyFromR",
  output_file_name = "colonyFromR",
  random_seed = 1234,
  update_allele_frequency = FALSE,
  inbreeding = FALSE,
  number_of_runs = 1,
  length_of_runs = 1,
  full_likelihood_precision = 3,
  monitor_interval = 10000,
  file = "colony.input.dat") {

  # get parameters
  dataset_name = paste0("'", dataset_name, "'")
  output_file_name = paste0("'", output_file_name, "'")
  n_offsprings = length(ids)
  stopifnot(n_offsprings == nrow(genotypes))
  n_loci = ncol(genotypes) / 2
  stopifnot(ncol(genotypes) %% 2 == 0)
  update_allele_frequency = as.numeric(update_allele_frequency)
  di_mono_species = 2
  inbreeding = as.numeric(inbreeding)
  diplo_haplodiplo_species = 0
  poly_monogamy_sire_dam = c(0, 0)
  clone_inference = 0
  sibship_prior = 0
  known_pop_allele_freq = 0
  monitor_results_method = 0
  non_windows_version = 1
  full_likelihood = 1
  stopifnot(length(marker_names) == 1/2 * ncol(genotypes))
  marker_types = rep(0, n_loci)
  allelic_dropout_rate = rep(0, n_loci)
  false_allele_rate = rep(0.0001, n_loci)
  prob_sire_dam_included = c(0, 0)
  n_candidates_sire_dam = c(0, 0)
  n_known_paternity = 0
  n_known_maternity = 0
  n_known_paternal_sibship = 0
  n_known_maternal_sibship = 0
  n_excluded_paternity = 0
  n_excluded_maternity = 0
  n_excluded_paternal_sibships = 0
  n_excluded_maternal_sibships = 0
  
  # open the file
  f = file(description = file, open = "w")

  # function to write to the file
  w = function(x, comment = "") {
    cat(as.character(x), paste0("! ", comment), "\n", file = f)
  }

  # function to write the genotypes
  w_genotypes = function(ids, genotypes) {
    o = ""
    genotypes[is.na(genotypes)] = 0
    g = data.frame(lapply(genotypes, as.character), stringsAsFactors = F)
    for (i in 1:length(ids)) {
      l = paste(paste0("O_", as.character(ids[i])),
        paste(g[i, ], collapse = " "), sep = " ")
      o = paste0(o, l, "\n")
    }
    cat(o, file = f)
  }
  
  # write to the file
  w(dataset_name, "Dataset name")
  w(output_file_name, "Output file name")
  w(n_offsprings, "Number of offspring on the sample")
  w(n_loci, "Number of loci")
  w(random_seed, "Seed for random number generator")
  w(update_allele_frequency, "0/1=Not updating/updating allele frequency")
  w(di_mono_species, "2/1=Dioecious/Monoecious species")
  w(inbreeding, "0/1=No inbreeding/inbreeding")
  w(diplo_haplodiplo_species, "0/1=Diploid species/HaploDiploid species")
  w(poly_monogamy_sire_dam, "0/1=Polygamy/Monogamy for males & females")
  w(clone_inference, "0/1=Clone inference =No/Yes")
  w(sibship_prior, "0,1,2,3=No,weak,medium,strong sibship size prior; mean paternal & meteral sibship size")
  w(known_pop_allele_freq, "0/1=Unknown/Known population allele frequency")
  w(number_of_runs, "Number of runs")
  w(length_of_runs, "Length of run")
  w(monitor_results_method, "0/1=Monitor method by Iterate#/Time in second")
  w(monitor_interval, "Monitor interval in Iterate# / in seconds")
  w(non_windows_version, "non-Windows version")
  w(full_likelihood, "Fulllikelihood")
  w(full_likelihood_precision, "1/2/3=low/medium/high Precision for Fulllikelihood")
  w(marker_names, "Marker names")
  w(marker_types, "Marker types, 0/1 = codominant/dominant")
  w(allelic_dropout_rate, "Allelic dropout rate")
  w(format(false_allele_rate, scientific = FALSE), "false allele rate")
  w_genotypes(ids = ids, genotypes = genotypes)
  w(prob_sire_dam_included, "prob. of dad/mum included in the candidates")
  w(n_candidates_sire_dam, "numbers of candidate males & females")
  w(n_known_paternity, "Number of known paternity")
  w(n_known_maternity, "Number of known maternity")
  w(n_known_paternal_sibship, "Number of known paternal sibship")
  w(n_known_maternal_sibship, "Number of known maternal sibship")
  w(n_excluded_paternity, "Number of offspring with known excluded paternity")
  w(n_excluded_maternity, "Number of offspring with known excluded maternity")
  w(n_excluded_paternal_sibships, "Number of offspring with known excluded paternal sibships")
  w(n_excluded_maternal_sibships, "Number of offspring with known excluded maternal sibships")
  
  # close file
  close(f)
  
}
