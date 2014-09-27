# from the ColonyUserGuide.pdf and colony2.dat provided with Colony, as
# downloaded on 2014-09-27
# @param length_of_runs Give a value of 1, 2, 3, 4 to indicate short, medium,
#   long, very long run.
# @param full_likelihood_precision 1/2/3=low/medium/high Precision for
#   Fulllikelihood 
# @param monitor_interval interval between records (number of iterations)


colPrepData = function(ids, genotypes,
  marker_names,
  dataset_name = "colonyFromR",
  output_file_name = "colonyFromR",
  random_seed = 1234,
  update_allele_frequency = TRUE,
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
