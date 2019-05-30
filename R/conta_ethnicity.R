# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Return CCGA1 ethnicity probabilities
#' @return Data frame of CCGA1 ethnicity probabilities
#'
#' @export
get_ethnicity_priors <- function() {
  ethnicity_probs <- list("prob_asian" = 0.009823183,
                          "prob_african" = 0.078585462,
                          "prob_hispanic" = 0.064833006,
                          "prob_white" = 0.827111984)
  return(ethnicity_probs)
}

#' Get sample name from gt.loh.tsv file
#'
#' @param gt_loh_file_loc sample gt file location
#' @return sample id
#'
#' @export
get_sample_name <- function(gt_loh_file_loc) {
  file_splits <- strsplit(gt_loh_file_loc, "/")
  sample_file <- sapply(file_splits,
                        function(x) suppressWarnings(as.character(x[[length(x)]])))
  sample_id <- gsub(".gt.loh.tsv", "", sample_file)
  return(sample_id)
}

#' Bounds allele frequencies based on a specified amount
#'
#' Adds in allele count frequency or removes an allele count frequency
#' based on specificied value.
#'
#' @param afs allele frequency
#' @param pseudocount_freq value to bound allele frequencies by
#' @return bounded allele frequencies
#'
#' @export
bounded_freqs <- function(afs, pseudocount_freq) {
  pseudocount_afs <- ifelse(afs < 1 - pseudocount_freq,
                            afs + pseudocount_freq,
                            1- pseudocount_freq)
  return(pseudocount_afs)
}

#' Spike in or remove an allele into a parsed 1000G vcf based on population allele frequency
#'
#' Spikes in an allele count or removes an allele count to SNPs with a 0 or 1
#' population allele frequency respectively so to not drive log likelihoods
#' to -Inf or 0.
#'
#' @param parsed_vcf parsed 1000G vcf generated from 'parse_1000g_vcf()'
#' @return Data frame of parsed vcf
#' @importFrom dplyr mutate
#'
#' @export
prep_1000g_vcf_af <- function(parsed_vcf) {
  single_allele_count_freq <- 1/(1+max(parsed_vcf$AN))

  prepped_vcf <- parsed_vcf %>%
    dplyr::mutate(EAS_AF = conta::bounded_freqs(EAS_AF, single_allele_count_freq),
                  AMR_AF = conta::bounded_freqs(AMR_AF, single_allele_count_freq),
                  AFR_AF = conta::bounded_freqs(AFR_AF, single_allele_count_freq),
                  EUR_AF = conta::bounded_freqs(EUR_AF, single_allele_count_freq),
                  SAS_AF = conta::bounded_freqs(SAS_AF, single_allele_count_freq))
  return(prepped_vcf)
}

#' Create data frame with the intersection of SNPs between sample and vcf
#'
#' @param gt_loh_file_loc sample gt.loh.tsv file
#' @param prepped_vcf 1000G vcf read in and prepped by 'prep_1000g_vcf_af()'
#' @return Data frame with the intersection of gt.loh.tsv file and prepped 1000G vcf
#' @importFrom readr read_tsv
#' @importFrom dplyr select inner_join arrange
#' @importFrom tidyr drop_na
#'
#' @export
get_snp_intersection <- function(gt_loh_file_loc, prepped_vcf) {

  # Read in gt.loh.tsv file and keep only rsid and bayesian genotype calls
  gt_loh_file <- conta::read_data_table(gt_loh_file_loc) %>%
    dplyr::select(rsid, bayes_gt)

  # Take intersect between gt_file and vcf
  snp_intersect <- data.frame(dplyr::inner_join(gt_loh_file, prepped_vcf,
                                                by = c("rsid"))) %>%
    dplyr::arrange(chr) %>%
    tidyr::drop_na()
  return(snp_intersect)
}

#' Calculate probabilities per genotype and allele frequency,
#' based off of Hardy–Weinberg equilibrium.
#'
#' @param genotype genotype, either "0/0", "0/1", or "1/1"
#' @param af allele frequency, should be 0 to 1
#' @return genotype likelihood
#'
#' @export
get_gt_likelihoods <- function(genotype, af) {
  gt_likelihoods <- ifelse(genotype == "0/0", (1-af)^2,
                            ifelse(genotype == "0/1", 2 * ((1-af)*(af)),
                                   ifelse(genotype == "1/1", af^2, NA)))
  return(gt_likelihoods)
}

#' Calculate genotype likelihoods per SNP per 1000G population
#'
#' @param snp_intersect Data frame of the intersection of snps between sample
#'                      and 1000G vcf by 'get_snp_intersection()'
#' @return Data frame with the initital probabilities
#' @importFrom dplyr mutate
#'
#' @export
get_snp_likelihoods <- function(snp_intersect) {
  # Calculate P(Data|Ethnicity)
  initial_snp_likelihoods <- snp_intersect %>%
    dplyr::mutate(EAS_likelihoods = conta::get_gt_likelihoods(bayes_gt, EAS_AF),
                  AMR_likelihoods = conta::get_gt_likelihoods(bayes_gt, AMR_AF),
                  AFR_likelihoods = conta::get_gt_likelihoods(bayes_gt, AFR_AF),
                  EUR_likelihoods = conta::get_gt_likelihoods(bayes_gt, EUR_AF),
                  SAS_likelihoods = conta::get_gt_likelihoods(bayes_gt, SAS_AF))
  return(initial_snp_likelihoods)
}

#' Predict ethnicity per chromosome
#'
#' This function produces ranked ethnicity predictions per chromosome. First,
#' takes in population likelihoods per snp as output from 'get_snp_likelihoods()'.
#' Second, calculates log likelihoods, and then log likelihood per chromosome.
#' Third, applies ethnicity priors to the likelihood values to generate the
#' posterior likelihoods. Lastly, calculates the P(data) and apply it to the posterior
#' likelihoods to generate the posterior probabilities per chromosome.
#'
#' @param initial_snp_likelihoods Data frame of the initial likelihoods per
#'                                snp per 1000G populations by
#'                                'get_snp_likelihoods()'
#' @return Data frame with the log likelihoods per chromosome
#' @importFrom dplyr mutate group_by summarise select transmute arrange rename row_number
#' @importFrom tidyr gather
#'
#' @export
calculate_ranked_predictions <- function(initial_snp_likelihoods) {
  # Calculate log of P(Data|Ethnicity)
  # Takes in initial liklihoods and does a log transformation
  ll <- initial_snp_likelihoods %>%
    dplyr::mutate(EAS_ll = log10(EAS_likelihoods),
                  AMR_ll = log10(AMR_likelihoods),
                  AFR_ll = log10(AFR_likelihoods),
                  EUR_ll = log10(EUR_likelihoods),
                  SAS_ll = log10(SAS_likelihoods))

  # Calculate sum(log(P(Data|Ethnicity))) per chromosome.
  # Takes the sum of the log likelihoods per ethnicity per chromosome
  ll_chr <- ll %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(n_snps = dplyr::row_number(),
                     EAS = sum(EAS_ll),
                     AMR = sum(AMR_ll),
                     AFR = sum(AFR_ll),
                     EUR = sum(EUR_ll),
                     SAS = sum(SAS_ll))

  # Capture the number of snps per chromosome so we know how much data was used
  # to generate these probabilities
  snps_per_chr <- ll_chr %>%
    dplyr::select(chr, n_snps)

  # Calculate exp(sum(log(P(Data|Ethnicity)))) * P(Ethnicity) per chromosome.
  # Read in ethnicity priors from CCG1, exponentiate the log likelihoods and
  # multiply by the ethnicity priors to calculate posterier likelihoods
  eth_priors <- conta::get_ethnicity_priors()

  posterior_ll <- ll_chr %>%
    dplyr::transmute(EAS = exp(EAS) * eth_priors$prob_asian,
                     AMR = exp(AMR) * eth_priors$prob_hispanic,
                     AFR = exp(AFR) * eth_priors$prob_african,
                     EUR = exp(EUR) * eth_priors$prob_white,
                     SAS = exp(SAS) * eth_priors$prob_asian)

  # Calculate exp(sum(log(P(Data|Ethnicity)))) * P(Ethnicity) / P(Data) per chromosome
  # Calculate the probability of data by taking the sum of the posterior likelihoods
  # across all possible ethnicities.
  posterior_ll_num <- posterior_ll %>%
    dplyr::mutate(prob_data = rowSums(posterior_ll))

  # Each of the ethnicity's posteiror likelihoods is normalized by the
  # probability of data, thus generating posterior probabilities per ethnicity
  probs_per_chr <- posterior_ll_num %>%
    dplyr::transmute(East_Asian = (EAS/prob_data),
                     American = (AMR/prob_data),
                     African = (AFR/prob_data),
                     European = (EUR/prob_data),
                     South_Asian = (SAS/prob_data))

  probs_per_chr_final <- cbind(snps_per_chr, probs_per_chr)
  probs_per_chr_final$chr <- factor(probs_per_chr_final$chr, levels = c(1:22))

  # Rank ethnicity predictions per chromosome
  pre_ranked_predictions <- probs_per_chr_final %>%
    dplyr::select(-n_snps) %>%
    tidyr::gather(column, value, -chr) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(pred_rank = rank(-value)) %>%
    dplyr::arrange(pred_rank, as.numeric(chr)) %>%
    dplyr::rename("ethnicity" = 2,
                  "probability" = 3)

  # Join SNPs per chromosome to ranked predictions data frame. This allows better
  # for better diagnostics through understanding how much power (SNPs) went into
  # calculating ethnicity probabilities per chromosome.
  pre_ranked_predictions$chr <- as.character(pre_ranked_predictions$chr)
  ranked_predictions <- dplyr::left_join(snps_per_chr,
                                         pre_ranked_predictions,
                                         by = "chr") %>%
    dplyr::arrange(pred_rank)

  return(ranked_predictions)
}

#' Make tile plot of ranked predictions per chromosome
#'
#' Makes tile plot of ranked predictions per chromosome using ggplot
#'
#' @param ranked_predictions Data frame of ranked predictions by `calculate_ranked_predictions()`
#' @return ggplot object of the plot
#'
#' @export
plot_ranked_predictions_per_chromosome <- function(ranked_predictions) {
  g <- ranked_predictions %>%
    ggplot(aes(x = chr, y = pred_rank, fill = ethnicity)) +
    geom_tile(color = "black", size = 0.5) +
    scale_y_reverse() +
    scale_x_discrete(limits=c(1:22, "X")) +
    labs(x="Chromosome",
         y="Prediction Rank") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14, face="bold")) +
    guides(fill=guide_legend("Ethnicity:"))
  return(g)
}

#' Run the conta ethnicity prediction workflow
#'
#' Runs the full conta ethnicity prediction workflow using a sample genotype
#' file with loh regions removed and 1000G vcf file. Workflow consists of
#' reading in, parsing, and preparing 1000G vcf file, followed by taking
#' the snp intersect with a samples gt.loh.tsv file. Ethnicity probabilities are
#' calculated per snp per chromosome. Generates tables and plots to the specified
#' output directory.
#'
#' @param gt_loh_file_loc sample gt.loh.tsv file
#' @param vcf 1000G vcf
#' @param outdir Path to directory to write outputs to.
#' @return None, writes out TSVs and plots to `outdir`.
#' @importFrom readr write_tsv
#' @importFrom dplyr filter group_by count arrange top_n select bind_cols bind_rows
#'
#' @export
conta_ethnicity <- function(gt_loh_file_loc,
                            vcf,
                            outdir){

  parsed_vcf <- conta::parse_1000g_vcf(vcf, remove_sex_chr = TRUE)

  prepped_vcf <- conta::prep_1000g_vcf_af(parsed_vcf)

  snp_intersect <- conta::get_snp_intersection(gt_loh_file_loc, prepped_vcf)

  initial_snp_likelihoods <- conta::get_snp_likelihoods(snp_intersect)

  ranked_predictions <- conta::calculate_ranked_predictions(initial_snp_likelihoods)

  # Write out table of each ranked prediction per chromosome with probabilities
  sample_id <- conta::get_sample_name(gt_loh_file_loc)
  readr::write_tsv(ranked_predictions,
                   file.path(outdir, paste(sample_id,
                                           "ranked_ethnicity_predictions_per_chromosome.tsv",
                                           sep = "_")))

  # Write out table of top 2 predictions and number of chromosomes which
  # support that prediction
  top_2_predictions <- ranked_predictions %>%
    dplyr::filter(pred_rank == 1) %>%
    dplyr::group_by(pred_rank) %>%
    dplyr::count(ethnicity) %>%
    dplyr::arrange(-n) %>%
    dplyr::top_n(2) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pred_rank)

  # Set blank data frame for output template
  output_template <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                              c("best_ethnicity_pred", "best_n_chr",
                                "second_ethnicity_pred", "second_n_chr"))

  # Bind predictions from sample to output template
  # "top_2_predictions" dataframe is either a 2x2 (two ethnicities with the
  # highest number of chromosome probabilities) or 1x2 table (when a single
  # ethnicity is predicted for all chromosomes).
  if (length(top_2_predictions == 2)) {
    top_2_predictions_transformed <- dplyr::bind_cols(top_2_predictions[1, ],
                                                      top_2_predictions[2, ])
    colnames(top_2_predictions_transformed) <- c("best_ethnicity_pred", "best_n_chr",
                                                 "second_ethnicity_pred", "second_n_chr")
    pre_top_2_predictions_table <- dplyr::bind_rows(output_template, top_2_predictions_transformed)
  } else {
    colnames(top_2_predictions) <- c("best_ethnicity_pred", "best_n_chr")
    pre_top_2_predictions_table <- dplyr::bind_rows(output_template, top_2_predictions)
  }

  # Bind sample id to top 2 prediction table
  top_2_predictions_table <- cbind(sample_id, pre_top_2_predictions_table)

  # Write out table of top 2 predictions and number of chromosomes that support the call
  readr::write_tsv(top_2_predictions_table,
                   file.path(outdir, paste(sample_id,
                                           "top_2_predictions.tsv",
                                           sep = "_")))

  # Plot ranked predictions per chromosome
  ranked_prediction_plot <- plot_ranked_predictions_per_chromosome(ranked_predictions)
  png(file.path(outdir, paste(sample_id, "ethnicity_predictions_tile_plot.png", sep = "_")),
      height=720,
      width=1280,
      res = 150)
  print(ranked_prediction_plot)
  dev.off()
}
