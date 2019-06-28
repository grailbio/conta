# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Run conta source detection
#'
#' Reads a base folder with each conta result under a folder. Each conta run
#' should contain one *.conta.tsv and *.gt.tsv that correspond to conta
#' results and genotypes/contaminants of that sample. Conta source detection
#' will identify the samples for which conta made a call, and then compare its
#' contaminant reads against genotypes of every other samples' genotypes.
#'
#' For each genotype comparison, maximum likelihood that is obtained by using
#' the genotype of the contaminant as the contamination probability
#' is compared against the original maximum likelihood that is obtained by
#' using the MAF as the prior. If this likelihood is higher, then the sample
#' with the highest likelihood is reported as the candidate source of
#' contamination. Second and third highest candidates are also reported.
#'
#' @param base character base directory with a batch of conta results
#' @param out_file character output file
#' @param batch_samples file that contains samples to limit the analysis to
#' @param subfolder subfolder name if conta results are in a subfolder
#' @param threshold re-call conta based on a new threshold
#' @param blackswan blackswan term for maximum likelihood estimation
#' @param outlier_frac fraction of outliers to remove
#' @param source_threshold difference between MAF llr and genotype llr to call
#' @param cores number of cores to use for calculations
#'
#' @return none
#'
#' @importFrom utils write.table
#' @export
conta_source <- function(base, out_file, batch_samples = NA,
                         subfolder = "", threshold = NA, blackswan = 1,
                         outlier_frac = 0.001, source_threshold = 0.01,
                         cores = 8, precision = 3) {
  options("mc.cores" = cores)

  # Add slash to the end of base if it wasn' specified
  base <- ifelse(!endsWith(base, "/"), paste(base, "/", sep = ""), base)

  # s3_ls reads both files and paths under a folder, we only need folders
  if (startsWith(base, "s3")) {
    if (requireNamespace("grails3r")) {
      files <- grails3r::s3_ls(base)$path
    } else {
      files <- get_s3_folders(base)
    }
    paths <- files[endsWith(files, "/")]
  } else {
    paths <- paste(list.dirs(base, recursive = FALSE, full.names = FALSE),
                  "/", sep = "")
  }

  # Read conta results and genotypes, good to cache them at this point
  subfolder <- ifelse(subfolder == "", subfolder,
                      ifelse( !endsWith(subfolder, "/"),
                              paste(subfolder, "/", sep = ""),
                              subfolder))
  # Limit to batch files
  if (!is.na(batch_samples)) {
    bs <- read_data_table(batch_samples, sep = "\t", header = F)
    bs <- rbind(bs, data.table( V1 = tools::file_path_sans_ext(
                                basename(batch_samples))))
    bs <- unique(bs)
    paths <- paths[sapply(1:length(paths),
                          function(i) sum(startsWith(paths[i], bs$V1)) > 0)]
  }
  run_names <- basename(paths)

  # Get all conta results
  lres <- parallel::mclapply(paste(base, paths, subfolder, run_names,
                                   ".conta.tsv", sep = ""), read_data_table)
  names(lres) <- run_names

  # Skip results with NA calls, and correct paths and names
  lres <- lres[unlist(lapply(lres, function(x) {!is.na(x$conta_call)}))]
  run_names <- names(lres)
  paths <- paste(run_names, "/", sep = "")

  # Re-call if necessary
  if (!is.na(threshold)) {
    for (i in 1:length(lres)) {
      lres[[i]]$conta_call <-
        ifelse(lres[[i]]$avg_log_lr >= threshold, TRUE, FALSE)
    }
  }

  # Keep only the same batch smples if batch files variable is specified
  gt_files <- paste(base, paths, subfolder, run_names,
                    ".gt.loh.tsv", sep = "")
  gt_source_files <- paste(dirname(out_file), "/", run_names,
                           ".gt.loh.source.tsv", sep = "")
  lgt <- parallel::mclapply(gt_files, read_data_table)

  # Define a data.frame for the output, one line for contaminated sample
  out <- data.frame(
    conta_version = character(),
    name = character(),
    cf = numeric(),
    source_call = logical(),
    avg_maf_lr = numeric(),
    best_gt_lr = numeric(),
    best_sample = character(),
    second_gt_lr = numeric(),
    second_sample = character(),
    third_gt_lr = numeric(),
    third_sample = character())

  # for each conta result file
  for (i in 1:length(lres)) {

    res <- lres[[i]]

    # find source regardless of whether it is called contaminated or not
    if (TRUE) {
      # contamination fraction, only test for originally discovered fraction
      cf <- res$cf

      # avg likelihood ratio obtained by using population frequencies
      avg_maf_lr <- res$avg_log_lr

      # read the genotypes and variant reads for this sample
      # we only need hom alleles for host sample
      gt1 <- lgt[[i]][gt != "0/1"]

      # calculate source likelihood scores for each of the other samples
      # TODO: Report number of SNPs used for source call
      scores <- parallel::mclapply(c(1:length(lres)), function(j) {

        # Keep the hets for candidate source
        gt2 <- lgt[[j]]

        conc <- get_genotype_concordance(gt1, gt2)

        # If no concordance or same sample return NA
        if (is.na(conc) | conc >= 0.95) {
          return(NA)
        }

        # Calculate likelihood ratio using genotypes as prior
        return(get_source_llr(gt1, gt2, cf, blackswan, outlier_frac))
        })

      # Combine and sort the scores, remove self results
      scores <- unlist(scores)
      score_df <- data.frame(score = scores, name = run_names)
      score_df <- score_df[order(score_df$score, decreasing = TRUE), ]
      score_df[is.na(score_df$score), ]$name <- NA
      score_df <- score_df[!is.na(score_df$name), ]

      # Make a source contamination call
      # Requirements are:
      # 1) A top scoring contaminant exists (non-NA result)
      # 2) Conta made a contamination call (without the use of sources)
      # 3) Contamination fraction is larger than zero
      # 4) Source score exceeds the original conta score by a specified fraction
      call_thr <- avg_maf_lr + source_threshold
      source_call <- !is.na(score_df[1, ]$score) &&
        res$conta_call && cf > 0 &&
        score_df[1, ]$score > call_thr

      # Create base output frame
      out_this <- data.frame(
        conta_version = as.character(packageVersion("conta")),
        name = run_names[i],
        cf = as.numeric(cf),
        source_call = source_call,
        avg_maf_lr = as.numeric(avg_maf_lr),
        best_gt_lr = as.numeric(NA),
        best_sample = NA,
        second_gt_lr = as.numeric(NA),
        second_sample = NA,
        third_gt_lr = as.numeric(NA),
        third_sample = NA)

      # Create genotype lr frame, copy of original data at first
      gtm <- lgt[[i]]

      if (nrow(score_df) >= 1) {
        out_this$best_gt_lr <- max(score_df[1, ]$score, 0)
        out_this$best_sample <- score_df[1, ]$name
        dat1 <- get_source_llr(gt1,
                               lgt[[which(run_names == score_df[1, ]$name)]],
                               cf, blackswan, outlier_frac,
                               detailed_results = TRUE)
        dat1s <- dat1[, c("rsid", "lr")] %>%
          dplyr::mutate(source_lr1 = specify_precision(lr, precision)) %>%
          dplyr::select(-lr)
        gtm <- merge(gtm, dat1s, by = "rsid", all.x = TRUE)
      }

      # Add second best result if there is a second result
      if (nrow(score_df) >= 2) {
        out_this$second_gt_lr <- max(score_df[2, ]$score, 0)
        out_this$second_sample <- score_df[2, ]$name
        dat2 <- get_source_llr(gt1,
                               lgt[[which(run_names == score_df[2, ]$name)]],
                               cf, blackswan, outlier_frac,
                               detailed_results = TRUE)
        dat2s <- dat2[, c("rsid", "lr")] %>%
          dplyr::mutate(source_lr2 = specify_precision(lr, precision)) %>%
          dplyr::select(-lr)
        gtm <- merge(gtm, dat2s, by = "rsid", all.x = TRUE)
      }

      # Add third best result if there is a third result
      if (nrow(score_df) >= 3) {
        out_this$third_gt_lr <- max(score_df[3, ]$score, 0)
        out_this$third_sample <- score_df[3, ]$name
        dat3 <- get_source_llr(gt1,
                               lgt[[which(run_names == score_df[3, ]$name)]],
                               cf, blackswan, outlier_frac,
                               detailed_results = TRUE)
        dat3s <- dat3[, c("rsid", "lr")] %>%
          dplyr::mutate(source_lr3 = specify_precision(lr, precision)) %>%
          dplyr::select(-lr)
        gtm <- merge(gtm, dat3s, by = "rsid", all.x = TRUE)
      }

      # Remove het snps and arrange rsid in ascending numeric order
      gtm <- gtm %>%
        dplyr::filter(gt != "0/1") %>%
        dplyr::mutate(id = as.numeric(stringr::str_replace(rsid, "rs", ""))) %>%
        dplyr::arrange(id) %>%
        dplyr::select(-id)

      # Add NA likelihoods for missing values
      if (!("source_lr1" %in% colnames(gtm))) {
        gtm <- gtm %>% dplyr::mutate(source_lr1 = NA)
      }
      if (!("source_lr2" %in% colnames(gtm))) {
        gtm <- gtm %>% dplyr::mutate(source_lr2 = NA)
      }
      if (!("source_lr3" %in% colnames(gtm))) {
        gtm <- gtm %>% dplyr::mutate(source_lr3 = NA)
      }

      # Merge this sample's result with overall results
      out <- rbind(out, out_this)

      gtm <- gtm %>%
        dplyr::mutate(cp = specify_precision(cp, precision))
      write_data_table(gtm, gt_source_files[i])
    }
  }

  # Format output precision
  out <- as.data.frame(out) %>%
    dplyr::mutate(cf = specify_precision(cf, precision)) %>%
    dplyr::mutate(avg_maf_lr = specify_precision(avg_maf_lr, precision)) %>%
    dplyr::mutate(best_gt_lr = specify_precision(best_gt_lr, precision)) %>%
    dplyr::mutate(second_gt_lr = specify_precision(second_gt_lr, precision)) %>%
    dplyr::mutate(third_gt_lr = specify_precision(third_gt_lr, precision))
  write_data_table(out, out_file)
}
