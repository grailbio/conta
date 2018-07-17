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
#' @param cores number of cores to use for calculations
#'
#' @return none
#'
#' @importFrom utils write.table
#' @export
conta_source <- function(base, out_file, batch_samples = NA,
                         subfolder = "", threshold = NA, blackswan = 0.05,
                         cores = 8) {

  options("digits" = 5)
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
  names <- basename(paths)

  lres <- mclapply(paste(base, paths, subfolder, names, ".conta.tsv", sep = ""),
                 read_data_table)

  # Re-call if necessary
  if (!is.na(threshold)) {
    for (i in 1:length(lres)) {
      lres[[i]]$conta_call <-
        ifelse(lres[[i]]$avg_log_lr >= threshold, TRUE, FALSE)
    }
  }

  # Keep only the same batch smples if batch files variable is specified
  lgt <- mclapply(paste(base, paths, subfolder, names, ".gt.loh.tsv", sep = ""),
                read_data_table)

  # Define a data.frame for the output, one line for contaminated sample
  out <- data.frame(name = character(),
                    cf = numeric(),
                    source_call = logical(),
                    avg_maf_lr = numeric(),
                    best_conf = numeric(),
                    best_sample = character(),
                    second_conf = numeric(),
                    second_sample = character(),
                    third_conf = numeric(),
                    third_sample = character())

  # for each conta result file
  for (i in 1:length(lres)) {

    res <- lres[[i]]

    # find source regardless of whether it is called contaminated or not
    if (TRUE) {
      # contamination fraction
      cf <- res$cf

      # avg likelihood ratio obtained by using population frequencies
      avg_maf_lr <- res$avg_log_lr

      # read the genotypes and variant reads for this sample
      # we only need hom alleles for host sample
      gt1 <- lgt[[i]][gt != "0/1"]

      # calculate source likelihood scores for each of the other samples
      # TODO: Report number of SNPs used for source call
      scores <- mclapply(c(1:length(lres)), function(j) {

        # Keep the hets for candidate source
        gt2 <- lgt[[j]]

        conc <- get_genotype_concordance(gt1, gt2)

        # If no concordance or same sample return NA
        if (is.na(conc) | conc >= 0.95) {
          return(NA)
        }

        # Calculate likelihood ratio using genotypes as prior
        avg_gt_lr <- get_source_lr(gt1, gt2, cf, blackswan)

        return(avg_gt_lr)
        })

      # Combine and sort the scores
      scores <- unlist(scores)
      score_df <- data.frame(score = scores, name = names)
      score_df <- score_df[order(score_df$score, decreasing = TRUE), ]
      score_df[is.na(score_df$score), ]$name <- NA

      # Make a source contamination call
      # Requirements are:
      # 1) A top scoring contaminant exists (non-NA result)
      # 2) Conta made a contamination call (without the use of sources)
      # 3) Contamination fraction is larger than zero
      # 4) Source score exceeds the original conta score by a specified fraction
      call_thr <- 1.1 * avg_maf_lr
      source_call <- !is.na(score_df[1, ]$score) &&
        res$conta_call && cf > 0 &&
        score_df[1, ]$score > call_thr

      # Create base output frame
      out_this <- data.frame(name = names[i],
                             cf = cf,
                             source_call = source_call,
                             avg_maf_lr = avg_maf_lr,
                             best_conf = score_df[1, ]$score,
                             best_sample = score_df[1, ]$name)

      # Add second best result if there is a second result
      if (nrow(score_df) >= 2) {
        out_this$second_conf <- score_df[2, ]$score
        out_this$second_sample <- score_df[2, ]$name
      }

      # Add third best result if there is a third result
      if (nrow(score_df) >= 2) {
        out_this$third_conf <- score_df[3, ]$score
        out_this$third_sample <- score_df[3, ]$name
      }

      # Merge this sample's result with overall results
      out <- rbind(out, out_this)

    }
  }

  out <- format(out, digits = 3, trim = TRUE)

  if (startsWith(base, "s3")) {
    if (requireNamespace("grails3r")) {
      grails3r::write_to_s3(out, out_file, write.table, sep = "\t", row.names = FALSE,
                            col.names = TRUE, quote = FALSE)
    } else {
      s3write_using(out, FUN = write.table, object = out_file, sep = "\t",
                    row.names = FALSE)
    }
  } else {
    write.table(out, out_file, sep = "\t", row.names = FALSE, col.names = TRUE,
                quote = FALSE)
  }
}
