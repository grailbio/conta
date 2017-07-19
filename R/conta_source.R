
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
#' is reported as with the highest likelihood. 2nd and 3rd highest scoring
#' source candidates are also reported.
#'
#' @param base character base directory with a batch of conta results
#' @param out_file character output file
#' 
#' @return none
#'
#' @export
conta_source <- function(base, out_file, cores = 8) {

  options("digits" = 5)
  options("mc.cores" = cores)
  
  # Add slash to the end of base if it wasn' specified
  base <- ifelse(!endsWith(base, "/"), paste(base, "/", sep = ""), base)

  # s3_ls reads both files and paths under a folder, we only need folders
  if (startsWith(base, "s3")) {
    files <- s3_ls(base)$path
    paths <- files[endsWith(files, "/")]
  } else {
    paths <- paste(list.dirs(base, recursive = FALSE, full.names = FALSE),
                  "/", sep = "")
  }

  names <- basename(paths)

  # Read conta results and genotypes, good to cache them at this point
  lres <- mclapply(paste(base, paths, names, ".conta.tsv", sep = ""),
                 read_data_table)

  lgt <- mclapply(paste(base, paths, names, ".gt.tsv", sep = ""),
                read_data_table)

  # Define a data.frame for the output, one line for contaminated sample
  out <- data.frame(name = character(),
                    cf = numeric(),
                    source_call = logical(),
                    avg_maf_lr = numeric(),
                    best_gt_score = numeric(),
                    best_sample = character(),
                    second_gt_score = numeric(),
                    second_sample = character(),
                    third_gt_score = numeric(),
                    third_sample = character())

  # for each conta result file
  for (i in 1:length(lres)) {

    res <- lres[[i]]

    # find source only if it was called contaminated
    if (res$conta_call) {

      # contamination fraction
      cf <- res$cf

      # avg likelihood ratio obtained by using population frequencies
      avg_maf_lr <- res$avg_log_lr

      # read the genotypes and variant reads for this sample
      # we only need hom alleles for source detection
      gt1 <- lgt[[i]][gt != "0/1"]

      # calculate source likelihood scores for each of the other samples
      scores <- mclapply(c(1:length(lres)), function(j) {
        
        gt2 <- lgt[[j]][gt != "0/1"]
        
        conc <- get_genotype_concordance(gt1, gt2)
        
        if (is.na(conc) | conc >= 0.95) {
          return(NA)
        }
        
        avg_gt_lr <- get_source_lr(gt1, gt2, cf)
        return(avg_gt_lr)
        })
      scores <- unlist(scores)

      # Find the best score, its sample as well as second and third scores
      # and their samples.
      # TODO: implement ranked list rather than best, second, third approach
      best <- max(c(scores, 0), na.rm = TRUE)
      best_loc <- which(best == scores)[1]
      second <- max(c(scores[-best_loc], 0), na.rm = TRUE)
      second_loc <- ifelse(second < best,
                           which(second == scores)[1],
                           which(second == scores)[2])
      third <- max(c(scores[-c(best_loc, second_loc)], 0), na.rm = TRUE)
      third_loc <- ifelse(third < best,
                              ifelse(third < second,
                                     which(third == scores)[1],
                                     which(third == scores)[2]),
                              which(third == scores)[3])

      # add the calls to output
      out <- rbind(out, data.frame(name = names[i],
                                   cf = cf,
                                   source_call = best > avg_maf_lr,
                                   avg_maf_lr = avg_maf_lr,
                                   best_gt_score = best,
                                   best_sample = names[best_loc],
                                   second_gt_score = second,
                                   second_sample = names[second_loc],
                                   third_gt_score = third,
                                   third_sample = names[third_loc]))
    }
  }

  out <- format(out, digits = 3)

  if (startsWith(base, "s3")) {
    write_to_s3(out, out_file, write.table, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE)
  } else {
    write.table(out, out_file, sep = "\t", row.names = FALSE, col.names = TRUE,
                quote = FALSE)
  }
}