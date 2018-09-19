#' Calculate error model
#'
#' Calculate error model for ref/ref and alt/alt alleles globally across
#' all positions that should have no trace of contamination.
#'
#' @param dat data.frame containing counts and metrics per SNP and only hom
#'   alleles (0/0 and 1/1). Het alleles should have already been removed.
#' @param save.dir character location to save the results
#' @param default_error_rate if absolutely no errors are observed in the data,
#'     use this number as the error rate
#' @param sample character sample name
#' @param context_mode whether to run in context mode
#'
#' @return data.frame with error model
#'
#' @importFrom utils write.table
#' @export
calculate_error_model <- function(dat, save.dir = NA, sample = NA,
                                  default_error_rate = 1e-5,
                                  context_mode = FALSE) {

  # Prepare error model data frame
  bases <- get_bases(context_mode)
  EE <- data.frame(ref = rep(bases, length(bases)),
                   subs = rep(bases, each = length(bases)),
                   er = 0, reads = 0, denom = 0)
  rownames(EE) <- paste(EE$ref, EE$subs, sep = ">")
  EE$ref <- as.character(EE$ref)

  # Rebuild EE for context from only compatible base changes
  for (i in bases) {
    for (j in bases) {

      # Skip mismatching contexts in the context mode
      if (context_mode && !is_compatible_context(i, j)) {
        subs <- paste(i, j, sep = ">")
        EE <- EE[-which(row.names(EE) == subs), ]
      }
    }
  }

  # Calculate number of reads observed for each substitution
  for (i in bases) {
    for (j in bases) {

      # Skip mismatching contexts in the context mode
      if (context_mode && !is_compatible_context(i, j))
        next

      subs <- paste(i, j, sep = ">")

      # Look at non-SNP base substitutions (major=i and minor!=j) and count the
      # occurence of each event i>j in the context of the genotype
      if (!context_mode) {
        EE[subs, ]$reads <-
          sum(as.matrix(dat[gt == "0/0" & major == i & minor != j, j,
                          with = FALSE]),
            as.matrix(dat[gt == "1/1" & minor == i & major != j, j,
                          with = FALSE]))
      } else {

        # Context reference and snp are already updated based on genotype
        # See update_context() function
        j2 <- substring(j, 2, 2)
        EE[subs, ]$reads <-
          sum(as.matrix(dat[gt != "0/1" & context == i & context_snp != j,
                            j2, with = FALSE]))
      }
    }
  }

  # Calculate denominator and error rate for each substitution
  for (i in bases) {
    for (j in bases) {

      # Skip mismatching contexts in the context mode
      if (context_mode && !is_compatible_context(i, j))
        next

      subs <- paste(i, j, sep = ">")

      # Denominator when calculating error rate on SNP position
      # is over-counted, therefore it needs to be corrected.
      # For example, to calculate C>T error rate, we count the number of
      # times the substitution C>T is observed on C>A and C>G SNPs. We cannot
      # use C>T SNPs for this purpose due to the possibility of contamination.
      # But, when calculating the denominator (how many times the C base is
      # observed), we calculate all C>C events that come from all 3 types of
      # SNPs described above. Thus we need to reduce this denominator to two
      # thirds. One caveat is that, Non-SNP alleles (maf = 0), should not
      # receive this correction since they are not supposed to have any
      # contamination and can be used to look at all types of error rates.
      # In other terms, denom_factor will be 0.66 if all bases in dat are
      # from dbSNP, otherwise it will be a larger number closer to 1.
      denom_factor <- 2/3 * dat[, mean(maf > 0)] + dat[, mean(maf == 0)]

      # Denominator in error model is all the times ref base was the original
      # base. For example, denom for A>A, A>T, A>G, and A>C are all the same,
      # and some of the reads representing these cases
      EE[subs, ]$denom <- round(denom_factor * sum(EE[EE$ref == i, ]$reads))

      # If no errors are observed, set error rate as NA if there are less than
      # 1000 reference observations, otherwise set it as 1 / observations.
      # If any number of errors are observed, set it as errors / observations.
      EE[subs, ]$er <- ifelse(EE[subs, ]$reads == 0,
                              ifelse(EE[subs, ]$denom < 1000, NA,
                                     1 / EE[subs, ]$denom),
                              EE[subs, ]$reads / EE[subs, ]$denom)
    }
  }

  # Remove ref->ref rates, they were only calculated for the denominator
  EE <- EE[EE$ref != EE$subs, ]

  # If there are missing error rates, set them as global average
  # First, handle edge cases if there are no errors
  if (sum(is.na(EE$er)) == nrow(EE) | mean(EE$er, na.rm = TRUE) == 0) {
    EE$er <- default_error_rate
  } else if (sum(is.na(EE$er)) > 0) {
    EE[is.na(EE$er), ]$er <- mean(EE$er, na.rm = TRUE)
  }

  # Write error model to file if save.dir and sample are specified
  if (!is.na(save.dir) && !is.na(sample)) {
    error.file <- file.path(save.dir, paste(sample, "error.tsv", sep = "."))
    write.table(format(EE, digits = 5, trim = TRUE), file = error.file,
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  return(EE)
}

#' Add error rates for each base based on ref/ref and alt/alt error rates
#'
#' @param dat data.table containing counts and metrics per SNP
#' @param EE data.frame errror model
#' @param context_mode
#'
#' @return data.frame with error model
#'
#' @export
add_error_rates <- function(dat, EE, context_mode = FALSE) {

  # Substitution type
  dat[, et := ""]
  if (context_mode) {
    dat[, et := paste(context, ">", context_snp, sep = "")]
  } else {
    dat[gt == "0/0", et := paste(major, ">", minor, sep = "")]
    dat[gt == "0/1", et := paste(major, ">", minor, sep = "")]
    dat[gt == "1/1", et := paste(minor, ">", major, sep = "")]
  }

  # Error rate corresponding to that substitution type
  dat[, er := EE[dat$et, ]$er]

  return(dat)
}
