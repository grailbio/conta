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
#'
#' @return data.frame with error model
#'
#' @importFrom utils write.table
#' @export
calculate_error_model <- function(dat, save.dir = NA, sample = NA,
                                  default_error_rate = 1e-6) {

  # Prepare error model data frame
  bases <- get_bases()
  EE <- data.frame(ref = rep(bases, length(bases)),
                   subs = rep(bases, each = length(bases)),
                   er = 0, reads = 0, denom = 0)
  rownames(EE) <- paste(EE$ref, EE$subs, sep = ">")

  # Calculate number of reads observed for each substitution
  for (i in bases) {
    for (j in bases) {
      subs <- paste(i, j, sep = ">")

      # Look at non-SNP positions (major=i and minor!=j) and count the occurence
      # of each event i>j in the context of the genotype
      EE[subs, ]$reads <-
        sum(as.matrix(dat[gt == "0/0" & major == i & minor != j, j,
                          with = FALSE]),
            as.matrix(dat[gt == "1/1" & minor == i & major != j, j,
                          with = FALSE]))
    }
  }

  # Calculate denominator and error rate for each substitution
  for (i in bases) {
    for (j in bases) {
      subs <- paste(i, j, sep = ">")

      # Since SNP alleles are masked, error rates are underestimated.
      # For example C>A error rate is calculated only from C>T and
      # C>G SNPs. This can be remedied by discounting two thirds of denominator.
      # Non-SNP alleles (maf = 0), should not receive this correction.
      denom_factor <- 2/3 * dat[, mean(maf > 0)] + dat[, mean(maf == 0)]
      EE[subs, ]$denom <- round(denom_factor * sum(EE[EE$ref == i, ]$reads))

      # If no errors are observed, set error rate as NA if there are less than
      # 1000 reference observations, otherwise set is as 1 / observations.
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
#'
#' @return data.frame with error model
#'
#' @export
add_error_rates <- function(dat, EE) {

  # Substitution type
  dat[, et := ""]
  dat[gt == "0/0", et := paste(major, ">", minor, sep = "")]
  dat[gt == "0/1", et := paste(major, ">", minor, sep = "")]
  dat[gt == "1/1", et := paste(minor, ">", major, sep = "")]

  # Error rate corresponding to that substitution type
  dat[, er := EE[dat$et, ]$er]

  return(dat)
}
