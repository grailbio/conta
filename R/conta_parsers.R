#' Read a data.table from given path
#'
#' Check if file name starts with s3 or not, use appropriate function to read.
#' If the file does not exist, exit with an error message if needed.
#'
#'
#' @param file character file name to read
#' @param header logical whether to read the header
#' @param sep character the separator
#' @param stop_if_missing logical whether to stop execution if file is missing
#' @return data.table
#'
#' @export
read_data_table <- function(file, header = TRUE, sep = "\t",
                            stop_if_missing = FALSE) {
  if (startsWith(file, "s3://")) {
    if (s3_file_exists(file)) {
      if (endsWith(file, ".gz")) {
        #TODO
        stop(paste("Reading zipped file from s3 not supported.", file))
      } else {
        return(read_from_s3(file, fread, header = header,
                          stringsAsFactors = FALSE, sep = sep))
      }
    }
  }
  else {
    if (file.exists(file)) {
      if (endsWith(file, ".gz")) {
        return(fread(sprintf("zcat %s", file), header = header,
                     stringsAsFactors = FALSE, sep = sep))
      } else {
        return(fread(file, header = header, stringsAsFactors = FALSE,
                     sep = sep))
      }
    }
  }

  if (stop_if_missing) {
    stop(paste("Input file not found:", file))
  } else {
    warning(paste("Input file not found:", file))
    return(NULL)
  }
}

#' Return true if file exists on s3
#'
#' @param file character file name to check existence
#' @export
s3_file_exists <- function(file) {
  return(suppressWarnings(nrow(s3_ls(file)) > 0))
}

#' Read a counts tsv file and prep it
#'
#' Supports reading file from s3 (to a temp file which is later removed)
#'
#' @param file character tsv file containing SNP info for the sample
#' @return data.table containing counts and metrics per SNP
#'
#' @export
read_and_prep <- function(file) {

  dat <- read_data_table(file, stop_if_missing = TRUE)

  # Find minor allele
  dat[, c("a1", "a2") := tstrsplit(alleles, "/", fixed = TRUE)]
  dat$minor <- ifelse(dat$a1 == dat$major, dat$a2, dat$a1)
  dat[, c("a1", "a2") := NULL]

  # Add some other fields
  dat <- ratio_and_counts(dat)

  return(dat)
}

#' Add in basic counts and ratios
#'
#' @param dat data.frame containing counts and metrics per SNP
#' @return data.table containing counts and metrics per SNP
#'
#' @export
ratio_and_counts <- function(dat) {
  
  # Recalculate depth
  dat[, depth := A + T + G + C]
  
  # TODO: implement with data.matrix and match function
  dat$major_count <- ifelse(dat$major == "A", dat$A,
                            ifelse(dat$major == "T", dat$T,
                                   ifelse(dat$major == "G", dat$G,
                                          ifelse(dat$major == "C", dat$C, 0))))
  
  dat$major_ratio <- dat$major_count/dat$depth
  
  dat$minor_count <- ifelse(dat$minor == "A", dat$A,
                            ifelse(dat$minor == "T", dat$T,
                                   ifelse(dat$minor == "G", dat$G,
                                          ifelse(dat$minor == "C", dat$C, 0))))
  
  dat$minor_ratio <- dat$minor_count/dat$depth
  
  dat$other_count <- dat$depth - dat$major_count - dat$minor_count
  dat$other_ratio <- round(1 - dat$major_ratio - dat$minor_ratio, 4)
  
  return(dat)
}

#' Add in more useful columns, annotation and filter sites.
#'
#' @param dat data.table containing counts and metrics per SNP
#' @param het_limit limit fraction to call a heterozygote site.
#      SNPs with variant frequency above this value and below 1 minus this value
#      will be called as heterozygotes. Rest of the SNPs will be called
#      homozygotes.
#' @param min_other_ratio If the ratio of highest depth non-SNP allele is at
#'     as high as this number, the SNP will be filtered. These SNPs either have
#'     unusually high error rate (or mutation), or has a hidden (multi-allelic)
#'     SNP. In either case, this position should be filtered from contamination
#'     analysis.
#' @param min_depth if a SNP has depth less than this value, it will be
#'     filtered. This metric is mainly for WGS where overall coverage is low.
#' @param max_sd_depth maximum number of standard deviations the depth
#'     of a given SNP is less than the mean depth across positions. If a
#'     position has depth less than that, it will be filtered. This metric is
#'     mainly for positions with high depth such as targeted or exome sequencing
#'     experiments.
#' @return data.frame containing counts and metrics per SNP
#'
#' @export
annotate_and_filter <- function(dat, het_limit = 0.25, min_other_ratio = 0.15,
                                min_depth = 5, max_sd_depth = 3) {

  # Genotype calls
  dat$gt <- ifelse(dat$minor_ratio < het_limit, "0/0",
                   ifelse(dat$major_ratio < het_limit, "1/1", "0/1"))

  # Noise (or contamination) levels, flipped for 1/1 alt alleles
  dat$vfn <- ifelse(dat$gt == "1/1", dat$major_ratio,
                    ifelse(dat$gt == "0/0", dat$minor_ratio, NA))

  # Noise number of reads
  dat$vr <- (dat[, vfn]) * (dat[, depth])

  # Contamination probability based on minor allele frequency
  dat$cp <- ifelse(dat$gt == "1/1", (1 - (dat$maf) ^ 2),
                   ifelse(dat$gt == "0/0", 1 - (1 - dat$maf) ^ 2, NA))

  # Remove otherAllele strange cases
  # This is some type of artifact that is seen in the tsv file
  # It may be due to alignments or tsv generation
  dat <- dat[ !(other_ratio >= min_other_ratio), ]

  # Remove low depth
  dat <- dat[ !is.na(depth) & depth > min_depth &
              depth > (mean(depth) - max_sd_depth * sd(depth)), ]
}

#' Get the fraction of bases covered by at least 1 read on Y chr.
#'
#' Input is a bedtools output where the first column shows the chromosome,
#' and each row stores statistics for a bin of certain size (e.g. 100,000 bp).
#' 7th column in this file represents the percentage of bases on this bin
#' covered by at least 1 read.
#'
#' @param bin_file character input tsv file name
#' @return Y chr counts
#'
#' @export
count_ychr <- function(bin_file) {
  if (is.na(bin_file) | bin_file == "") return(NA)
  bins <- read_data_table(bin_file, header = FALSE)
  return(ifelse(is.null(bins), NA, bins[V1 == "chrY", mean(V7)]))
}

#' Get the final stein score from CNV qc file
#'
#' @param cnv_file character input tsv file name
#' @return Cnv z-score
#'
#' @export
get_final_stein <- function(cnv_file) {
  if (is.na(cnv_file) | cnv_file == "") return(NA)
  cnv <- read_data_table(cnv_file, header = T, sep = "\t", stop_if_missing = F)
  return(ifelse(is.null(cnv), NA, cnv[keys == "final_stein", ]$values))
}

#' Get the final mapd score from CNV qc file
#'
#' @param cnv_file character input tsv file name
#' @return Cnv z-score
#'
#' @export
get_final_mapd <- function(cnv_file) {
  if (is.na(cnv_file) | cnv_file == "") return(NA)
  cnv <- read_data_table(cnv_file)
  return(ifelse(is.null(cnv), NA, cnv[keys == "final_mapd", ]$values))
}
