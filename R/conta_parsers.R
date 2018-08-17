# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Read a data.table from given path
#'
#' Check if file name starts with s3 or not, use appropriate function to read.
#' If the file does not exist, exit with an error message if needed.
#'
#'
#' @param file file name to read
#' @param header logical whether to read the header
#' @param sep the separator
#' @param stop_if_missing logical whether to stop execution if file is missing
#' @param ... additional arguments
#' @return data.table
#'
#' @export
read_data_table <- function(file, header = TRUE, sep = "\t",
                            stop_if_missing = FALSE, ...) {

  if (startsWith(file, "s3://")) {
    if (requireNamespace("grails3r")) {
      if (grails3r::s3_file_exists(file)) {
        if (endsWith(file, ".gz")) {
          stop(paste("Reading zipped file from s3 not supported.", file))
      } else {
        return(read_from_s3(file, fread, header = header,
                          stringsAsFactors = FALSE, sep = sep, ...))
      }
    }} else {
        # Use aws.s3 package to read if grails3r is not found
        creds <- aws.signature::locate_credentials()
        if (as.logical(do.call(aws.s3::head_object, c(file, creds)))) {
          if (endsWith(file, ".gz")) {
            stop(paste("Reading zipped file from s3 not supported.", file))
            } else {
              return(s3read_using(fread, object = file, header = header,
                                  stringsAsFactors = FALSE, sep = sep, ...))
              }
          }
        }
    } else {
    if (file.exists(file)) {
      if (endsWith(file, ".gz")) {
        return(fread(sprintf("zcat %s", file), header = header,
                     stringsAsFactors = FALSE, sep = sep, ...))
      } else {
        return(fread(file, header = header, stringsAsFactors = FALSE,
                     sep = sep, ...))
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

#' Read a counts tsv file and prep it
#'
#' @param file tsv file containing SNP info for the sample
#' @param file_rev tsv file containing SNP counts for reverse strand
#' @return data.table containing counts and metrics per SNP
#'
#' @export
read_and_prep <- function(file, file_rev = NA) {

  dat <- read_data_table(file, stop_if_missing = TRUE)

  # If reverse strand counts are provided, merge it with the first
  # set of counts which is supposed to come from the positive strand
  if (!is.na(file_rev)) {
    dat2 <- read_data_table(file_rev, stop_if_missing = TRUE)
    dat <- combine_counts(dat, dat2)
  }

  # Return if data.table is empty
  if (nrow(dat) == 0)
    return(dat)

  # Find minor allele
  dat[, c("a1", "a2") := data.table::tstrsplit(alleles, "/", fixed = TRUE)]
  dat$minor <- ifelse(dat$a1 == dat$major, dat$a2, dat$a1)
  dat[, c("a1", "a2") := NULL]

  # Add some other fields
  dat <- ratio_and_counts(dat)

  return(dat)
}

#' Combine two count data tables
#'
#' @param dat data.table counts from first file
#' @param dat2 data.table counts from second file
#' @return data.table containing merged counts
#'
#' @export
combine_counts <- function(dat, dat2) {

  dat <- suppressWarnings(set_numeric_chrs(dat))
  dat2 <- suppressWarnings(set_numeric_chrs(dat2))

  data.table::setkey(dat, chrom_int, pos, ref, major, alleles, rsid, maf, chrom)
  data.table::setkey(dat2, chrom_int, pos, ref, major, alleles, rsid, maf, chrom)
  datm <- merge(dat, dat2, sort = TRUE, all = TRUE)
  datm$A <- rowSums(datm[, c("A.x", "A.y")], na.rm = TRUE)
  datm$T <- rowSums(datm[, c("T.x", "T.y")], na.rm = TRUE)
  datm$G <- rowSums(datm[, c("G.x", "G.y")], na.rm = TRUE)
  datm$C <- rowSums(datm[, c("C.x", "C.y")], na.rm = TRUE)
  datm$N <- rowSums(datm[, c("N.x", "N.y")], na.rm = TRUE)
  return(datm[, -c("A.x", "A.y", "T.x", "T.y", "G.x", "G.y",
                    "C.x", "C.y", "N.x", "N.y", "chrom_int")])
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

  dat$major_ratio <- dat$major_count / dat$depth

  dat$minor_count <- ifelse(dat$minor == "A", dat$A,
                            ifelse(dat$minor == "T", dat$T,
                                   ifelse(dat$minor == "G", dat$G,
                                          ifelse(dat$minor == "C", dat$C, 0))))

  dat$minor_ratio <- dat$minor_count / dat$depth

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
#' @param max_depth if a SNP has depth more than this value, it will be
#'     filtered. This metric is mainly for RNA where some genes have extreme
#'     coverage.
#' @param out_frac remove SNPs that are outliers to depth distribution.
#'     This filter is applied post min_depth and max_depth filters to further
#'     remove any outliers that also tend to generate unexpected noise.
#' @return data.frame containing counts and metrics per SNP
#'
#' @importFrom stats quantile
#' @export
annotate_and_filter <- function(dat, het_limit = 0.25, min_other_ratio = 0.15,
                                min_depth = 5, max_depth = 10000,
                                out_frac = 0) {

  # Return if data.table is empty
  if (nrow(dat) == 0)
    return(dat)

  # Remove non-ATGC alleles
  dat <- dat[major %in% c(get_bases(), "N")]
  dat <- dat[minor %in% c(get_bases(), "N")]

  # Genotype calls
  dat$gt <- ifelse(dat$minor_ratio < het_limit, "0/0",
                   ifelse(dat$major_ratio < het_limit, "1/1", "0/1"))

  # Noise (or contamination) levels, flipped for 1/1 alt alleles
  dat$vfn <- ifelse(dat$gt == "1/1", dat$major_ratio,
                    ifelse(dat$gt == "0/0", dat$minor_ratio,
                           ifelse(dat$gt == "0/1", pmin(dat$major_ratio,
                                                        dat$minor_ratio), NA)))

  # Noise number of reads
  dat$vr <- (dat[, vfn]) * (dat[, depth])

  # Contamination probability based on minor allele frequency
  dat$cp <- ifelse(dat$gt == "1/1", (1 - (dat$maf) ^ 2),
                   ifelse(dat$gt == "0/0", 1 - (1 - dat$maf) ^ 2,
                          ifelse(dat$gt == "0/1",
                                 1 - 2 * dat$maf * (1 - dat$maf), NA)))

  # Remove otherAllele strange cases
  # This is some type of artifact that is seen in the tsv file
  # It may be due to alignments or tsv generation
  dat <- dat[ !(other_ratio >= min_other_ratio), ]

  # Remove low and high depth
  dat <- dat[ !is.na(depth) & depth >= min_depth & depth <= max_depth, ]

  # Remove depth outliers
  dat <- dat[depth <= quantile(depth, 1 - out_frac)
             & depth >= quantile(depth, out_frac), ]

  # Add chunks
  dat <- add_chunks(dat)
}

#' Get the stats for a given chromosome from long biometrics file, including:
#'  - Mapq60 read count
#'  - Normalized mapq60 read count (normalized by total reads for this sample)
#'  - Fraction of bases covered by at least 1 read
#'
#' @param biometrics_file input tsv file name
#' @param chr_name name of the chromosome to be counted
#'
#' @return list of chromosome metrics
#'
#' @export
chr_stats <- function(biometrics_file, chr_name) {

  # Read file
  if (file.exists(biometrics_file)) {
    k1 <- read_data_table(biometrics_file, sep = "\t")
  }
  else {
    return(list( count = NA, normalized_count = NA, fraction_covered = NA))
  }

  # Calculate total reads
  mapq60_reads <- k1[(key == "mapq60_count") &
                       (metric_key == "fragment_counts")]$value

  # Cast metrics of interest to long tables
  m11 <- k1[metric_key == "fraction_chr_covered", ]
  d11 <- dcast(m11, group_id ~ key, value.var = "value")
  m12 <- k1[metric_key == "mapq60_fragment_counts_per_chr"]
  d12 <- dcast(m12, group_id ~ key, value.var = "value")

  # Normalize counts
  d12$normalized_count <- as.numeric(d12$count) / as.numeric(mapq60_reads)

  # Return list with NA elements if the specified chr does not exist
  if (nrow(d12[ chr == chr_name ]) == 0) {
    return(list( count = NA, normalized_count = NA, fraction_covered = NA))
  }

  # Return a list of values
  return(list( count = as.numeric(d12[chr == chr_name, count]),
               normalized_count = d12[chr == chr_name, normalized_count],
               fraction_covered = as.numeric(d11[chr == chr_name,
                                                 fraction_covered])))
}

#' Read a vcf file as a data table
#'
#' @param vcf_file name of the input file
#' @param n max number of header lines to skip
#'
#' @return data.table containing vcf data
#'
#' @export
read_vcf_dt <- function(vcf_file, n = 100000) {
  if (startsWith(vcf_file, "s3://")) {
    if (requireNamespace("grails3r"))
      lines <- grails3r::read_from_s3(vcf_file, readLines, n = n)
    else
      lines <- aws.s3::s3read_using(readLines, object = vcf_file, n = n)
  } else {
    lines <- readLines(vcf_file, n)
  }
  skip_lines <- 0
  for (i in 1:length(lines)) {
    if (startsWith(lines[i], "##"))
      next
    else {
      skip_lines <- i - 1
      break
    }
  }
  vcf_dt <- read_data_table(vcf_file, header = TRUE, sep = "\t",
                            stop_if_missing = FALSE, skip = skip_lines)

  caf <- parse_field(vcf_dt$INFO, "CAF")
  vcf_dt$maf <- as.numeric(substring(caf, regexpr(",", caf) + 1,
                                     nchar(caf) - 2))

  return(vcf_dt)
}

#' Retrieve a specified field from a vcf info string
#'
#' @param info vcf info string
#' @param field a specific field from info
#' @return string field
#'
#' @export
parse_field <- function(info, field) {
  start_loc <- regexpr("CAF=", info)
  end_loc <- start_loc +
    regexpr(";", substring(info, start_loc, nchar(info))) - 1
  caf_string <- substring(info, start_loc, end_loc)
}

#' Partition each chromosome to bins and tag each SNP with its bin
#'
#' @param dat data.table containing SNPs
#' @param max_portions number of portions to split each chromosome
#' @return data.table containing SNPs with tagged chunks
#'
#' @export
add_chunks <- function(dat, max_portions = 10) {

  if (nrow(dat) == 0) return(data.table())

  portions <- max_portions # partions per chr
  dat[, chunk := 0]
  for (j in dat[, unique(chrom)]) {
    bin_size <- ceiling(as.numeric(dat[chrom == j, .(.N / portions)]))
    dat[chrom == j, chunk := .(ceiling(.I / bin_size))]
  }
  dat$chrom <- factor(dat$chrom, levels = unique(dat$chrom))

  return(dat)
}

#' Read simulation results from provided path
#'
#' Reads conta simulation results from a local or S3 path. Checks
#' that the required columns are present and properly formatted,
#' and returns the simulation data as a data frame.
#'
#' @param file file name to read
#' @return data.table
#'
#' @export
read_sim_results <- function(file) {
  raw_data <- read_data_table(file,
                              stop_if_missing = TRUE)

  req_cols <- c("sample", "avg_log_lr")

  diff <- setdiff(req_cols, colnames(raw_data))
  if (length(diff) != 0) {
    stop(paste0("ERROR: Simulation results file is missing required columns.\n",
                "file: ", file, "\nmissing cols: ", diff))
  }

  sample_splits <- strsplit(raw_data$sample, "_")

  sample_names <- sapply(sample_splits,
                         function(x) x[[1]])
  conta_level <- sapply(sample_splits,
                        function(x) as.numeric(x[[2]]))
  if (any(is.na(conta_level))) {
    stop(paste0("ERROR: Expect simulation `sample` column to be of the form ",
                "{sample_id}_{conta_level}"))
  }

  raw_data$sample <- sample_names
  raw_data$conta_level <- conta_level

  return(raw_data)
}
