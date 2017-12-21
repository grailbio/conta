#'
#' Run swap tool
#'
#' Compare genotypes across specified sample pairings (specific to GRAIL CCGA
#' data). Both concordance and distance metrics are reported along with other
#' informative features such as coverage, CNV indicators and contamination.
#'
#' @param pairing_table_name name of the mapping table
#' @param out_table_name name of the results table to output
#' @param shiny_loc location of the shiny database to retrieve metrics
#' @param dbsnp_targeted targeted dbsnp vcf file same as used for targeted seq
#' @param corrections file with swap corrections that were confirmed
#' @param randomize generate random pairs instead of expected ones
#' @param cutoff cutoff to make a cfdna targeted vs wgs swap call
#' @param max_conta if contamination is above this level, do not call swap
#' @param signatures location of file containing signature SNPs
#'
#' @export
swap <- function(pairing_table_name, out_table_name, shiny_loc, dbsnp_targeted,
                 corrections, randomize = FALSE,
                 cutoff = 0.7, max_conta = 0.2, signatures) {

  if (randomize) {
    out_table_name <- paste(gsub(".csv", "", out_table_name),
                           ".random.csv", sep = "")
  }

  # Read mapping table. This table is an input that stores the pairings between
  # cfDNA and gDNA samples, along with a fake pair ID (DID) and StudyBatch
  pairing_table <- read_data_table(pairing_table_name, showProgress = FALSE)

  # Read a set of corrections from table
  corrections_table <- read_data_table(corrections)

  # Read signature SNPs
  signature_SNPs <- read_data_table(signatures)

  # Results are stored in a match table. Existing results are not re-calculated.
  # To re-calculate the output file needs to be removed.
  mt <- init_match_table(out_table_name, pairing_table)

  # Read targeted SNPs file. This is a list of SNPs in dbSNP vcf format that
  # cover the targeted panel.
  snps <- read_vcf_dt(dbsnp_targeted)

  # SeqSampleIDs and common metrics for each sample are stored in shiny data.
  shiny_dt <- read_shiny_data(shiny_loc)

  # Run swap concordance for each sample if it was not already calculated
  for (i in 1:nrow(mt)) {

    # This part is input dependent, it makes a table with locations and metrics
    set_locations_and_metrics(shiny_dt, mt, i, randomize)

    # Load data.tables for each of the samples (data.tables are loaded only if
    # necessary. If all of its associated scores were already calculated, it
    # will not be loaded.
    d <- get_data_tables(mt, i, snps)

    # Calculate swap concordance for each of the samples
    calculate_concordance(mt, i, d, cutoff, signature_SNPs)

    # Do not call swap if contamination is above a threshold
    if (mt[i, swap == "true"] & ((!is.na(mt[i, ctar_max_lh])
                                  && mt[i, ctar_max_lh] >= max_conta) |
                                 (!is.na(mt[i, gtar_max_lh])
                                  && mt[i, gtar_max_lh] >= max_conta)))
      mt[i, swap := "contaminated"]

    # Do not call swap if it is previously detected and corrected
    if (mt[i, swap == "true"] &
        (mt[i, gtar_SeqSampleID] %in% corrections_table$orig_id ||
         mt[i, gwgs_SeqSampleID] %in% corrections_table$orig_id ||
         mt[i, ctar_SeqSampleID] %in% corrections_table$orig_id ||
         mt[i, cwgs_SeqSampleID] %in% corrections_table$orig_id ||
         mt[i, wgbs_SeqSampleID] %in% corrections_table$orig_id))
      mt[i, swap := "corrected"]

    # Update results
    write_results(d, mt, out_table_name)

  }
}


#' Create results table
#'
#' If it is first time: copy mapping table and add new fields as NA
#' If it already exists, match it with mapping table to add any new
#' samples that were added since last time.
#'
#' @param out_table_name location of output file
#' @param pairing_table data.table for mapping file
#'
#' @return initialized result table
#'
#' @export
init_match_table <- function(out_table_name, pairing_table) {

  mt <- read_data_table(out_table_name, showProgress = FALSE)

  # If first time creating the output, copy it from match original table,
  # and add new fields.
  if (is.null(mt)) {
    nl <- nrow(pairing_table)
    mt <- data.table(DID = pairing_table$DID,
                     StudyBatch = pairing_table$StudyBatch,
                     gdna_PID = pairing_table$gdna_PID,
                     cfdna_PID = pairing_table$cfdna_PID,
                     wgbs_PID = pairing_table$wgbs_PID,
                     swap = character(length = nl),
                     ctar_SeqSampleID = character(length = nl),
                     cwgs_SeqSampleID = character(length = nl),
                     gtar_SeqSampleID = character(length = nl),
                     gwgs_SeqSampleID = character(length = nl),
                     wgbs_SeqSampleID = character(length = nl),
                     ctar_mtc = numeric(length = nl),
                     cwgs_cov = numeric(length = nl),
                     gtar_mtc = numeric(length = nl),
                     gwgs_cov = numeric(length = nl),
                     wgbs_cov = numeric(length = nl),
                     ctar_num_amplified_genes = double(length = nl),
                     cwgs_final_stein = double(length = nl),
                     cwgs_final_mapd = double(length = nl),
                     gtar_num_amplified_genes = double(length = nl),
                     gwgs_final_stein = double(length = nl),
                     gwgs_final_mapd = double(length = nl),
                     wgbs_final_stein = double(length = nl),
                     wgbs_final_mapd = double(length = nl),
                     conc_ctar_cwgs = double(length = nl),
                     conc_gtar_gwgs = double(length = nl),
                     conc_ctar_gtar = double(length = nl),
                     conc_cwgs_gwgs = double(length = nl),
                     conc_cwgs_wgbs = double(length = nl),
                     conc_gwgs_wgbs = double(length = nl),
                     ctar_conta_call = logical(length = nl),
                     ctar_max_lh = double(length = nl),
                     ctar_lh_diff = double(length = nl),
                     gtar_conta_call = logical(length = nl),
                     gtar_max_lh = double(length = nl),
                     gtar_lh_diff = double(length = nl),
                     ctar_signature = character(length = nl),
                     gtar_signature = character(length = nl),
                     cwgs_signature = character(length = nl),
                     gwgs_signature = character(length = nl),
                     wgbs_signature = character(length = nl),
                     ctar_loc = character(length = nl),
                     gtar_loc = character(length = nl),
                     cwgs_loc = character(length = nl),
                     gwgs_loc = character(length = nl),
                     wgbs_loc = character(length = nl))
    mt$swap <- "false"
  } else {
    mt <- merge(pairing_table, mt,
                         by = c("DID", "gdna_PID", "cfdna_PID", "StudyBatch"),
                         all.x = TRUE, all.y = TRUE)
  }
  return(mt)
}

#' Return pre-calculated shiny record
#'
#' @param shiny_dt shiny data table
#' @param tid trimmed patient id with last digit missing
#' @param aid CFDNA_CLINICAL or GDNA_CLINICAL
#' @param atype wgs or targeted
#' @return matching shiny record
#' @export
get_shiny_record <- function(shiny_dt, tid, aid, atype) {
  sds <- shiny_dt[(assay_id == aid) & (assay == atype)]
  return(sds[which(tid == substring(sds$SeqSampleID, 1, 6))][1])
}

#' Get conta location from a shiny record
#'
#' @param record record needs to have fields Results.Location and SeqSampleID
#' @param targeted whether the file set is targeted (vfn) or wgs (gt)
#' @return conta file with genotype and allele frequencies
#'
#' @export
get_conta_file_loc <- function(record, targeted) {
  if (is.null(record) | is.na(record$Results.Location)) {
    return(NA)
  } else if (targeted) {
    manifest_loc <- paste(record$Result.Location, "/MANIFEST", sep = "")
    if (s3_file_exists(manifest_loc)) {
      manifest <- read_data_table(manifest_loc)
      return(paste(record$Results.Location,
                   manifest[FileType == "conta_vfn_table", MatchPattern],
                   sep = "/"))
    } else {
      return(paste(record$Results.Location, "/report/conta/",
                   record$SeqSampleID, ".vfn.csv", sep = ""))
    }
  } else {
    return(paste(record$Results.Location, "/contamination/",
                 record$SeqSampleID, ".gt.tsv", sep = ""))
  }
}


#' Numerize columns from shiny data if they are set as string
#'
#' @param shiny_dt data.table
#' @param cols columns to numerize
#'
#' @export
numerize_columns <- function(shiny_dt, cols) {

  # Find their location across the record
  locs <- match(cols, colnames(shiny_dt))

  # Set them as numeric
  for (loc in seq_along(locs)) {
    set(shiny_dt, NULL, locs[loc], as.numeric(shiny_dt[[locs[loc]]]))
  }
}

#' Load conta files if necessary
#'
#' For a given record, mt and snps, read vfn or gt file, and return it.
#'
#' @param conta_loc location of conta vfn or gt file
#' @param snps targeted dbsnp data.table
#' @param targeted whether to retrieve vfn or gt file
#' @return data.table with conta genotypes
#'
#' @export
load_conta_file <- function(conta_loc, snps = NULL, targeted = TRUE) {

  if ((!is.null(snps)) & !is.na(conta_loc) &
      (file.exists(conta_loc) || s3_file_exists(conta_loc))) {

    # Read targeted genotypes and rsid from vfn file format
    conta_loc_dt <- read_data_table(conta_loc, sep = ",", showProgress = FALSE)

    # Add rsid. Targeted conta output doesn't have this info
    if (!is.null(conta_loc_dt))
      conta_loc_dt[, rsid := snps$ID]

  } else if (!is.na(conta_loc) &
             (file.exists(conta_loc) || s3_file_exists(conta_loc))) {

    # Read conta genotypes and rsid from gt file format
    conta_loc_dt <- read_data_table(conta_loc, showProgress = FALSE)

    # Vf is the variant allele frequency, needs to be set of wgs files
    conta_loc_dt[, vf := vr / dp]

  }
  else {

    conta_loc_dt <- NULL

  }

  return(conta_loc_dt)
}

#' Calculate concordance between two samples' genotypes and make a swap call
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param dt1 genotype set 1
#' @param dt2 genotype set 2
#' @param col1 column name 1 for results
#' @param conc_cutoff concordance cutoff to make a swap call
#' @param min_dp minimum average depth to include a sample for concordance test
#' @param min_maf minimum population frequency for a SNP to be included
#'
#' @export
swap_concordance <- function(mt, i, dt1, dt2, col1, conc_cutoff,
                             min_dp = 10, min_maf = 0.25) {

  # Check if any of the datasets are null, return original table
  if (is.null(dt1) | is.null(dt2))
    return()

  # Do not calculate concordance if depth is less than 10 for any of the samples
  if (mean(dt1$dp, na.rm = TRUE) < min_dp | mean(dt2$dp, na.rm = TRUE) < min_dp)
    return()

  # Calculate concordance
  conc <- genotype_concordance(dt1, dt2, min_maf)

  # Write concordance to results table
  mt[i, eval(col1) := conc]

  # Call swap if the concordance is less than cut-off
  if (!is.na(conc) & (conc <= conc_cutoff))
    mt[i, swap := "true"]
}

#' Calculate concordance between two samples' genotypes
#'
#' @param dt1 genotype set 1
#' @param dt2 genotype set 2
#' @param min_maf minimum population frequency to include a SNP
#'
#' @export
genotype_concordance <- function(dt1, dt2, min_maf = 0.25) {

  # If maf is specified (only targeted conta results have it), filter them
  maf_field1 <- which(colnames(dt1) %in% "maf")
  if (length(maf_field1) > 0)
    dt1 <- dt1[maf >= min_maf, ]

  # Merge the two data tables to put each SNP in a single row
  m1 <- merge(dt1, dt2, by = "rsid")

  # Calculate genotype concordance which is the fraction of matching genotypes
  return(m1[, mean(gt.x == gt.y, na.rm = TRUE)])
}

#' Read shiny data to a table
#' @param shiny_loc s3 or local path to shiny table
#' @return data.table
#'
#' @export
read_shiny_data <- function(shiny_loc) {
  shiny_dt <- read_data_table(shiny_loc, sep = ",", showProgress = FALSE)
  shiny_dt <- shiny_dt[ !is.na(shiny_dt$SeqSampleID), ]
  shiny_dt <- shiny_dt[ SampleType == "sample", ]

  # Set following columns as numeric
  cols <- c("collapsed_fragment_mean_target_coverage", "num_amplified_genes",
            "DEDUP_FILTERED_BAM_coverage", "final_stein", "max_lh", "lh_diff",
            "final_mapd")
  numerize_columns(shiny_dt, cols)

  return(shiny_dt)
}

#' Get signature for a given set of SNPs from a conta vfn or gt file
#'
#' @param conta_dt targeted vfn data table
#' @param signature_SNPs table containing rsid, ref and alt of signature SNPs
#' @return data.table with conta genotypes
#'
#' @export
get_signature <- function(conta_dt, signature_SNPs) {

  m1 <- merge(conta_dt, signature_SNPs, by = "rsid", all.y = TRUE)
  return(paste(ifelse(is.na(m1$gt), "N",
                      ifelse(m1$gt == "0/0", m1$ref,
                             ifelse(m1$gt == "0/1", get_IUPAC(m1$ref, m1$alt),
                                    m1$alt))), collapse = ""))
}

#' Get ambiguity code for a given combination
#'
#' TODO: support full IUPAC table:
#' https://droog.gs.washington.edu/parc/images/iupac.html
#'
#' @param base1 first base
#' @param base2 second base
#' @return iupac code for the combination of bases
#'
#' @export
get_IUPAC <- function(base1, base2) {
  iupac_code <- dplyr::case_when(
    (base1 == base2) ~ as.character(base1),
    (base1 == "A" & base2 == "T") ~ "W",
    (base1 == "T" & base2 == "A") ~ "W",
    TRUE ~ "N")
  return(iupac_code)
}

#' Get records for each of the targeted and wgs datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param shiny_dt shiny data table to get the records from
#' @param mt swap table to write to
#' @param i location on swap table
#' @param randomize whether to randomize SNP order in order to generate mismatch
#'
#' @export
get_records <- function(shiny_dt, mt, i, randomize) {

  # Retrieve matching shiny records
  cfid <- mt[i]$cfdna_PID
  gid <- mt[i]$gdna_PID
  bsid <- mt[i]$wgbs_PID

  if (randomize) {
    # Keep current cfdna id
    cfid <- mt[i]$cfdna_PID

    # Select a random number that is not the same as current one
    set.seed(1359)
    ri <- sample(setdiff(1:nrow(mt), i), 1)
    gid <- mt[ri]$gdna_PID

    ri2 <- sample(setdiff(1:nrow(mt), c(i, ri)), 1)
    bsid <- mt[ri2]$wgbs_PID
  }

  ctar <- get_shiny_record(shiny_dt, cfid, "CFDNA_CLINICAL", "targeted")
  cwgs <- get_shiny_record(shiny_dt, cfid, "CFDNA_CLINICAL", "wgs")
  gtar <- get_shiny_record(shiny_dt, gid, "GDNA_CLINICAL", "targeted")
  gwgs <- get_shiny_record(shiny_dt, gid, "GDNA_CLINICAL", "wgs")
  wgbs <- get_shiny_record(shiny_dt, bsid, "METHYL_CLINICAL", "wgbs")

  return(list(ctar = ctar, cwgs = cwgs, gtar = gtar, gwgs = gwgs, wgbs = wgbs))
}

#' Get data tables for all conta datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param snps targeted snps data table
#'
#' @return list of conta genotype locations
#' @export
get_data_tables <- function(mt, i, snps) {

  # Read files if necessary (i.e. results are not yet calculated)
  ctar <- NULL
  cwgs <- NULL
  gtar <- NULL
  gwgs <- NULL
  wgbs <- NULL

  # Read cfDNA targeted if necessary (i.e. concordances are not calculated)
  if ((mt[i, conc_ctar_cwgs] == 0) | (mt[i, conc_ctar_gtar] == 0))
    ctar <- load_conta_file(mt[i, ctar_loc], snps)

  # Read cfDNA wgs if necessary (i.e. concordances are not calculated)
  if ((mt[i, conc_ctar_cwgs] == 0) | (mt[i, conc_cwgs_gwgs] == 0))
    cwgs <- load_conta_file(mt[i, cwgs_loc])

  # Read gDNA targeted if necessary (i.e. concordances are not calculated)
  if ((mt[i, conc_gtar_gwgs] == 0) | (mt[i, conc_ctar_gtar] == 0))
    gtar <- load_conta_file(mt[i, gtar_loc], snps)

  # Read gDNA wgs if necessary (i.e. concordances are not calculated)
  if ((mt[i, conc_gtar_gwgs] == 0) | (mt[i, conc_cwgs_gwgs] == 0))
    gwgs <- load_conta_file(mt[i, gwgs_loc])

  # Read wgbs if necessary (i.e. concordances are not calculated)
  if ((mt[i, conc_cwgs_wgbs] == 0) | (mt[i, conc_gwgs_wgbs] == 0))
    wgbs <- load_conta_file(mt[i, wgbs_loc])

  return(list(ctar = ctar,
              cwgs = cwgs,
              gtar = gtar,
              gwgs = gwgs,
              wgbs = wgbs))
}

#' Set signatures for all datasets
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param d list of data tables
#' @param signature_SNPs table containing rsid, ref and alt of signature SNPs
#'
#' @export
set_signatures <- function(mt, i, d, signature_SNPs) {

  # Calculate cfDNA targeted signature
  if (!is.null(d$ctar))
    mt[i, ctar_signature := get_signature(d$ctar, signature_SNPs)]

  # Calculate cfDNA wgs signature
  if (!is.null(d$cwgs))
    mt[i, cwgs_signature := get_signature(d$cwgs, signature_SNPs)]

  # Calculate gDNA targeted signature
  if (!is.null(d$gtar))
    mt[i, gtar_signature := get_signature(d$gtar, signature_SNPs)]

  # Calculate gDNA wgs signature
  if (!is.null(d$gwgs))
    mt[i, gwgs_signature := get_signature(d$gwgs, signature_SNPs)]

  # Calculate wgbs signature
  if (!is.null(d$wgbs))
    mt[i, wgbs_signature := get_signature(d$wgbs, signature_SNPs)]
}


#' Calculate concordance for all datasets
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param d list of data tables
#' @param cutoff cutoff for concordance
#' @param signature_SNPs table containing rsid, ref and alt of signature SNPs
#'
#' @export
calculate_concordance <- function(mt, i, d, cutoff, signature_SNPs) {

  # Calculate genotype signatures for each sample type
  set_signatures(mt, i, d, signature_SNPs)

  # Calculate concordances only if both files were read
  # cfDNA targeted vs cfDNA WGS concordance
  swap_concordance(mt, i, d$ctar, d$cwgs, "conc_ctar_cwgs", cutoff)

  # gDNA targeted vs gDNA WGS concordance
  swap_concordance(mt, i, d$gtar, d$gwgs, "conc_gtar_gwgs", cutoff)

  # cfDNA targeted vs gDNA targeted concordance
  swap_concordance(mt, i, d$ctar, d$gtar, "conc_ctar_gtar", cutoff)

  # cfDNA wgs vs gDNA wgs concordance
  swap_concordance(mt, i, d$cwgs, d$gwgs, "conc_cwgs_gwgs", cutoff)

  # cfDNA wgs vs wgbs concordance
  swap_concordance(mt, i, d$cwgs, d$wgbs, "conc_cwgs_wgbs", cutoff)

  # gDNA wgs vs wgbs concordance
  swap_concordance(mt, i, d$gwgs, d$wgbs, "conc_gwgs_wgbs", cutoff)
}


#' Set result locations and metrics from shiny_dt
#'
#' @param shiny_dt data.table with shiny results
#' @param mt swap table to write to
#' @param i location on swap table
#' @param randomize whether to randomize the data
#'
#' @export
set_locations_and_metrics <- function(shiny_dt, mt, i, randomize) {

  # Retrieve ids, and randomize if necessary
  r <- get_records(shiny_dt, mt, i, randomize)

  # Get s3 locations for results
  set_result_locations(mt, i, r)

  # Transfer metrics from shiny to swap table
  set_metrics(mt, i, r)
}

#' Set various metrics from shiny to swap table
#' @param mt data.table with swap data
#' @param i location on the table
#' @param r list of records containing ctar, gtar, cwgs and gwgs results
#' @param min_max_lh minimum max_lh to make a contamination call
#' @param min_lh_diff minimum lh_diff to make a contamination call
#'
#' @return mt with metrics added
#'
#' @export
set_metrics <- function(mt, i, r,
                        min_max_lh = 0.00025,
                        min_lh_diff = 50) {

  # Set cfdna targeted record IDs
  if (!is.na(r$ctar$SeqSampleID)) {
    mt[i, ctar_SeqSampleID := r$ctar$SeqSampleID]
    mt[i, ctar_mtc :=
         r$ctar$collapsed_fragment_mean_target_coverage]
    mt[i, ctar_conta_call :=
         (r$ctar$max_lh >= min_max_lh) &
         (r$ctar$lh_diff >= min_lh_diff)]
    mt[i, ctar_max_lh := r$ctar$max_lh]
    mt[i, ctar_lh_diff := r$ctar$lh_diff]
    mt[i, ctar_num_amplified_genes := r$ctar$num_amplified_genes]
  }

  # Set cfdna wgs record IDs
  if (!is.na(r$cwgs$SeqSampleID)) {
    mt[i, cwgs_SeqSampleID := r$cwgs$SeqSampleID]
    mt[i, cwgs_cov := r$cwgs$DEDUP_FILTERED_BAM_coverage]
    mt[i, cwgs_final_stein := r$cwgs$final_stein]
    mt[i, cwgs_final_mapd := r$cwgs$final_mapd]
  }

  # Set gdna targeted record IDs
  if (!is.na(r$gtar$SeqSampleID)) {
    mt[i, gtar_SeqSampleID := r$gtar$SeqSampleID]
    mt[i, gtar_mtc :=
         r$gtar$collapsed_fragment_mean_target_coverage]

    mt[i, gtar_conta_call :=
         (r$gtar$max_lh >= min_max_lh) &
         (r$gtar$lh_diff >= min_lh_diff)]
    mt[i, gtar_max_lh := r$gtar$max_lh]
    mt[i, gtar_lh_diff := r$gtar$lh_diff]
    mt[i, gtar_num_amplified_genes := r$gtar$num_amplified_genes]
  }

  # Set gdna wgs record IDs
  if (!is.na(r$gwgs$SeqSampleID)) {
    mt[i, gwgs_SeqSampleID := r$gwgs$SeqSampleID]
    mt[i, gwgs_cov := r$gwgs$DEDUP_FILTERED_BAM_coverage]
    mt[i, gwgs_final_stein := r$gwgs$final_stein]
    mt[i, gwgs_final_mapd := r$gwgs$final_mapd]
  }

  # Set wgbs record IDs
  if (!is.na(r$wgbs$SeqSampleID)) {
    mt[i, wgbs_SeqSampleID := r$wgbs$SeqSampleID]
    mt[i, wgbs_cov := r$wgbs$DEDUP_FILTERED_BAM_coverage]
    mt[i, wgbs_final_stein := r$wgbs$final_stein]
    mt[i, wgbs_final_mapd := r$wgbs$final_mapd]
  }
}


#' Get records for each of the targeted and wgs datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param mt data.table with swap data
#' @param i location on the table
#' @param r result_locations for each sample type
#'
#' @export
set_result_locations <- function(mt, i, r) {
  mt[i, ctar_loc := get_conta_file_loc(r$ctar, targeted = TRUE)]
  mt[i, cwgs_loc := get_conta_file_loc(r$cwgs, targeted = FALSE)]
  mt[i, gtar_loc := get_conta_file_loc(r$gtar, targeted = TRUE)]
  mt[i, gwgs_loc := get_conta_file_loc(r$gwgs, targeted = FALSE)]
  mt[i, wgbs_loc := get_conta_file_loc(r$wgbs, targeted = FALSE)]
}

#' Write results to file if any updates were made
#'
#' @param d list of data tables
#' @param mt swap table to write to
#' @param out_table_name location to write (either s3 or local)
#'
#' @export
write_results <- function(d, mt, out_table_name) {

  # If any new calculations were made, then the data tables read should be
  # non-null, in that case, write to results
  if (!is.null(d$ctar) | !is.null(d$cwgs)
      | !is.null(d$gtar) | !is.null(d$gwgs)) {
    if (startsWith(out_table_name, "s3")) {
      write_to_s3(mt, out_table_name, write.table, sep = "\t",
                  col.names = TRUE, row.names = FALSE, quote = FALSE)
    } else {
      write.table(mt, out_table_name, sep = "\t",
                  col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  }
}
