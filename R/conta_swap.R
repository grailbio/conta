#'
#' Run swap tool
#'
#' Compare genotypes across specified sample pairings (specific to GRAIL CCGA
#' data). Both concordance and distance metrics are reported along with other
#' informative features such as coverage, CNV indicators and contamination.
#'
#' @param pairing_table_name name of the mapping table
#' @param out_table_name name of the results table to output
#' @param shiny_loc location of the shiny database to retrieve sequencing metrics
#' @param dbsnp_targeted targeted dbsnp vcf file same as used for targeted seq
#' @param randomize generate random pairs instead of expected ones
#' @param ctar_cwgs_cutoff cutoff to make a cfdna targeted vs wgs swap call
#' @param gtar_gwgs_cutoff cutoff to make a gdna targeted vs wgs swap call
#' @param ctar_gtar_cutoff cutoff to make a cfdna targeted vs gdna swap call
#' @param cwgs_gwgs_cutoff cutoff to make a wgs targeted vs gdna wgs swap call
#' @param max_conta if contamination is above this level, do not call swap
#'
#' @export
swap <- function(pairing_table_name, out_table_name, shiny_loc, dbsnp_targeted,
                 randomize = FALSE,
                 ctar_cwgs_cutoff = 0.7, gtar_gwgs_cutoff = 0.7,
                 ctar_gtar_cutoff = 0.7, cwgs_gwgs_cutoff = 0.7,
                 max_conta = 0.2) {

  if (randomize) {
    out_table_name <- paste(gsub(".csv", "", out_table_name),
                           ".random.csv", sep = "")
  }

  # Read mapping table. This table is an input that stores the pairings between
  # cfDNA and gDNA samples, along with a fake pair ID (DID) and StudyBatch
  pairing_table <- read_data_table(pairing_table_name, sep = ",", showProgress = FALSE)

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
    calculate_concordance(mt, i, d, ctar_cwgs_cutoff, gtar_gwgs_cutoff,
                          ctar_gtar_cutoff, cwgs_gwgs_cutoff)

    # Do not call swap if contamination is above a threshold
    if (mt[i, swap] & ((mt[i, ctar_MaxLh] >= max_conta) |
                       (mt[i, gtar_MaxLh] >= max_conta)))
      mt[i, swap := NA]

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

  mt <- read_data_table(out_table_name, sep = ",",
                                 showProgress = FALSE)

  # If first time creating the output, copy it from match original table,
  # and add new fields.
  if (is.null(mt)) {
    nl <- nrow(pairing_table)
    mt <- data.table(DID = pairing_table$DID,
                     StudyBatch = pairing_table$StudyBatch,
                     gdna_PID = pairing_table$gdna_PID,
                     cfdna_PID = pairing_table$cfdna_PID,
                     swap = logical(length = nl),
                     ctar_SeqSampleID = character(length = nl),
                     cwgs_SeqSampleID = character(length = nl),
                     gtar_SeqSampleID = character(length = nl),
                     gwgs_SeqSampleID = character(length = nl),
                     ctar_MTC = numeric(length = nl),
                     cwgs_COV = numeric(length = nl),
                     gtar_MTC = numeric(length = nl),
                     gwgs_COV = numeric(length = nl),
                     ctar_num_amplified_genes = double(length = nl),
                     cwgs_final_stein = double(length = nl),
                     cwgs_final_mapd = double(length = nl),
                     gtar_num_amplified_genes = double(length = nl),
                     gwgs_final_stein = double(length = nl),
                     gwgs_final_mapd = double(length = nl),
                     conc_ctar_cwgs = double(length = nl),
                     conc_gtar_gwgs = double(length = nl),
                     conc_ctar_gtar = double(length = nl),
                     conc_cwgs_gwgs = double(length = nl),
                     dist_ctar_cwgs = double(length = nl),
                     dist_gtar_gwgs = double(length = nl),
                     dist_ctar_gtar = double(length = nl),
                     dist_cwgs_gwgs = double(length = nl),
                     ctar_ContaminationCall = logical(length = nl),
                     ctar_MaxLh = double(length = nl),
                     ctar_LhDiff = numeric(length = nl),
                     gtar_ContaminationCall = logical(length = nl),
                     gtar_MaxLh = double(length = nl),
                     gtar_LhDiff = numeric(length = nl),
                     ctar_signature = character(length = nl),
                     gtar_signature = character(length = nl),
                     cwgs_signature = character(length = nl),
                     gwgs_signature = character(length = nl),
                     ctar_loc = character(length = nl),
                     gtar_loc = character(length = nl),
                     cwgs_loc = character(length = nl),
                     gwgs_loc = character(length = nl))
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

#' Logicize columns from shiny data if they are set as string
#'
#' @param shiny_dt data.table
#' @param cols columns to logicize
#'
#' @return shiny_dt that is numerized
#'
#' @export
logicize_columns <- function(shiny_dt, cols) {

  # Columns to set as logical
  locs <- match(cols, colnames(shiny_dt))
  for (loc in seq_along(locs)) {
    set(shiny_dt, NULL, locs[loc], as.logical(shiny_dt[[locs[loc]]]))
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

  if ((!is.null(snps)) & !is.na(conta_loc) & s3_file_exists(conta_loc)) {

    # Read targeted genotypes and rsid from vfn file format
    conta_loc_dt <- read_data_table(conta_loc, sep = ",", showProgress = FALSE)

    # Add rsid. Targeted conta output doesn't have this info
    if (!is.null(conta_loc_dt))
      conta_loc_dt[, rsid := snps$ID]

  } else if (!is.na(conta_loc) & s3_file_exists(conta_loc)) {

    # Read conta genotypes and rsid from gt file format
    conta_loc_dt <- read_data_table(conta_loc, showProgress = FALSE)

    # Vf is the variant allele frequency, needs to be set of wgs files
    conta_loc_dt[, vf := vr / dp]

  } else {

    conta_loc_dt <- NULL

  }

  return(conta_loc_dt)
}

#' Calculate concordance between two samples' genotypes
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param dt1 genotype set 1
#' @param dt2 genotype set 2
#' @param col1 column name 1
#' @param col2 column name 2
#' @param conc_cutoff concordance cutoff to make a swap call
#'
#' @export
swap_concordance <- function(mt, i, dt1, dt2, col1, col2, conc_cutoff,
                             min_dp = 10, min_maf = 0.25) {

  # Check if any of the datasets are null, return original table
  if (is.null(dt1) | is.null(dt2))
    return()

  # Do not calculate concordance if depth is less than 10 for one of the samples
  if (mean(dt1$dp, na.rm = TRUE) < min_dp | mean(dt2$dp, na.rm = TRUE) < min_dp)
    return()

  # If maf is specified (only targeted conta results have it), filter them
  maf_field1 <- which(colnames(dt1) %in% "maf")
  if (length(maf_field1) > 0)
    dt1 <- dt1[maf >= min_maf, ]

  # Merge the two data tables to put each SNP in a single row
  m1 <- merge(dt1, dt2, by = "rsid")

  # Calculate genotype concordance which is the fraction of matching genotypes
  conc <- m1[, mean(gt.x == gt.y, na.rm = TRUE)]

  # Calculate euclidian distance between variant allele frequencies
  dist <- c(dist(rbind(m1$vf.x, m1$vf.y), method = "euclidian"))

  # Write concordance to results table
  mt[i, eval(col1) := conc]

  # Write distance to results table
  mt[i, eval(col2) := dist]

  # Call swap if the concordance is less than cut-off
  if (!is.na(conc) & (conc <= conc_cutoff))
    mt[i, swap := TRUE]
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
  cols <- c('collapsed_fragment_MEAN_TARGET_COVERAGE', 'num_amplified_genes',
            'DEDUP_FILTERED_BAM_coverage', 'final_stein', 'MaxLh', 'LhDiff',
            'final_mapd')
  numerize_columns(shiny_dt, cols)

  # Set following columns as logical
  cols2 <- c('ContaminationCall')
  logicize_columns(shiny_dt, cols2)

  return(shiny_dt)
}

#' Get signature for a given set of SNPs from a conta vfn or gt file
#'
#' @param conta_dt targeted vfn data table
#' @return data.table with conta genotypes
#'
#' @export
get_signature <- function(conta_dt) {

  # Select 10 SNPs that have minor allele frequency around 50%
  signature_SNPs <- data.table( rsid = c("rs2237290", "rs6474353", "rs1335873",
                                         "rs9551465", "rs354439", "rs1454361",
                                         "rs2494737", "rs2836358", "rs987640"),
                                ref = c("T", "T", "T", "A", "A",
                                        "T", "T", "A", "T"),
                                alt = c("A", "A", "A", "T", "T",
                                        "A", "A", "T", "A"))

  m1 <- merge(conta_dt, signature_SNPs, by = "rsid", all.y = TRUE)
  return(paste(ifelse(is.na(m1$gt), 'N',
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
    (base1 == 'A' & base2 == 'T') ~ 'W',
    (base1 == 'T' & base2 == 'A') ~ 'W',
    TRUE ~ 'N')
  return(iupac_code)
}

#' Get records for each of the targeted and wgs datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param shiny_dt shiny data table to get the records from
#' @param mt swap table to write to
#' @param i location on swap table
#'
#' @export
get_records <- function(shiny_dt, mt, i, randomize) {

  # Retrieve matching shiny records
  cfid <- mt[i]$cfdna_PID
  gid <- mt[i]$gdna_PID

  if (randomize) {
    # Select a random number that is not the same as current one
    set.seed(1359)
    ri <- sample(setdiff(1:nrow(mt), i), 1)
    cfid <- mt[i]$cfdna_PID
    gid <- mt[ri]$gdna_PID
  }

  ctar <- get_shiny_record(shiny_dt, cfid, "CFDNA_CLINICAL", "targeted")
  cwgs <- get_shiny_record(shiny_dt, cfid, "CFDNA_CLINICAL", "wgs")
  gtar <- get_shiny_record(shiny_dt, gid, "GDNA_CLINICAL", "targeted")
  gwgs <- get_shiny_record(shiny_dt, gid, "GDNA_CLINICAL", "wgs")

  return(list(ctar = ctar, cwgs = cwgs, gtar = gtar, gwgs = gwgs))
}

#' Get data tables for all conta datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param co previously calculated concordances
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

  return(list(ctar = ctar,
              cwgs = cwgs,
              gtar = gtar,
              gwgs = gwgs))
}

#' Set signatures for all datasets
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param d list of data tables
#'
#' @export
set_signatures <- function(mt, i, d) {

  # Calculate cfDNA targeted signature
  if (!is.null(d$ctar))
    mt[i, ctar_signature := get_signature(d$ctar)]

  # Calculate cfDNA wgs signature
  if (!is.null(d$cwgs))
    mt[i, cwgs_signature := get_signature(d$cwgs)]

  # Calculate gDNA targeted signature
  if (!is.null(d$gtar))
    mt[i, gtar_signature := get_signature(d$gtar)]

  # Calculate gDNA wgs signature
  if (!is.null(d$gwgs))
    mt[i, gwgs_signature := get_signature(d$gwgs)]
}


#' Calculate concordance for all datasets
#'
#' @param mt swap table to write to
#' @param i location on swap table
#' @param d list of data tables
#' @param ctar_cwgs_cutoff cutoff for ctar vs cwgs concordance
#' @param gtar_gwgs_cutoff cutoff for gtar vs gwgs concordance
#' @param ctar_gtar_cutoff cutoff for ctar vs gtar concordance
#' @param cwgs_gwgs_cutoff cutoff for cwgs vs gwgs concordance
#'
#' @export
calculate_concordance <- function(mt, i, d, ctar_cwgs_cutoff, gtar_gwgs_cutoff,
                                  ctar_gtar_cutoff, cwgs_gwgs_cutoff) {

  # Calculate genotype signatures for each sample type
  set_signatures(mt, i, d)

  # Calculate concordances only if both files were read
  # cfDNA targeted vs cfDNA WGS concordance
  swap_concordance(mt, i, d$ctar, d$cwgs, "conc_ctar_cwgs", "dist_ctar_cwgs",
                   ctar_cwgs_cutoff)

  # gDNA targeted vs gDNA WGS concordance
  swap_concordance(mt, i, d$gtar, d$gwgs, "conc_gtar_gwgs", "dist_gtar_gwgs",
                   gtar_gwgs_cutoff)

  # cfDNA targeted vs gDNA targeted concordance
  swap_concordance(mt, i, d$ctar, d$gtar, "conc_ctar_gtar", "dist_ctar_gtar",
                   ctar_gtar_cutoff)

  # cfDNA wgs vs gDNA wgs concordance
  swap_concordance(mt, i, d$cwgs, d$gwgs, "conc_cwgs_gwgs", "dist_cwgs_gwgs",
                   cwgs_gwgs_cutoff)
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
#' @param min_MaxLh minimum MaxLh to make a contamination call
#' @param min_LhDiff minimum LhDiff to make a contamination call
#'
#' @return mt with metrics added
#'
#' @export
set_metrics <- function(mt, i, r,
                        min_MaxLh = 0.00025,
                        min_LhDiff = 50) {

  # Set cfdna targeted record IDs
  if (!is.na(r$ctar$SeqSampleID)) {
    mt[i, ctar_SeqSampleID := r$ctar$SeqSampleID]
    mt[i, ctar_MTC :=
         r$ctar$collapsed_fragment_MEAN_TARGET_COVERAGE]
    mt[i, ctar_ContaminationCall :=
         (r$ctar$MaxLh >= min_MaxLh) &
         (r$ctar$LhDiff >= min_LhDiff)]
    mt[i, ctar_MaxLh := r$ctar$MaxLh]
    mt[i, ctar_LhDiff := r$ctar$LhDiff]
    mt[i, ctar_num_amplified_genes := r$ctar$num_amplified_genes]
  }

  # Set cfdna wgs record IDs
  if (!is.na(r$cwgs$SeqSampleID)) {
    mt[i, cwgs_SeqSampleID := r$cwgs$SeqSampleID]
    mt[i, cwgs_COV := r$cwgs$DEDUP_FILTERED_BAM_coverage]
    mt[i, cwgs_final_stein := r$cwgs$final_stein]
    mt[i, cwgs_final_mapd := r$cwgs$final_mapd]
  }

  # Set gdna targeted record IDs
  if (!is.na(r$gtar$SeqSampleID)) {
    mt[i, gtar_SeqSampleID := r$gtar$SeqSampleID]
    mt[i, gtar_MTC :=
         r$gtar$collapsed_fragment_MEAN_TARGET_COVERAGE]

    mt[i, gtar_ContaminationCall :=
         (r$gtar$MaxLh >= min_MaxLh) &
         (r$gtar$LhDiff >= min_LhDiff)]
    mt[i, gtar_MaxLh := r$gtar$MaxLh]
    mt[i, gtar_LhDiff := r$gtar$LhDiff]
    mt[i, gtar_num_amplified_genes := r$gtar$num_amplified_genes]
  }

  # Set gdna wgs record IDs
  if (!is.na(r$gwgs$SeqSampleID)) {
    mt[i, gwgs_SeqSampleID := r$gwgs$SeqSampleID]
    mt[i, gwgs_COV := r$gwgs$DEDUP_FILTERED_BAM_coverage]
    mt[i, gwgs_final_stein := r$gwgs$final_stein]
    mt[i, gwgs_final_mapd := r$gwgs$final_mapd]
  }
}


#' Get records for each of the targeted and wgs datasets
#'
#' If randomize is selected, randomly select the entry
#'
#' @param r result_locations for each sample type
#'
#' @export
set_result_locations <- function(mt, i, r) {
  mt[i, ctar_loc := get_conta_file_loc(r$ctar, targeted = TRUE)]
  mt[i, cwgs_loc := get_conta_file_loc(r$cwgs, targeted = FALSE)]
  mt[i, gtar_loc := get_conta_file_loc(r$gtar, targeted = TRUE)]
  mt[i, gwgs_loc := get_conta_file_loc(r$gwgs, targeted = FALSE)]
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
      write_to_s3(mt, out_table_name, write.table, sep = ",",
                  col.names = TRUE, row.names = FALSE, quote = FALSE)
    } else {
      write.table(mt, out_table_name, sep = ",",
                  col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  }
}
