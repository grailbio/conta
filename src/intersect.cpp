// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <memory>
#include "genome.h"

using namespace Rcpp;

// Convert given chromosome to integer representation.
// Conta does not expect SNPs from additional chromosomes in its input,
// so they are not guaranteed to be sorted correctly in this output.
// TODO: make this code sorted correctly for the extra decoy chromosomes
int numericChr(std::string chr) {

  if (chr[0] =='X') return 23;
  if (chr[0] =='Y') return 24;
  if (chr[0] =='M') return 25;

  int chr_num = atoi(chr.c_str());
  if (chr_num <= 0 || chr_num > 25)
    return (abs(chr_num) + 26);
  else
    return chr_num;
}

// Compare two genomic positions.
int compare(int this_chr, int this_pos, int other_chr, int other_pos) {
  if (this_chr < other_chr) return -1;
  if (this_chr > other_chr) return 1;

  // chrs are the same.
  return (this_pos - other_pos);
}

// A generic VCF record for a VCF without genotype field.
class VcfRecord {
 public:
  explicit VcfRecord(std::string line) {
    this->line = line;
    std::istringstream ss(line);
    ss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info;
    if (chr.find("chr") == 0 || chr.find("Chr") == 0 || chr.find("CHR") == 0)
      chr = chr.substr(3);
    num_chr = numericChr(chr);

    int cafStart = info.find("CAF=");

    if (cafStart != -1) {
      int cafEnd = info.find(',', cafStart);
      // Since we skipped multi-allelic it is guaranteed there is only one MAF.
      if (cafEnd != -1) {
        CAF = atof(info.substr(cafStart+4, cafEnd-cafStart-4).c_str());
      } else {
        CAF = atof(info.substr(cafStart+4).c_str());
      }
      MAF = 1 - CAF;
    } else {
      int afStart = info.find("AF=");

      if (afStart == -1) {
        int afEnd = info.find(',', afStart);
        // Since we skipped multi-allelic it is guaranteed there is only one MAF.
        if (afEnd != -1) {
          MAF = atof(info.substr(afStart, afEnd).c_str());
        } else {
          MAF = atof(info.substr(afStart).c_str());
        }
      }
    }

    // Set whether this is an indel or multi-allelic SNP.
    indel_or_multi_allelic = alt.length() > 1;
  }

  std::string toString() {
    return (line);
  }

  bool indel_or_multi_allelic;
  int pos, num_chr;
  double CAF = 0, MAF = 0;
  std::string line, chr, rsid, qual, filter, info, ref, alt;
};

// A container that streams VCF records.
class VcfReader {
 public:
  explicit VcfReader(std::string file) {
    fin_.open(file.c_str());
    if (!fin_) {
      std::cerr << "Could not read file" << file <<  "\n";
      exit(EXIT_FAILURE);
    }

    // Skip VCF headers.
    getline(fin_, line_);
    while (line_.find("#") == 0)
      getline(fin_, line_);
  }

  std::shared_ptr<VcfRecord> getNext() {
    if (line_.empty() | fin_.eof())
      return nullptr;

    std::shared_ptr<VcfRecord> vcf_record (new VcfRecord(line_));
    total++;
    getline(fin_, line_);

    // Skip if the new record is multi allelic.
    while (vcf_record->indel_or_multi_allelic) {
      skipped++;
      if (line_.empty() | fin_.eof())
        return nullptr;
      vcf_record.reset(new VcfRecord(line_));
      total++;
      getline(fin_, line_);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prev_record_ is null for the first record only.
    if (prev_record_ != nullptr &&
        compare(vcf_record->num_chr, vcf_record->pos,
                prev_record_->num_chr, prev_record_->pos) < 0 ) {
      std::cerr << "Records out of order: " << prev_record_->toString() << "\n"
                << "larger than: " << vcf_record->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prev_record_ = vcf_record;
    return vcf_record;
  }

  int total = 0, skipped = 0;

 private:
  std::shared_ptr<VcfRecord> prev_record_;
  std::string line_;
  std::ifstream fin_;
};

// A TSV record holds allele frequency oberved for a sample (custom format)
class TsvRecord {
 public:
  explicit TsvRecord(std::string line, bool PILEUP = FALSE) {
    this->line = line;
    std::istringstream ss(line);
    if (PILEUP) {
      ss >> chr >> pos >> depth >> ref >> A >> C >> G >> T >> N >> INS >> DEL;
    } else {
      ss >> chr >> pos >> A >> T >> G >> C >> N >> ref >> major >> alts >> rsid;
    }
    if (chr.find("chr") == 0 || chr.find("Chr") == 0 || chr.find("CHR") == 0)
      chr = chr.substr(3);
    num_chr = numericChr(chr);

    // Set whether this is a multi-allelic SNP.
    multi_allelic = alts.length() > 3;
  }

  std::string toString() {
    return (line);
  }

  std::string toTsvFormat(std::shared_ptr<VcfRecord> vcf_record,
                          std::shared_ptr<Genome> genome) {
    std::stringstream pp;
    if (vcf_record != nullptr) {
      std::string alts = vcf_record->ref + "/" + vcf_record->alt;
      pp << chr << "\t" << pos << "\t" << A
        << "\t" << T << "\t" << G << "\t" << C << "\t" << N << "\t"
        << vcf_record->ref << "\t" << vcf_record->ref << "\t"
        << alts << "\t" << vcf_record->rsid << "\t" << vcf_record->MAF;
    } else {
      std::string alts = ref + "/N";
      pp << chr << "\t" << pos << "\t" << A
         << "\t" << T << "\t" << G << "\t" << C << "\t" << N << "\t"
         << ref << "\t" << ref << "\t"
         << alts << "\t" << "NA" << "\t" << "0";
    }
    if (genome != nullptr) {
      // genome is 0-based, then start from previous position, get 3 bases
      pp << "\t" << genome->getPos(chr, pos-2, 3);
    }
    return pp.str();
  }

  std::string toPileupFormat() {
    // Re-calculate depth if 0 (happens if file was tsv format).
    if (depth == 0) {
      int depth = A + T + G + C + N;
      int INS, DEL = 0;
    }
    std::stringstream pp;
    pp << chr << "\t" << pos << "\t" << depth << "\t" << ref << "\t"
       << A << "\t" << C << "\t" << G << "\t" << T << "\t"
       << N << "\t" << INS << "\t" << DEL;
    return pp.str();
  }

  std::string line;
  std::string chr, alts, rsid;
  int depth = 0;
  int pos, A, T, G, C, N, num_chr, INS, DEL;
  std::string ref, major;
  bool multi_allelic;
};

// Object that holds a vector of TSV records.
class TsvReader {
 public:
  explicit TsvReader(std::string file) {
    fin_.open(file.c_str());
    if (!fin_) {
      std::cerr << "Could not read file" << file <<  "\n";
      exit(EXIT_FAILURE);
    }

    pileup_header = "CHROM\tPOS\tDEPTH\tREF\tA\tC\tG\tT\tN\tINS\tDEL";
    tsv_header = "chrom\tpos\tA\tT\tG\tC\tN\tref\tmajor\talleles\trsid";
    tsv_header_with_maf = tsv_header + "\tmaf";
    tsv_header_with_genome = tsv_header_with_maf + "\tcontext";
    total = 0;
    skipped = 0;

    // Skip header of TSV file.
    getline(fin_, header);
    if (header == pileup_header) {
      PILEUP = TRUE;
    } else if (header == tsv_header | header == tsv_header_with_maf |
               header == tsv_header_with_genome) {
      PILEUP = FALSE;
    } else {
      std::cerr << "\nHeader does not match expected: " << header <<  "\n";
      exit(EXIT_FAILURE);
    }

    // Move to next line.
    getline(fin_, line_);
  }

  std::shared_ptr<TsvRecord> getNext() {
    if (line_.empty() | fin_.eof())
      return nullptr;

    std::shared_ptr<TsvRecord> tsv_record(new TsvRecord(line_, PILEUP));
    total++;
    getline(fin_, line_);

    // Skip if the new record is multi allelic
    while (tsv_record->multi_allelic) {
      skipped++;
      if (line_.empty() | fin_.eof())
        return nullptr;
      tsv_record.reset(new TsvRecord(line_, PILEUP));
      total++;
      getline(fin_, line_);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prev_record_ is null for the first record only.
    if (prev_record_ != nullptr &&
        compare(tsv_record->num_chr, tsv_record->pos,
                prev_record_->num_chr, prev_record_->pos) < 0 ) {
      std::cerr << "Records out of order: " << prev_record_->toString() << "\n"
                << "larger than: " << tsv_record->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prev_record_ = tsv_record;
    return tsv_record;
  }

  bool PILEUP = FALSE;
  int total = 0, skipped = 0;
  std::string pileup_header, tsv_header, tsv_header_with_maf,
    tsv_header_with_genome, header;

 private:
  std::string line_;
  std::ifstream fin_;
  std::shared_ptr<TsvRecord> prev_record_;
};

//' Intersect a TSV file with VCF file and write a TSV output
//' that contains maf values for each SNP in the TSV file.
//' In the process, remove multiple-allele SNPs if the option
//' is turned on.
//'
//' @param tsv_filename string name of input TSV file
//' @param out_tsv_filename string name of output TSV file
//' @param vcf_filename string dbsnp VCF file name to be intersected
//' @param non_dbSNP allow non dbSNP positions to be included
//' @param genome_filename string name of the genome.fa
//' @param DEBUG boolean whether to print out messages
//' @return null
//'
//' @export
// [[Rcpp::export]]
void intersect_snps(const char* tsv_filename, const char* out_tsv_filename,
                    const char* vcf_filename, bool non_dbSNP = false,
                    const char* genome_filename = "NA", bool DEBUG = false) {
    clock_t start_time, end_time;
    start_time = clock();

    if (DEBUG)
      std::cout << "\nReading genome file from " << genome_filename << "\n";

    std::shared_ptr<Genome> genome = nullptr;
    if (std::strcmp(genome_filename, "NA") != 0)
      genome = (std::shared_ptr<Genome>) (new Genome(genome_filename, DEBUG));

    if (DEBUG)
      std::cout << "genome is null: " << (genome == nullptr) << "\n";

    if (DEBUG)
      std::cout << "\nReading TSV file from " << tsv_filename << "\n";
    TsvReader tsv_reader(tsv_filename);

    // Set up out file to write
    std::ofstream outfile;
    outfile.open(out_tsv_filename);
    if (genome != nullptr)
      outfile << tsv_reader.tsv_header_with_genome << "\n";
    else
      outfile << tsv_reader.tsv_header_with_maf << "\n";

    if (DEBUG)
      std::cout << "Reading VCF file from " << vcf_filename << "\n";
    VcfReader vcf_reader(vcf_filename);

    // Get initial records from the files.
    int written = 0;
    std::shared_ptr<VcfRecord> vcf_record = vcf_reader.getNext();
    std::shared_ptr<TsvRecord> tsv_record = tsv_reader.getNext();

    if (DEBUG)
      std::cout << "Running the linear join\n";

    while (1) {
      if (tsv_record == nullptr || vcf_record == nullptr) {
        while (non_dbSNP && tsv_record != nullptr) {
          outfile << tsv_record->toTsvFormat(nullptr, genome) << "\n";
          tsv_record = tsv_reader.getNext();
          written++;
        }
        break;
      }

      if (DEBUG) {
        std::cout << "Comparing:" << "\n";
        std::cout << "tsv_record: " << tsv_record->toString() << "\n";
        std::cout << "vcf_record: " << vcf_record->toString() << "\n";
      }

      int cmp = compare(tsv_record->num_chr, tsv_record->pos,
                        vcf_record->num_chr, vcf_record->pos);

      if (cmp < 0) {  // VCF is ahead
        if (non_dbSNP)
          outfile << tsv_record->toTsvFormat(nullptr, genome) << "\n";
        tsv_record = tsv_reader.getNext();
      } else if (cmp > 0) {  // TSV is ahead
        vcf_record = vcf_reader.getNext();
      } else {  // a matching record
        outfile << tsv_record->toTsvFormat(vcf_record, genome) << "\n";
        tsv_record = tsv_reader.getNext();
        vcf_record = vcf_reader.getNext();
        written++;
      }
    }

    if (DEBUG) {
      std::cout << "Read " << vcf_reader.total << " VCF records, and "
                << tsv_reader.total << " TSV records and wrote " << written
                << " TSV records.\nTotal of " << tsv_reader.skipped
                << " TSV records were skipped (multi-allelic).\n"
                << "Total of " << vcf_reader.skipped
                << " VCF records were skipped (indel or multi-allelic)\n";
    }
    outfile.close();

    end_time = clock();
    float elapsed_time = static_cast<float> (end_time - start_time);
    float elapsed_secs = elapsed_time / CLOCKS_PER_SEC;

    if (DEBUG)
      std::cout << "Elapsed seconds: " << elapsed_secs << "\n";
}

//' Return sequence specified by the chromosome start and lengths
//'
//' @param genome_filename string name of the genome.fa
//' @param chr chromosome name
//' @param start sequence start base
//' @param length length of sequence
//' @param DEBUG print out debug messages
//' @return sequence string
//'
//' @export
// [[Rcpp::export]]
std::string get_genomic_seq(const char* genome_filename = "NA",
                            const char* chr = "NA",
                            const int start = 0,
                            const int length = 0,
                            bool DEBUG = false) {

  Genome genome(genome_filename, false);

  return(genome.getPos(chr, start, length));
}

//' Return sequences specified by the chromosome start and lengths
//'
//' @param genome_filename string name of the genome.fa
//' @param chrs chromosome name
//' @param starts sequence start base
//' @param lengths length of sequence
//' @param DEBUG print out debug messages
//' @return sequence string
//'
//' @export
// [[Rcpp::export]]
CharacterVector get_genomic_seqs(const char* genome_filename,
                                 CharacterVector chrs,
                                 NumericVector starts,
                                 NumericVector lengths,
                                 bool DEBUG = false) {

  Genome genome(genome_filename, false);

  int n = chrs.size();
  CharacterVector seqs(n);

  for (int i=0; i<n; i++)
    seqs[i] = genome.getPos((const char *) chrs[i], starts[i], lengths[i]);

  return(seqs);
}
