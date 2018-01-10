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

using namespace Rcpp;

// Convert given chromosome to integer representation.
int numericChr(std::string chr) {
  if (chr.find("chr") == 0 || chr.find("Chr") == 0 || chr.find("CHR") == 0)
    chr = chr.substr(3);

  if (chr[0] =='X') return 23;
  if (chr[0] =='Y') return 24;
  if (chr[0] =='M') return 25;

  return atoi(chr.c_str());
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
    fin.open(file.c_str());
    if (!fin) {
      std::cerr << "Could not read file" << file <<  "\n";
      exit(EXIT_FAILURE);
    }

    // Skip VCF headers.
    getline(fin, line);
    while (line.find("#") == 0)
      getline(fin, line);
  }

  std::shared_ptr<VcfRecord> getNext() {
    if (line.empty() | fin.eof())
      return nullptr;

    std::shared_ptr<VcfRecord> vcf_record (new VcfRecord(line));
    total++;
    getline(fin, line);

    // Skip if the new record is multi allelic.
    while (vcf_record->indel_or_multi_allelic) {
      skipped++;
      if (line.empty() | fin.eof())
        return nullptr;
      vcf_record.reset(new VcfRecord(line));
      total++;
      getline(fin, line);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prev_record is null for the first record only.
    if (prev_record != nullptr &&
        compare(vcf_record->num_chr, vcf_record->pos,
                prev_record->num_chr, prev_record->pos) < 0 ) {
      std::cerr << "Records out of order: " << prev_record->toString() << "\n"
                << "larger than: " << vcf_record->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prev_record = vcf_record;
    return vcf_record;
  }

  int total = 0, skipped = 0;

 private:
  std::shared_ptr<VcfRecord> prev_record;
  std::string line;
  std::ifstream fin;
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
    num_chr = numericChr(chr);

    // Set whether this is a multi-allelic SNP.
    multi_allelic = alts.length() > 3;
  }

  std::string toString() {
    return (line);
  }

  std::string toTsvFormat(std::shared_ptr<VcfRecord> vcf_record) {
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
    fin.open(file.c_str());
    if (!fin) {
      std::cerr << "Could not read file" << file <<  "\n";
      exit(EXIT_FAILURE);
    }

    pileup_header = "CHROM\tPOS\tDEPTH\tREF\tA\tC\tG\tT\tN\tINS\tDEL";
    tsv_header = "chrom\tpos\tA\tT\tG\tC\tN\tref\tmajor\talleles\trsid";
    tsv_header_with_maf = tsv_header + "\tmaf";
    total = 0;
    skipped = 0;

    // Skip header of TSV file.
    getline(fin, header);
    if (header == pileup_header) {
      PILEUP = TRUE;
    } else if (header == tsv_header | header == tsv_header_with_maf) {
      PILEUP = FALSE;
    } else {
      std::cerr << "\nHeader does not match expected: " << header <<  "\n";
      exit(EXIT_FAILURE);
    }

    // Move to next line.
    getline(fin, line);
  }

  std::shared_ptr<TsvRecord> getNext() {
    if (line.empty() | fin.eof())
      return nullptr;

    std::shared_ptr<TsvRecord> tsv_record(new TsvRecord(line, PILEUP));
    total++;
    getline(fin, line);

    // Skip if the new record is multi allelic
    while (tsv_record->multi_allelic) {
      skipped++;
      if (line.empty() | fin.eof())
        return nullptr;
      tsv_record.reset(new TsvRecord(line, PILEUP));
      total++;
      getline(fin, line);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prev_record is null for the first record only.
    if (prev_record != nullptr &&
        compare(tsv_record->num_chr, tsv_record->pos,
                prev_record->num_chr, prev_record->pos) < 0 ) {
      std::cerr << "Records out of order: " << prev_record->toString() << "\n"
                << "larger than: " << tsv_record->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prev_record = tsv_record;
    return tsv_record;
  }

  bool PILEUP = FALSE;
  int total = 0, skipped = 0;
  std::string pileup_header, tsv_header, tsv_header_with_maf, header;

 private:
  std::string line;
  std::ifstream fin;
  std::shared_ptr<TsvRecord> prev_record;
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
//' @param DEBUG boolean whether to print out messages
//' @return null
//'
//' @export
// [[Rcpp::export]]
void intersect_snps(const char* tsv_filename, const char* out_tsv_filename,
               const char* vcf_filename, bool non_dbSNP = false,
               bool DEBUG = false) {
    clock_t start_time, end_time;
    start_time = clock();

    if (DEBUG)
      std::cout << "\nReading TSV file from " << tsv_filename << "\n";
    TsvReader tsv_reader(tsv_filename);

    // Set up out file to write
    std::ofstream outfile;
    outfile.open(out_tsv_filename);
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
          outfile << tsv_record->toTsvFormat(nullptr) << "\n";
          tsv_record = tsv_reader.getNext();
          written++;
        }
        break;
      }

      int cmp = compare(tsv_record->num_chr, tsv_record->pos,
                        vcf_record->num_chr, vcf_record->pos);

      if (cmp < 0) {  // VCF is ahead
        if (non_dbSNP)
          outfile << tsv_record->toTsvFormat(nullptr) << "\n";
        tsv_record = tsv_reader.getNext();
      } else if (cmp > 0) {  // TSV is ahead
        vcf_record = vcf_reader.getNext();
      } else {  // a matching record
        outfile << tsv_record->toTsvFormat(vcf_record) << "\n";
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
