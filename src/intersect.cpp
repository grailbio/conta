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

  if (chr[0] =='M') return 0;
  if (chr[0] =='X') return 23;
  if (chr[0] =='Y') return 24;

  return atoi(chr.c_str());
}

// Compare two genomic positions.
int compare(int thisChr, int thisPos, int otherChr, int otherPos) {
  if (thisChr < otherChr) return -1;
  if (thisChr > otherChr) return 1;
  
  // chrs are the same.
  return (thisPos - otherPos);
}

// A generic VCF record for a VCF without genotype field.
class VcfRecord {
 public:
  explicit VcfRecord(std::string line) {
    this->line = line;
    std::istringstream ss(line);
    ss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info;
    numChr = numericChr(chr);

    int cafStart = info.find("CAF=");
    if (cafStart != -1) {
      int cafEnd = info.find(',', cafStart);
      // Since we skipped multi-allelic it is guaranteed there is only one MAF.
      CAF = atof(info.substr(cafStart+4, cafEnd-cafStart-4).c_str());
      MAF = 1 - CAF;
    }

    // Set whether this is an indel or multi-allelic SNP.
    indelOrmultiAllelic = alt.length() > 1;
  }

  std::string toString() {
    return (line);
  }

  bool indelOrmultiAllelic;
  int pos, numChr;
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

    std::shared_ptr<VcfRecord> vcfRecord (new VcfRecord(line));
    total++;
    getline(fin, line);

    // Skip if the new record is multi allelic.
    while (vcfRecord->indelOrmultiAllelic) {
      skipped++;
      if (line.empty() | fin.eof())
        return nullptr;
      vcfRecord.reset(new VcfRecord(line));
      total++;
      getline(fin, line);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prevRecord is null for the first record only.
    if (prevRecord != nullptr &&
        compare(vcfRecord->numChr, vcfRecord->pos,
                prevRecord->numChr, prevRecord->pos) < 0 ) {
      std::cerr << "Records out of order: " << prevRecord->toString() << "\n"
                << "larger than: " << vcfRecord->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prevRecord = vcfRecord;
    return vcfRecord;
  }

  int total = 0, skipped = 0;

 private:
  std::shared_ptr<VcfRecord> prevRecord;
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
    numChr = numericChr(chr);

    // Set whether this is a multi-allelic SNP.
    multiAllelic = alts.length() > 3;
  }

  std::string toString() {
    return (line);
  }

  std::string toTsvFormat(std::shared_ptr<VcfRecord> vcfRecord) {
    std::stringstream pp;
    if (vcfRecord != nullptr) {
      std::string alts = vcfRecord->ref + "/" + vcfRecord->alt;
      pp << chr << "\t" << pos << "\t" << A
        << "\t" << T << "\t" << G << "\t" << C << "\t" << N << "\t"
        << vcfRecord->ref << "\t" << vcfRecord->ref << "\t"
        << alts << "\t" << vcfRecord->rsid << "\t" << vcfRecord->MAF;
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
  int pos, A, T, G, C, N, numChr, INS, DEL;
  std::string ref, major;
  bool multiAllelic;
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

    pileupHeader = "CHROM\tPOS\tDEPTH\tREF\tA\tC\tG\tT\tN\tINS\tDEL";
    tsvHeader = "chrom\tpos\tA\tT\tG\tC\tN\tref\tmajor\talleles\trsid";
    tsvHeaderWithMaf = tsvHeader + "\tmaf";
    total = 0;
    skipped = 0;

    // Skip header of TSV file.
    getline(fin, header);
    if (header == pileupHeader) {
      PILEUP = TRUE;
    } else if (header == tsvHeader | header == tsvHeaderWithMaf) {
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

    std::shared_ptr<TsvRecord> tsvRecord(new TsvRecord(line, PILEUP));
    total++;
    getline(fin, line);

    // Skip if the new record is multi allelic
    while (tsvRecord->multiAllelic) {
      skipped++;
      if (line.empty() | fin.eof())
        return nullptr;
      tsvRecord.reset(new TsvRecord(line, PILEUP));
      total++;
      getline(fin, line);
    }

    // Test if the previous record was smaller than this one, otherwise file
    // is not sorted. prevRecord is null for the first record only.
    if (prevRecord != nullptr &&
        compare(tsvRecord->numChr, tsvRecord->pos,
                prevRecord->numChr, prevRecord->pos) < 0 ) {
      std::cerr << "Records out of order: " << prevRecord->toString() << "\n"
                << "larger than: " << tsvRecord->toString() << "\n";
      exit(EXIT_FAILURE);
    }

    prevRecord = tsvRecord;
    return tsvRecord;
  }

  bool PILEUP = FALSE;
  int total = 0, skipped = 0;
  std::string pileupHeader, tsvHeader, tsvHeaderWithMaf, header;

 private:
  std::string line;
  std::ifstream fin;
  std::shared_ptr<TsvRecord> prevRecord;
};

//' Intersect a TSV file with VCF file and write a TSV output
//' that contains maf values for each SNP in the TSV file.
//' In the process, remove multiple-allele SNPs if the option
//' is turned on.
//'
//' @param tsvFileName string name of input TSV file
//' @param outTsvFileName string name of output TSV file
//' @param vcfFileName string dbsnp VCF file name to be intersected
//' @param DEBUG boolean whether to print out messages
//' @return null
//'
//' @export
// [[Rcpp::export]]
void intersect(const char* tsvFileName, const char* outTsvFileName,
               const char* vcfFileName, bool nonDbSnp = false,
               bool DEBUG = false) {
    clock_t start_time, end_time;
    start_time = clock();

    if (DEBUG)
      std::cout << "\nReading TSV file from " << tsvFileName << "\n";
    TsvReader tsvReader(tsvFileName);

    // Set up out file to write
    std::ofstream outfile;
    outfile.open(outTsvFileName);
    outfile << tsvReader.tsvHeaderWithMaf << "\n";

    if (DEBUG)
      std::cout << "Reading VCF file from " << vcfFileName << "\n";
    VcfReader vcfReader(vcfFileName);

    // Get initial records from the files.
    int written = 0;
    std::shared_ptr<VcfRecord> vcfRecord = vcfReader.getNext();
    std::shared_ptr<TsvRecord> tsvRecord = tsvReader.getNext();

    if (DEBUG)
      std::cout << "Running the linear join\n";

    while (1) {
      if (tsvRecord == nullptr || vcfRecord == nullptr) {
        while (nonDbSnp && tsvRecord != nullptr) {
          outfile << tsvRecord->toTsvFormat(nullptr) << "\n";
          tsvRecord = tsvReader.getNext();
          written++;
        }
        break;
      }

      int cmp = compare(tsvRecord->numChr, tsvRecord->pos,
                        vcfRecord->numChr, vcfRecord->pos);

      if (cmp < 0) {  // VCF is ahead
        if (nonDbSnp)
          outfile << tsvRecord->toTsvFormat(nullptr) << "\n";
        tsvRecord = tsvReader.getNext();
      } else if (cmp > 0) {  // TSV is ahead
        vcfRecord = vcfReader.getNext();
      } else {  // a matching record
        outfile << tsvRecord->toTsvFormat(vcfRecord) << "\n";
        tsvRecord = tsvReader.getNext();
        vcfRecord = vcfReader.getNext();
        written++;
      }
    }

    if (DEBUG) {
      std::cout << "Read " << vcfReader.total << " VCF records, and "
                << tsvReader.total << " TSV records and wrote " << written
                << " TSV records.\nTotal of " << tsvReader.skipped
                << " TSV records were skipped (multi-allelic).\n"
                << "Total of " << vcfReader.skipped
                << " VCF records were skipped (indel or multi-allelic)\n";
    }
    outfile.close();

    end_time = clock();
    float elapsed_time = static_cast<float> (end_time - start_time);
    float elapsed_secs = elapsed_time / CLOCKS_PER_SEC;

    if (DEBUG)
      std::cout << "Elapsed seconds: " << elapsed_secs << "\n";
}
