// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

#include <Rcpp.h>
#include "genome.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace Rcpp;

Genome::Genome(std::string file, bool DEBUG = false) {

  clock_t start_time, end_time;
  start_time = clock();
  std::string line;
  std::ifstream fin;
  fin.open(file.c_str());
  if (!fin) {
    std::cerr << "Could not read file" << file <<  "\n";
    exit(EXIT_FAILURE);
  }

  std::string seq = "";
  std::string name = "";
  while (!fin.eof()) {
    if (line.substr(0, 1) == ">") {
      if (seq.length() > 0) {
        genome[name] = seq;
      }
      seq = "";
      name = line.substr(1, line.length());

      // Remove chr from beginning of fasta id
      if (name.find("chr") == 0 || name.find("Chr") == 0 ||
          name.find("CHR") == 0)
        name = name.substr(3);

      getline(fin, line);
    } else {
        seq += line;
        getline(fin, line);
    }
  }
  genome[name] = seq;

  end_time = clock();
  float elapsed_time = static_cast<float> (end_time - start_time);
  float elapsed_secs = elapsed_time / CLOCKS_PER_SEC;

  if (DEBUG)
    std::cout << "Read genome in: " << elapsed_secs << " seconds\n";
}

std::string Genome::getPos(std::string contig, int pos, int length) {
  if (pos >= 0 && genome[contig].length() > pos) {
    return genome[contig].substr(pos, length);
  } else {
    return std::string(length, 'N');
  }
}

