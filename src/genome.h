#ifndef ANALYSIS_CONTA_SRC_GENOME_H_
#define ANALYSIS_CONTA_SRC_GENOME_H_

// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

#include <Rcpp.h>
#include <map>
#include <iostream>
#include <string>

// Genome.h
class Genome {
 public:
    Genome(std::string file, bool DEBUG);
    std::string getPos(std::string contig, int pos, int length);
    std::map<std::string, std::string> genome;
};

#endif  // ANALYSIS_CONTA_SRC_GENOME_H_
