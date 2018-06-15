#!/bin/sh

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

Rscript ../scripts/conta_run.R --dbSNP_file dbSNP_common_panV2overlap_dupRemoved_multiRemoved_20170613.vcf \
	--lr_th 0.01 --min_cf 5e-05 --sample test --save_dir test --tsv_file \
	test.collapsed.nobaq.tsv --min_maf 0.1
