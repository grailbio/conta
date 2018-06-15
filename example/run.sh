#!/bin/sh

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

dbSNP_file=common_all_20180423.vcf
dbSNP_ftp=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
if [ ! -f "$dbSNP_file" ]
then
	wget "$dbSNP_ftp"
	gunzip "$dbSNP_file.gz"
fi

for tsv_file in 170501_cfDNABRH1308323_Conta_*.tsv
do
	name="${tsv_file%.tsv}"
	echo Running "$name"
	Rscript ../scripts/conta_run.R --dbSNP_file "$dbSNP_file" --min_cf 5e-05 --sample "$name" \
	--save_dir "$name" --tsv_file "$tsv_file" --min_maf 0.01 --cf_correction 0
done
