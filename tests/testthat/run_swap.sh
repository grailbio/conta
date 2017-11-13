#!/bin/sh

map_table_name="s3://grail-proc-monitoring/analysis/data/CCGA_BST/mapping_cfdna_n_gdna_20170831.csv"
out_table_name="s3://grail-proc-monitoring/analysis/data/CCGA_BST/mapping_cfdna_n_gdna_20170831_swap.csv"
shiny_loc="s3://grail-proc-monitoring/data/shiny_data.csv"
dbsnp_targeted="s3://grail-onur/dbsnp/dbSNP_common_panV2overlap_dupRemoved_multiRemoved.vcf"

Rscript /grail/src/grail/analysis/conta/scripts/swap_run.R -m $map_table_name \
  -o $out_table_name -s $shiny_loc -d $dbsnp_targeted
