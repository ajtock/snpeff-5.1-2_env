#!/bin/zsh

# Usage:
# ./snpeff_annotate_clinvar.sh clinvar_20231015 > clinvar_20231015.ann.log 2>&1 &
# echo $! > clinvar_20231015.ann.pid

PREFIX=$1

source ~/.zshrc
conda activate snpeff

snpeff \
  -nodownload \
  -verbose \
  GRCh38.mane.1.0.ensembl \
  ${PREFIX}.vcf.gz \
  > ${PREFIX}.ann.vcf

conda deactivate
