#!/usr/bin/env python3

# Author: Andy Tock
# Date: 19/10/2023

# Usage:
# conda activate python_3.11
# ./parse_vcf_get_TopEffect.py --prefix clinvar_20231015.ann
# conda deactivate

import os
import sys
import glob
import argparse
import re
import pandas as pd
import numpy as np
import gc

def create_parser():
    """
    Capture input command-line arguments
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--prefix", type=str, default="clinvar_20231015.ann",
        help="VCF file prefix. Default: clinvar_20231015.ann")
    return parser

parser = create_parser().parse_args()
print(parser)

# Get symbols for genes in panel
bed = pd.read_csv(
    "/Users/ajtock/dnanexus/snakemake_vep_plugins/resources/Origin_Assets_EU/Nonacus_target_merged.bed",
    header=None,
    sep="\t")
gene_list = sorted( list( set( 
    [re.sub("_.+", "", x) for x in list(bed[3])]
) ) )
print(gene_list)
print("Number of genes in gene_list: %s" % len(gene_list))
# NOTE: No single-gene variants in HRAS found on ClinVar

# Parse SnpEff-annotated ClinVar VCF and write to TSV
def parse_vcf(prefix, info_cols):
    """
    Load and concatenate VCF files from specified directories.
    :param prefix: VCF file prefix
    :param info_cols: dictionary of field:dtype for INFO fields to
    store as distinct columns.
    :return: DataFrame containing SnpEff-annotated ClinVar VCF information
    """
    #prefix=parser.prefix
    #info_cols={
    #    "ALLELEID": "int",
    #    "CLNDISDB": "str",
    #    "CLNDN": "str",
    #    "CLNHGVS": "str",
    #    "CLNREVSTAT": "str",
    #    "CLNSIG": "str",
    #    "CLNVC": "str",
    #    "CLNVCSO": "str",
    #    "GENEINFO": "str",
    #    "MC": "str",
    #    "ORIGIN": "str",
    #    "ANN": "str"
    #}
    ## End args
    header = "CHROM POS ID REF ALT QUAL FILTER INFO".split()
    # Load VCF data from subdir_list
    vcf = pd.read_csv(
        prefix + ".vcf",
        sep="\t",
        comment="#",
        names=header,
        dtype={
            "CHROM": str})

    # Convert INFO into a dictionary of subfields
    vcf["INFO"] = vcf["INFO"].str.split(";") \
        .apply(lambda x: dict([y.split("=") for y in x]))

    # Add given INFO subfields as columns
    if info_cols is not None:
        for field, dtype in info_cols.items():
            vcf[field] = vcf["INFO"].apply(lambda x: x.get(field, None))
            vcf[field] = vcf[field].astype(dtype)

    # Convert ANN into a list in which each element is one annotation containing multiple subfields
    vcf["ANN_list"] = vcf["ANN"].str.split(",")

    # Convert each ANN_list element into a list of subfields
    vcf["ANN_list_subfields"] = vcf["ANN_list"] \
        .apply(lambda x: [y.split("|") for y in x])
    
    # Get "TopEffect" annotations as in TwinStrand pipeline by selecting the first reported effect
    # On SnpEff effect sort order, see http://pcingola.github.io/SnpEff/snpeff/inputoutput/
    vcf["ANN_TopEffect"] = vcf["ANN_list_subfields"] \
        .apply(lambda x: x[0])
    topeffect_cols = [
        "TopEffectAllele",
        "TopEffectAnnotation",
        "TopEffectImpact",
        "TopEffectGeneName",
        "TopEffectGeneId",
        "TopEffectFeatureType",
        "TopEffectFeatureId",
        "TopEffectTranscriptBiotype",
        "TopEffectExonIntronRankTotal",
        "TopEffectHGVSc",
        "TopEffectHGVSp",
        "TopEffectcDNAposcDNAlen",
        "TopEffectCDSposCDSlen",
        "TopEffectProteinposProteinlen",
        "TopEffectDistanceToFeature",
        "TopEffectInfo"
    ]
    for y in range(0, len(topeffect_cols)):
        vcf[topeffect_cols[y]] = vcf["ANN_TopEffect"]. \
            apply(lambda x: x[y] if len(x) == len(topeffect_cols) else '')

    # Remove variants without SnpEff annotations
    vcf = vcf.loc[~(vcf["ANN"] == "None")]

    # Remove duplicate rows based on CHROM, POS, REF and ALT
    vcf.drop_duplicates(
        subset=["CHROM", "POS", "REF", "ALT"],
        keep="first",
        inplace=True,
        ignore_index=False)

    # Remove unnecessary columns
    vcf.drop(
        columns=["INFO", "ANN", "ANN_list", "ANN_list_subfields", "ANN_TopEffect"],
        inplace=True)

    # Get variants affecting genes in panel
    vcf_genes = vcf.loc[
        (vcf["TopEffectGeneName"].str.fullmatch("|".join(gene_list)))]
    
    # Write to TSV
    vcf.to_csv(
        prefix + "_parsed_genomewide.tsv",
        na_rep="NaN", sep="\t", header=True, index=False)
    vcf_genes.to_csv(
        prefix + "_parsed_50genepanel.tsv",
        na_rep="NaN", sep="\t", header=True, index=False)

    return None


def main():
    # Parse SnpEff-annotated ClinVar VCF and write to TSV
    parse_vcf(
        prefix=parser.prefix,
        info_cols={
            "ALLELEID": "int",
            "CLNDISDB": "str",
            "CLNDN": "str",
            "CLNHGVS": "str",
            "CLNREVSTAT": "str",
            "CLNSIG": "str",
            "CLNVC": "str",
            "CLNVCSO": "str",
            "GENEINFO": "str",
            "MC": "str",
            "ORIGIN": "str",
            "ANN": "str"
        })


if __name__ == "__main__":
    main()
