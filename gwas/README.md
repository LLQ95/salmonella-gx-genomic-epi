# Fluoroquinolone GWAS with pyseer

This pipeline runs a genome-wide association analysis for fluoroquinolone
resistance while controlling for population structure. The analysis is
designed to run on the biflab-1 server.

## Inputs

The joint core alignment and VCF provide the variant data. PointFinder
supplies chromosomal QRDR mutations in gyrA, gyrB, parC and parE, and
ResFinder supplies plasmid-mediated qnrS and oqxAB carriage. The unified
metadata JSON defines the sample set.

## Phenotypes

Five binary phenotype columns are generated. The qrdr column marks any
chromosomal QRDR mutation, while gyrA and parC separate the two loci. The
qnrS1 and oqxAB columns mark plasmid-mediated gene carriage. Samples
without a detected mutation or gene are coded as zero.

## Population structure

A pairwise core-genome similarity matrix is computed from the alignment.
pyseer uses this matrix with multidimensional scaling to account for
clonal relatedness, which limits spurious associations driven by lineage
structure.

## Running

Set WORKDIR to the joint analysis directory and run the shell driver. The
script checks for pyseer, prepares the phenotype table and similarity
matrix, and runs pyseer separately for each phenotype. Each output table
contains effect estimates, standard errors and p values, with a
Benjamini-Hochberg adjusted column for multiple testing.

## Interpretation

The QRDR analysis addresses chromosomal resistance determinants, whereas
the qnrS and oqxAB analyses address plasmid-mediated elements that on
their own usually confer reduced susceptibility rather than full clinical
resistance. Results should be reported as associations with adjusted p
values and effect sizes, not as causal evidence.
