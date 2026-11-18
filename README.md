# Genomic epidemiology of Salmonella Typhimurium and its monophasic variant

This repository contains the analysis code and derived data tables used to
investigate resistance-region structure, fluoroquinolone resistance
associations and policy-period trends in the joint genomic collection.

## Contents

The rr_analysis directory extracts continuous long-read loci, classifies
the resistance regions (RR) and scans the short-read assemblies. It maps
the data-driven structures to the RR1 through RR8 nomenclature and includes
the reference and structure probes.

The gwas directory builds binary phenotypes from PointFinder and
ResFinder, computes a core-genome similarity matrix and runs pyseer with
multidimensional scaling to control population structure.

The interrupted_time_series directory aggregates annual prevalence and
fits a segmented regression around the 2017 intervention, with
Newey-West standard errors.

The linkage directory contains the pairwise core-SNP comparison between
Guangxi isolates and the reference collections.

The data directory holds the derived long-read typing table after
exclusion, the GWAS phenotype table, the annual series and the segmented
regression coefficients.

## Reproducibility

The shell drivers use environment variables for all input paths and run
dependency checks. Raw sequencing reads and the full reference databases
are not included. Statistical associations are reported with adjusted p
values and effect estimates, and observational results are not interpreted
as causal evidence.
