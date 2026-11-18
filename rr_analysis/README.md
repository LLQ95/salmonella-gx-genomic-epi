# RR and novel RR analysis across Guangxi, adjacent provinces, and Southeast Asia

This pipeline tests whether the long-read genomes contain a resistance region

(RR) outside RR1 to RR8 as described in Wang et al. 2023, and extends the

same gene-content scan to every short-read assembly in the joint tree.

## Inputs

The default paths match biflab-1 and can be overridden through environment

variables.



| Variable      | Default                                       |
| ------------- | --------------------------------------------- |
| WORKDIR       | /data/lilanqi/salmonella/joint\_gx\_sea\_tree |
| LONGREAD\_DIR | \${WORKDIR}/longreads                         |
| FASTA\_DIR    | \${WORKDIR}/input\_fasta                      |
| TRAIT\_TABLE  | \${WORKDIR}/trait\_table.tsv                  |
| OUT\_BASE     | \${WORKDIR}/rr\_novel\_analysis               |

The trait table must contain node, dataset, and region columns. The probes

directory is self-contained and includes LT2 boundary genes, robust anchor

genes, the Tn21 mer probe, the ResFinder reference sequences, and the finished

RR3 reference KR856283 with its structural probes.

## Reference material

download\_references.py fetches the RR3 GenBank record KR856283 through NCBI

efetch and is idempotent, so an existing file is left untouched. The script

also lists the Wang et al. 2023 representative runs. PacBio RS II runs from

2013 are stored as tar.gz archives of legacy HDF5 files and are skipped unless

\--with-pacbio is given. Illumina runs for RR2 and RR8 can be added with

\--with-illumina and assembled separately. prepare\_structure\_probes.py extracts

the IS26, IS1, IS10, MAK1, and meth markers from the RR3 reference.

## Running

The pipeline needs blastn, BioPython, and pandas. Activate an environment that

provides them, then run:



```
bash run_rr_analysis.sh
```

## Steps

Step 1 extracts one continuous fljAB locus from each long-read genome using a

left and right anchor chain. Step 2 types the resolved loci, recognizes the

known RR3, RR5, RR6, RR8, and RR4/RR7 cargo patterns, and flags only

unmatched patterns as RR\_variant. Step 5 adds structural evidence from the

RR3 reference, counts IS26 copies, and separates RR4 from RR7 because RR4

keeps the MAK1 and meth region while RR7 has lost it. Step 3 scans every

short-read assembly for the same boundary genes and resistance cargo without

assuming one continuous locus. Step 4 summarizes cargo patterns by dataset and

reports short-read strains that share a long-read-confirmed novel pattern.

## Outputs

All outputs are written under OUT\_BASE. longread\_rr\_calls.csv gives the

cargo-based calls, and longread\_precise\_calls.csv gives the precise structural

calls with the left boundary, IS26 count, and MAK1 and meth markers.

shortread\_rr\_scan.csv gives the strain-level gene content for the full cohort.

region\_summary contains region\_cargo\_patterns.csv and

novel\_rr\_shortread\_candidates.csv.

## Resolution limit

The adjacent-province and Southeast-ASEA strains are short-read assemblies

with fragmented contigs, and the fljAB repeat region is often split across

contigs. Their matches are reported as short-read candidates and cannot be

confirmed as novel RR without long-read structural evidence.