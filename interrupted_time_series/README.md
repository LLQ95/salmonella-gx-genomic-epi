# Interrupted time series analysis

This pipeline estimates how annual resistance prevalence changed around the
Ministry of Agriculture announcement 2292. The announcement withdrew four
human and veterinary fluoroquinolones from food-producing animals at the end
of 2016, while enrofloxacin remained available under prescription.

## Data

The binary phenotype table produced for the GWAS is aggregated by sample
year. The analysis uses the continuous clinical monitoring window from 2014
through 2024, so earlier reference isolates and the single later sample are
excluded. For each year the table reports sample size and the prevalence of
QRDR mutations, gyrA, parC, qnrS and oqxAB.

## Model

A segmented regression is fitted for each determinant. The model includes
a pre-intervention time trend, a term for the immediate level change in
2017 and a term for the change in slope after the intervention. Annual
sample sizes are used as weights, and Newey-West standard errors with a
single lag are reported to account for serial correlation.

## Interpretation

The level change is the estimated one-year shift in 2017, and the slope
change is the difference between the post-intervention and pre-intervention
trends. The series has only eleven points, so estimates and p values should
be read as descriptive evidence. Because enrofloxacin was retained and
several other policy actions occurred in the same period, an observed shift
cannot be attributed to announcement 2292 alone, and results should be
reported as associations rather than policy causation.
