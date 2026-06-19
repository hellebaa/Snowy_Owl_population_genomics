Snowy Owl Diversity Analyses Using pixy v2.0.0

========================================
OVERVIEW
========

This directory contains analyses of genome-wide genetic diversity
in Snowy Owls using pixy v2.0.0.

pixy was chosen because it provides unbiased estimates of population
genetic summary statistics from datasets containing missing data and
requires both variant and invariant callable sites.

Input VCF:

* SnowyOwl_autosomes_allSites_noContam.vcf.gz

This VCF contains:

* 1,324,912 filtered autosomal SNPs
* 36,529,297 filtered callable invariant sites

for a total of:

* 37,854,209 callable autosomal sites

Two contaminated individuals were removed prior to analysis:

* GRL04_merged_X14
* SKW02_merged_X22

All analyses were performed in non-overlapping 50 kb windows.

========================================
ANALYSIS 1: ALL INDIVIDUALS
===========================

Directory:

* 01_pixy_all/

Population file:

* population_all.txt

This analysis treats all Snowy Owl samples as belonging to a single
population.

Purpose:

* Estimate species-wide nucleotide diversity (π)
* Estimate species-wide Tajima's D
* Characterize overall levels of genetic variation across the genome

This analysis is appropriate because previous analyses indicated very
weak population structure across the sampled range.

Statistics calculated:

* π
* Tajima's D

Outputs:

* SnowyOwl_ALL_pi.txt
* SnowyOwl_ALL_tajima_d.txt

These estimates represent genome-wide diversity for Snowy Owls as a
single population.

========================================
ANALYSIS 2: REGIONAL POPULATIONS
================================

Directory:

* 02_pixy_regions/

Population file:

* population_regions.txt

Individuals were assigned to populations based on sampling location
prefixes in sample names (e.g. FNM, GRL, WRG, SKW).

Purpose:

* Compare diversity among regions
* Test for geographic population structure
* Quantify differentiation among regions

Statistics calculated:

* π
* Tajima's D
* FST
* dXY

Outputs:

* SnowyOwl_regions_pi.txt
* SnowyOwl_regions_tajima_d.txt
* SnowyOwl_regions_fst.txt
* SnowyOwl_regions_dxy.txt

Interpretation:

π:

* Within-population genetic diversity

Tajima's D:

* Deviations from neutral equilibrium expectations
* Sensitive to demographic history and selection

FST:

* Relative genetic differentiation among regions
* Values near zero indicate little population structure

dXY:

* Absolute sequence divergence among regions
* Less affected by within-population diversity than FST

========================================

RATIONALE FOR RUNNING BOTH ANALYSES
===================================

The ALL analysis provides the best estimate of species-wide diversity
because all individuals contribute to the calculation.

The REGIONAL analysis addresses a different question:

"Do different parts of the Snowy Owl distribution differ in genetic
diversity or show evidence of population structure?"

Together, these analyses allow estimation of both:

1. Overall species-wide diversity
2. Geographic variation in diversity and differentiation

========================================

SUMMARY OF RESULTS
==================

All individuals, Genome-wide estimates:

* Mean π = 0.00721304
* Mean Tajima's D = 0.220997

Regional differentiation:

* Mean FST = 0.000238
* Mean dXY = 0.007217

These values indicate high genetic diversity and extremely weak
population structure across the sampled range.

========================================

SOFTWARE
========

pixy v2.0.0

Reference:

Korunes KL & Samuk K (2021)
pixy: Unbiased estimation of nucleotide diversity and divergence in
the presence of missing data.
Molecular Ecology Resources 21:1359–1368.
