# Major analysis changes, September 2026

This note summarizes the main analysis and metadata changes made during the final manuscript cleanup. It is intended as a quick orientation document for collaborators and for checking that the manuscript, supporting tables, online data tables, and repository metadata all describe the same analysis state.

## 1. pixy all-sites dataset was rebuilt

The original pixy all-sites VCF contained too few invariant sites because invariant sites had been retained with a strict `QUAL >= 30` filter. This produced only about 36.5 million invariant records, which was much lower than expected for an autosomal all-sites dataset.

We rebuilt the invariant-site VCF from chromosome-level GATK all-sites VCFs without filtering invariant sites by `QUAL`. Genotypes with depth below 5 were set to missing before pixy missing-data accounting. The final all-sites VCF used for pixy contains:

- invariant records: 1,145,122,611
- variant records: 1,783,056
- total records before pixy missing-data accounting: 1,146,905,667
- callable sites after pixy missing-data accounting: about 1.10 billion

This corrected the denominator for diversity analyses and substantially lowered the genome-wide estimates of nucleotide diversity compared with the earlier run.

## 2. pixy was rerun without a MAF filter for pi and Tajima's D

The final diversity analyses use a no-MAF SNP set together with invariant sites. This was done because nucleotide diversity and especially Tajima's D should not be estimated from a dataset where rare variants have been removed.

The final pixy run therefore used:

- invariant autosomal sites, no `QUAL` filter, genotype depth filter `DP >= 5`
- SNPs after genotype-level filtering and missingness filtering, but without MAF filtering
- the curated 50-sample set described below

The final genome-wide estimates are:

- weighted mean pi: 0.000247
- weighted mean Tajima's D: -0.699
- macrochromosome weighted mean pi: 0.000246
- microchromosome weighted mean pi: 0.000249
- macrochromosome weighted mean Tajima's D: -0.692
- microchromosome weighted mean Tajima's D: -0.721

## 3. Final pixy sample set was changed to curated50

The final pixy analyses use 50 curated individuals. The sample set excludes high-missingness/problem samples and close relatives, retaining FNM12 as the representative individual from the related Fennoscandia group.

Final pixy sample counts:

- FN: 13
- GR: 8
- NY: 13
- SK: 11
- WR: 5
- total: 50

Important sample decisions:

- retained: FNM12
- excluded from pixy: FNM11, FNM15, WRG01, GRL04, SKW02
- additional exclusions follow the curated QC rule in `data/metadata/samples_curated.csv`

The `include_pixy` column in the sample metadata and supporting Table S2 was updated to match this final curated50 dataset.

## 4. Singleton diagnostic confirmed the curated set was needed

SUPER_10 sensitivity analyses showed that the all-sample no-MAF dataset gave strongly negative Tajima's D, driven by problematic singleton patterns. Removing WRG01 had a large effect, and the curated50 dataset gave stable estimates.

SUPER_10 sensitivity summary:

- all66 no-MAF: pi = 0.000267, Tajima's D = -1.846
- noWRG01_65 no-MAF: pi = 0.000238, Tajima's D = -0.754
- curated50 no-MAF: pi = 0.000236, Tajima's D = -0.715

Within the curated set, singleton heterozygotes and common heterozygotes had similar allele-balance and genotype-quality distributions, supporting the final curated dataset.

## 5. FST and dXY were rerun from the same curated50 pixy output

FST and dXY were also rerun using the final curated50 pixy all-sites dataset. Genome-wide pairwise FST remains effectively zero among sampling regions. Negative mean FST estimates are interpreted as zero in figures and text.

The Figure 1B heatmap was updated from pixy-derived pairwise FST values, with population order:

FN, GR, NY, SK, WR

The heatmap now displays only the upper triangle and uses negative estimates as zero for visualization.

## 6. FST outlier analyses were updated

FST outlier tables and plots were regenerated from the curated50 pixy results.

Current interpretation:

- no convincing genome-wide peaks indicating local adaptation
- 187 windows were high-FST outliers in at least one pairwise comparison
- 35 windows were outliers in more than one comparison
- many high-FST windows involve WR comparisons, which should be interpreted cautiously because WR has only five individuals in the curated pixy set
- recurrent/high-FST windows often overlap repeats, especially LINE/CR1 elements, and 122 overlap annotated genes

The outlier results should be described as exploratory rather than evidence for local adaptation.

## 7. Centromere and pericentromere diversity summaries were updated

Centromere/pericentromere analyses were rerun with the corrected curated50 no-MAF pixy results.

Final weighted pi estimates:

- centromere: 0.000071
- pericentromere: 0.000251
- other/background: 0.000250

Final weighted Tajima's D estimates:

- centromere: -0.749
- pericentromere: -0.813
- other/background: -0.696

These results support lower diversity in annotated putative centromeric windows, but interpretation remains cautious because centromere annotations are available only for a subset of chromosomes and these regions may have mapping/callability biases.

## 8. Mean FST versus centromere regions was regenerated

The mean-FST-by-centromere-region plot had not been updated in the first rerun. It was regenerated from the curated50 FST results.

For this plot, mean pairwise FST across comparisons was calculated per window, and negative values were set to zero for display.

Summary:

- centromere mean FST: 0.00502
- pericentromere mean FST: 0.00422
- other/background mean FST: 0.00424
- Kruskal-Wallis p = 0.969

There is no evidence that mean FST differs among centromere, pericentromere, and background windows.

## 9. Supporting tables updated

The supporting tables with actual content changes are:

- Table S2: sequencing/genotype QC and analysis inclusion
- Table S5: pairwise FST and dXY among sampling regions
- Table S6: genome-wide diversity by macro- and microchromosomes
- Table S7: centromeric, pericentromeric, and background diversity

Table S1 was retitled as "Sample metadata overview" because the analysis inclusion flags are in Table S2.

## 10. Online data table 1 and pixy metadata were synchronized

The full sample metadata table was updated so that the `include_pixy` column matches the final curated50 pixy dataset.

The pixy population files were also updated:

- `data/metadata/pixy/population_all.txt`
- `data/metadata/pixy/population_regions.txt`

Both now contain 50 samples.

## 11. Main interpretation after the corrections

The main biological interpretation is broadly unchanged for population structure: snowy owls show very little genome-wide population differentiation across sampled regions.

The diversity/genomic-health interpretation changed more substantially. With the corrected all-sites denominator and no-MAF diversity analysis, genome-wide pi is lower than previously estimated. This makes the results more consistent with a history of past and recent bottlenecks, while the ROH results still suggest little evidence for very recent severe inbreeding in the sampled individuals.

