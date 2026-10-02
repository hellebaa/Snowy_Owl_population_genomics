SnowyOwl_autosomes_allSites_curated50_noQUAL_DP5_noMAF.vcf.gz

Contains:
- autosomes only
- 50 curated high-quality individuals
- high-missingness, contamination-review/remove, and close-related samples
  removed, except FNM12 retained as the Fennoscandia related-group representative
- SNP sites from:
  SnowyOwl_autosomes_snps_curated50_missing20_noMAF.recalc.vcf.gz
- corrected invariant sites from:
  SnowyOwl_autosomes_invariant_curated50_noQUAL_DP5.vcf.gz

Invariant sites were retained as ALT=. records from the GATK all-sites VCF.
No site-level QUAL filter was applied to invariant sites. Individual invariant
genotypes with DP < 5 were set to missing before running pixy.

Record counts:
- autosomal SNPs retained without MAF filtering: 1,783,056
- autosomal invariant-site records: 1,145,122,611
- autosomal all-sites records: 1,146,905,667
- pixy ALL callable sites after missing-data accounting: 1,101,113,709

Used for pixy nucleotide diversity, Tajima's D, FST, and dXY.
