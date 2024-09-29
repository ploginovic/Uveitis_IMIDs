#!/bin/bash

# Define base directories for BGEN and SAMPLE files
bgen_base="/slade/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/imputed_data"
sample_file="${bgen_base}/ukb9072_imp_autosomes.sample"
score_file="/slade/home/pl450/Uveitis/GRS/Psoriasis/Psoriasis_109significant_sumstats.tsv"

# Loop over all chromosomes (1-22)
for chr in {1..22}
do

  # BGEN file for the current chromosome
  bgen_file="${bgen_base}/ukb_imp_chr${chr}_v3.bgen"

  # Output file prefix
  output_prefix="GRS_0207/Psoriasis_GRS_chr${chr}_020724"

  # Run PLINK2 for scoring
  plink2 --bgen $bgen_file ref-first \
         --sample $sample_file \
         --score $score_file 1 4 6 header cols=+scoresums \
         --out $output_prefix

  echo "Scoring completed for chromosome ${chr}"
done

echo "All scoring operations completed."