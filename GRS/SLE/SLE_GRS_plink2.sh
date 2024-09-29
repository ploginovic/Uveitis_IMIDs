#!/bin/bash

# Define base directories for BGEN and SAMPLE files
bgen_base="/slade/projects/Research_Project-MRC158833/UKBiobank/500K_Genetic_data/imputed_data"
sample_file="${bgen_base}/ukb9072_imp_autosomes.sample"
score_file="/slade/home/pl450/Uveitis/GRS/SLE/filtered_SLE_QC_clumped"

# Loop over all chromosomes (1-22) excluding 14, 15, and 21
for chr in {1..22}
do
  # Skip chromosomes 14, 15, and 21
  if [ "$chr" -eq 14 ] || [ "$chr" -eq 15 ] || [ "$chr" -eq 21 ]; then
    echo "Skipping chromosome $chr"
    continue
  fi

  # BGEN file for the current chromosome
  bgen_file="${bgen_base}/ukb_imp_chr${chr}_v3.bgen"

  # Output file prefix
  output_prefix="GRS_0705/SLE_GRS_chr${chr}_070524"

  # Run PLINK2 for scoring
  plink2 --bgen $bgen_file ref-first \
         --sample $sample_file \
         --score $score_file 3 5 7 header cols=+scoresums \
         --out $output_prefix

  echo "Scoring completed for chromosome ${chr}"
done

echo "All scoring operations completed."
