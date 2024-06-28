#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal         # Queue (partition) name
#SBATCH -t 24:00:00       # Run time (hh:mm:ss)
#SBATCH -A OTH21060       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-type=all   # Send email at begin and end of job
#SBATCH --mail-user=peter.tanksley@austin.utexas.edu

#=======================================================================#
# Set up working environment
#=======================================================================#

export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS

# Define input/output directories
INPUT="../input/ah_genotypes/ah_imputed_set1_1kg"
TEMP="../temp"

# Create the temp directory if it doesn't exist
mkdir -p $TEMP

# Array of chromosomes
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

#=======================================================================#
# Convert VCF files to PLINK binary format
#=======================================================================#

# Create merge list for PLINK
vcf_list="${TEMP}/vcf_list.txt"
> "$vcf_list"
for CHR in "${CHROMOSOMES[@]}"; do
    plink --vcf "${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz" --make-bed --out "${TEMP}/chr${CHR}"
    echo "${TEMP}/chr${CHR}" >> "$vcf_list"
done

#=======================================================================#
# Merge all binary files using PLINK (v1.9)
#=======================================================================#

# Merge using merge list
plink --merge-list "$vcf_list" --make-bed --out "${TEMP}/merged"

#=======================================================================#
# Filter to biallelic SNPs and convert to PLINK2 format with ACGT alleles using PLINK2
#=======================================================================#

# Filter to biallelic SNPs, convert to PLINK2 format, and filter to ACGT alleles
plink2 --bfile "${TEMP}/merged" --biallelic-only strict --snps-only just-acgt --make-pgen --out "${TEMP}/ah_biallelic_acgt"

