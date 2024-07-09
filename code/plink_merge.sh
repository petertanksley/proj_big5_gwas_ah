#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal         # Queue (partition) name
#SBATCH -t 48:00:00       # Run time (hh:mm:ss)
#SBATCH -A OTH21060       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-type=all   # Send email at begin and end of job
#SBATCH --mail-user=peter.tanksley@austin.utexas.edu

#=======================================================================#
# Set up working environment
#=======================================================================#

export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS
export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS/bcftools-1.20

# Define input/output directories
INPUT="../input/phg001099.v1.AddHealth.genotype-imputed-data.c1.GRU-IRB-PUB-GSO.set1"
TEMP="../temp"

# Create the temp directory if it doesn't exist
mkdir -p $TEMP

# Array of chromosomes
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

#=======================================================================#
# Normalize VCF files to join multiallelic variants
#=======================================================================#

# Normalize VCF files
for CHR in "${CHROMOSOMES[@]}"; do
    vcf_file="${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz"
    normalized_vcf_file="${TEMP}/chr${CHR}.normalized.vcf.gz"

    if [ ! -f "$normalized_vcf_file" ]; then
        bcftools norm -m -both "$vcf_file" -O z -o "$normalized_vcf_file"
    fi
done

#=======================================================================#
# Convert normalized VCF files to PLINK binary format
#=======================================================================#

# Create pgen_list file for PLINK
pgen_list="${TEMP}/pgen_list.txt"
> "$pgen_list"

# Convert normalized VCF files
for CHR in "${CHROMOSOMES[@]}"; do
    normalized_vcf_file="${TEMP}/chr${CHR}.normalized.vcf.gz"
    pgen_file="${TEMP}/chr${CHR}.pgen"

    if [ ! -f "$pgen_file" ]; then
        plink2 --vcf "$normalized_vcf_file" --make-pgen --out "${TEMP}/chr${CHR}"
    fi
    echo "${TEMP}/chr${CHR}" >> "$pgen_list"
done

#=======================================================================#
# Merge PLINK binary files
#=======================================================================#

plink2 --pmerge-list "$pgen_list" --make-pgen --merge-max-allele-ct 2 --out "${TEMP}/merged"

#=======================================================================#
# Filter to biallelic SNPs and convert to PLINK2 format with ACGT alleles using PLINK2
#=======================================================================#

# Filter to biallelic SNPs, convert to PLINK2 format, and filter to ACGT alleles
plink2 --pfile "${TEMP}/merged" --max-alleles 2 --snps-only just-acgt --make-pgen --out "${TEMP}/ah_biallelic_acgt"

