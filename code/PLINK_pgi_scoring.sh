#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 22             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal         # Queue (partition) name
#SBATCH -t 15:00:00       # Run time (hh:mm:ss)
#SBATCH -A OTH21060       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-type=all   # Send email at begin and end of job
#SBATCH --mail-user=peter.tanksley@austin.utexas.edu

#=======================================================================#
# Set up working environment
#=======================================================================#

export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS
export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS/bcftools-1.20

# Define input/output directories
INPUT="../input/ah_genotypes/ah_imputed_set1_1kg"
TEMP="../temp"

# Create the temp directory if it doesn't exist
mkdir -p $TEMP

# Array of chromosomes
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

#=======================================================================#
# Index VCF files
#=======================================================================#

for CHR in "${CHROMOSOMES[@]}"; do
    VCF_FILE="${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz"
    if [ ! -f "${VCF_FILE}.tbi" ] && [ ! -f "${VCF_FILE}.csi" ]; then
        bcftools index $VCF_FILE
    fi
done

#=======================================================================#
# Merge all VCF files
#=======================================================================#

# Generate a list of input VCF files
vcf_list="${TEMP}/vcf_list.txt"
> $vcf_list
for CHR in "${CHROMOSOMES[@]}"; do
    echo "${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz" >> $vcf_list
done

# Merge VCF files if the merged VCF file does not exist
merged_vcf="${TEMP}/merged.vcf.gz"
if [ ! -f $merged_vcf ]; then
    bcftools merge -l $vcf_list --force-samples -o $merged_vcf -O z
fi

#=======================================================================#
# Filter multi-allelic variants
#=======================================================================#

filtered_vcf="${TEMP}/filtered.vcf.gz"

# Remove multi-allelic variants and keep only bi-allelic SNPs if the filtered VCF file does not exist
if [ ! -f $filtered_vcf ]; then
    bcftools view -m2 -M2 -v snps $merged_vcf -o $filtered_vcf -O z
fi

#=======================================================================#
# Convert filtered VCF to PLINK2 format
#=======================================================================#

# Define the output prefix for PLINK2 conversion
plink_prefix="${TEMP}/ah_biallelic"

# Check if PLINK2 output files exist before running the conversion
if [ ! -f "${plink_prefix}.pgen" ] || [ ! -f "${plink_prefix}.pvar" ] || [ ! -f "${plink_prefix}.psam" ]; then
    plink2 --vcf $filtered_vcf --make-pgen --out $plink_prefix
fi

#=======================================================================#
# Filter to ACGT alleles using PLINK2
#=======================================================================#

# Define the output prefix for the final PLINK2 files with ACGT alleles
plink_acgt_prefix="${TEMP}/ah_biallelic_acgt"

# Check if the final PLINK2 output files exist before running the ACGT allele filter
if [ ! -f "${plink_acgt_prefix}.pgen" ] || [ ! -f "${plink_acgt_prefix}.pvar" ] || [ ! -f "${plink_acgt_prefix}.psam" ]; then
    plink2 --pfile $plink_prefix --snps-only just-acgt --make-pgen --out $plink_acgt_prefix
fi

#=======================================================================#
# Clean up intermediate files
#=======================================================================#

# Uncomment the line below if you want to remove intermediate files
# rm $filtered_vcf $merged_vcf $vcf_list

#=======================================================================#
# Done
#=======================================================================#

