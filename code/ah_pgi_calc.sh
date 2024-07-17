#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal         # Queue (partition) name
#SBATCH -t 01:00:00       # Run time (hh:mm:ss)
#SBATCH -A OTH21060       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-type=all   # Send email at begin and end of job
#SBATCH --mail-user=peter.tanksley@austin.utexas.edu

#=======================================================================#
# Set up working environment
#=======================================================================#

export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS
export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS/bcftools-1.20

INPUT="../input/phg001099.v1.AddHealth.genotype-imputed-data.c1.GRU-IRB-PUB-GSO.set1"
TEMP="../temp"
SNP_WEIGHTS="../input/snp_weights"

# Create the temp directory if it doesn't exist
mkdir -p $TEMP

# Define the list of chromosomes
CHROMOSOMES=$(seq 1 22)

#=======================================================================#
# Step 1: Convert VCF files to binary format
#=======================================================================#

for CHR in $CHROMOSOMES; do
    if [ ! -f ${TEMP}/chr${CHR}.pgen ]; then
        plink2 --vcf ${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz --make-pgen --out ${TEMP}/chr${CHR}
    fi
done

#=======================================================================#
# Step 2: Filter for biallelic SNPs with A, C, G, T alleles
#=======================================================================#

for CHR in $CHROMOSOMES; do
    if [ ! -f ${TEMP}/chr${CHR}_snps.pgen ]; then
        plink2 --pfile ${TEMP}/chr${CHR} --snps-only just-acgt --make-pgen --out ${TEMP}/chr${CHR}_snps
    fi
done

#=======================================================================#
# Step 3: Fix IDs using chromosome, position, ref, alt (match score file)
#=======================================================================#

for CHR in $CHROMOSOMES; do
    if [ ! -f ${TEMP}/chr${CHR}_snps_allids.pgen ]; then
        plink2 --pfile ${TEMP}/chr${CHR}_snps --set-all-var-ids @:#\$r,\$a --new-id-max-allele-len 23 missing --make-pgen --out ${TEMP}/chr${CHR}_snps_allids
    fi
done

#=======================================================================#
# Step 4: remove duplicates
#=======================================================================#

for CHR in $CHROMOSOMES; do
    if [ ! -f ${TEMP}/chr${CHR}_snps_allids_nodups.pgen ]; then
        plink2 --pfile ${TEMP}/chr${CHR}_snps_allids --rm-dup force-first --make-pgen --out ${TEMP}/chr${CHR}_snps_allids_nodups
    fi
done

if [ ! -f "${TEMP}/ah_snps_allids_nodups_merged.pgen" ]; then

	# Create empty merge_list file for PLINK
	merge_list="${TEMP}/merge_list.txt"
	> "$merge_list"

	# generate merge list
	for CHR in $CHROMOSOMES; do
		pfile_prefix="${TEMP}/chr${CHR}_snps_allids_nodups"
	echo "$pfile_prefix" >> "$merge_list"
	done

	#merge
	plink2 --pmerge-list "$merge_list" --multiallelics-already-joined --make-pgen --out "${TEMP}/ah_snps_allids_nodups_merged"
fi

#=======================================================================#
# Step 4: Compute scores using --score with no-sum modifier
#=======================================================================#


# List of SNP weight files by Big5 factor
factors=("agr" "con" "ext" "neu" "ope")

# Loop through each item
for factor in "${factors[@]}"; do
    # Define the output file path
    output_file="${TEMP}/ah_${factor}_pgi"

    # Check if the output file exists
    if [ ! -f "${output_file}.sscore" ]; then
        plink2 --pfile "${TEMP}/ah_snps_allids_nodups_merged" --score "${SNP_WEIGHTS}/all_${factor}.PRS_input.snpRes" 12 5 8 header-read center --out "${output_file}"
    fi
done


