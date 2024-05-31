#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 22             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal         # Queue (partition) name
#SBATCH -t 08:00:00       # Run time (hh:mm:ss)
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
# Prepare merge list file
#=======================================================================#

# Create or overwrite the merge list file
merge_list_file="${TEMP}/merge_list.txt"
> $merge_list_file

# Generate launcher task file for smaller chromosomes
task_file="${TEMP}/launcher_tasks.txt"
> $task_file

#=======================================================================#
# Process large chromosomes sequentially
#=======================================================================#

for CHR in {1..4}; do
    # Define input VCF file name
    VCF_FILE="${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz"

    # Define intermediate and output file names
    NORM_VCF="${TEMP}/chr${CHR}.normalized.vcf.gz"
    FILTERED_VCF="${TEMP}/chr${CHR}.filtered.vcf.gz"
    PLINK_PREFIX="${TEMP}/chr${CHR}_biallelic"

    # Check if PLINK files already exist
    if [ ! -f "${PLINK_PREFIX}.pgen" ] || [ ! -f "${PLINK_PREFIX}.pvar" ] || [ ! -f "${PLINK_PREFIX}.psam" ]; then
        # Run bcftools and plink2 commands sequentially
        bcftools norm -m -both -o ${NORM_VCF} -O z ${VCF_FILE}
        bcftools view -m2 -M2 -v snps ${NORM_VCF} -o ${FILTERED_VCF} -O z
        plink2 --vcf ${FILTERED_VCF} --make-pgen --out ${PLINK_PREFIX}
    fi

    # Add the PLINK file prefix to the merge list
    echo "${PLINK_PREFIX}" >> $merge_list_file
done

#=======================================================================#
# Prepare tasks for smaller chromosomes using launcher
#=======================================================================#

for CHR in {5..22}; do
    # Define input VCF file name
    VCF_FILE="${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz"

    # Define intermediate and output file names
    NORM_VCF="${TEMP}/chr${CHR}.normalized.vcf.gz"
    FILTERED_VCF="${TEMP}/chr${CHR}.filtered.vcf.gz"
    PLINK_PREFIX="${TEMP}/chr${CHR}_biallelic"

    # Check if PLINK files already exist
    if [ ! -f "${PLINK_PREFIX}.pgen" ] || [ ! -f "${PLINK_PREFIX}.pvar" ] || [ ! -f "${PLINK_PREFIX}.psam" ]; then
        # Add bcftools and plink2 commands to the task file
        echo "bcftools norm -m -both -o ${NORM_VCF} -O z ${VCF_FILE} && bcftools view -m2 -M2 -v snps ${NORM_VCF} -o ${FILTERED_VCF} -O z && plink2 --vcf ${FILTERED_VCF} --make-pgen --out ${PLINK_PREFIX}" >> $task_file
    fi

    # Add the PLINK file prefix to the merge list
    echo "${PLINK_PREFIX}" >> $merge_list_file
done

#=======================================================================#
# Run the launcher for smaller chromosomes if task file is not empty
#=======================================================================#

if [ -s $task_file ]; then
    module load launcher
    export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
    export LAUNCHER_RMI=SLURM
    export LAUNCHER_JOB_FILE=$task_file

    $LAUNCHER_DIR/paramrun
fi

#=======================================================================#
# Wait for all tasks to complete
#=======================================================================#

wait

#=======================================================================#
# Merge the PLINK2 files using the merge list
#=======================================================================#

plink2 --pmerge-list $merge_list_file --make-pgen --out ${TEMP}/ah_merged

