#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1		  # Total # of nodes (must be 1 for serial)
#SBATCH -n 22		  # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal   	  # Queue (partition) name
#SBATCH -t 12:00:00	  # Run time (hh:mm:ss)
#SBATCH -A OTH21060	  # Project/Allocation name (req'd if you have more than 1)
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

# Create or overwrite the merge list file
merge_list_file="${TEMP}/merge_list.txt"
> $merge_list_file

# Generate launcher task file
task_file="${TEMP}/launcher_tasks.txt"
> $task_file

# Loop through each chromosome and generate tasks
for CHR in "${CHROMOSOMES[@]}"; do
    # Define input VCF file name
    VCF_FILE="${INPUT}/chr${CHR}.dbGaP.dose.vcf.gz"

    # Define intermediate and output file names
    NORM_VCF="${TEMP}/chr${CHR}.normalized.vcf.gz"
    FILTERED_VCF="${TEMP}/chr${CHR}.filtered.vcf.gz"
    PLINK_PREFIX="${TEMP}/chr${CHR}_biallelic"

    # Initialize an empty command string
    command=""

    # Check if normalization is needed
    if [ ! -f "${NORM_VCF}" ]; then
        command="bcftools norm -m -both -o ${NORM_VCF} -O z ${VCF_FILE}"
    fi

    # Check if filtering is needed
    if [ ! -f "${FILTERED_VCF}" ]; then
        if [ -n "$command" ]; then
            command="${command} && "
        fi
        command="${command}bcftools view -m2 -M2 -v snps ${NORM_VCF} -o ${FILTERED_VCF} -O z"
    fi

    # Check if PLINK conversion is needed
    if [ ! -f "${PLINK_PREFIX}.pgen" ] || [ ! -f "${PLINK_PREFIX}.pvar" ] || [ ! -f "${PLINK_PREFIX}.psam" ]; then
        if [ -n "$command" ]; then
            command="${command} && "
        fi
        command="${command}plink2 --vcf ${FILTERED_VCF} --make-pgen --out ${PLINK_PREFIX}"
    fi

    # Add the command to the task file if there are tasks to be done
    if [ -n "$command" ]; then
        echo "$command" >> $task_file
    fi

    # Add the PLINK file prefix to the merge list
    echo "${PLINK_PREFIX}" >> $merge_list_file
done

# Submit tasks using the launcher
module load launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$task_file

$LAUNCHER_DIR/paramrun

# Merge the PLINK2 files using the merge list
plink2 --pmerge-list $merge_list_file --make-pgen --out ${TEMP}/ah_merged



#plink 	--bfile
#	--score all_neu.PRS_input.snpRes 12 5 8 header sum center
#	--out regpc_neu.PGI
