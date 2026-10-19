```bash
#!/bin/bash -l
#SBATCH --job-name=EldsDeer_Het
#SBATCH --time=54:00:00
#SBATCH --cpus-per-task=45
#SBATCH --mem=100G
#SBATCH --partition=shared
#SBATCH --array=0-34
#SBATCH --output=/anvil/scratch/x-bbist/Eld_Deer/ANGSD_new/Genome_wide_hetereozygosity/logs/het_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# ------------------------------------------------------------------------------
# 1. ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
module load anaconda/2024.02-py311

ANGSD="$HOME/.conda/envs/angsd_env/bin/angsd"
REALSFS="$HOME/.conda/envs/angsd_env/bin/realSFS"

# ------------------------------------------------------------------------------
# 2. DEFINE MASTER PATHS
# ------------------------------------------------------------------------------
WORKING_DIR="/anvil/scratch/x-bbist/ANGSD_new"
REFERENCE="${WORKING_DIR}/mRucEld1.hap1_genomic.fna"
OUTPUT_DIR="${WORKING_DIR}/Genome_wide_hetereozygosity"
MASTER_FILE="${OUTPUT_DIR}/All_35_Samples_Master_Heterozygosity.txt"

mkdir -p "$OUTPUT_DIR"
mkdir -p "${WORKING_DIR}/logs"

# ------------------------------------------------------------------------------
# 3. DYNAMIC BAM SELECTION VIA ARRAY ID
# ------------------------------------------------------------------------------
# Read sample BAM paths directly from the directory
BAM_FILES=($(ls ${WORKING_DIR}/*.bam | sort))
BAM_INPUT="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$BAM_INPUT" ]; then
    echo "ERROR: No BAM file found mapped to Slurm Array ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract Sample Name (e.g., CRR773090_final_dedup.bam -> CRR773090)
FILENAME=$(basename "$BAM_INPUT")
SAMPLE_NAME=$(echo "$FILENAME" | cut -d'_' -f1)
OUT_PREFIX="${OUTPUT_DIR}/Elds_Deer_${SAMPLE_NAME}"

# ------------------------------------------------------------------------------
# 4. INITIALIZE MASTER FILE HEADER (Only Task 0 does this)
# ------------------------------------------------------------------------------
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo -e "Species\tSample\tHeterozygosity" > "$MASTER_FILE"
fi

# ------------------------------------------------------------------------------
# 5. STEP 1: GENERATE SAF FILE
# ------------------------------------------------------------------------------
echo "[$(date)] Starting Step 1: Generating SAF file for Sample: $SAMPLE_NAME..."

$ANGSD -P 12 \
    -i "$BAM_INPUT" \
    -anc "$REFERENCE" \
    -ref "$REFERENCE" \
    -out "$OUT_PREFIX" \
    -dosaf 1 \
    -gl 1 \
    -C 50 \
    -minQ 20 \
    -minMapQ 30 \
    -baq 1 \
    -uniqueOnly 1 \
    -remove_bads 1

# ------------------------------------------------------------------------------
# 6. STEP 2: CALCULATE FOLDED SFS
# ------------------------------------------------------------------------------
echo "[$(date)] Starting Step 2: Running realSFS for Sample: $SAMPLE_NAME..."

if [ -f "${OUT_PREFIX}.saf.idx" ]; then
    $REALSFS "${OUT_PREFIX}.saf.idx" \
        -P 12 \
        -fold 1 \
        -maxIter 100 \
        > "${OUT_PREFIX}.est.ml"
else
    echo "ERROR: SAF index file (${OUT_PREFIX}.saf.idx) not found."
    exit 1
fi

# ------------------------------------------------------------------------------
# 7. CALCULATE AND SAVE DUAL OUTPUTS
# ------------------------------------------------------------------------------
echo "------------------------------------------------"
echo "[$(date)] Saving heterozygosity metrics..."

# Calculate Heterozygosity: Heterozygous sites / (Homozygous + Heterozygous sites)
HET=$(awk '{print $2/($1+$2)}' "${OUT_PREFIX}.est.ml")

# Save to individual sample file
INDIVIDUAL_FILE="${OUT_PREFIX}_final_heterozygosity.txt"
echo -e "Species\tSample\tHeterozygosity" > "$INDIVIDUAL_FILE"
echo -e "Elds_Deer\t${SAMPLE_NAME}\t$HET" >> "$INDIVIDUAL_FILE"

# Append to Unified Master File using flock to prevent race conditions
exec 9>>"$MASTER_FILE"
flock -x 9
echo -e "Elds_Deer\t${SAMPLE_NAME}\t$HET" >> "$MASTER_FILE"
exec 9>&-

echo "Saved data row to master file: $MASTER_FILE"
echo "Sample ${SAMPLE_NAME} Heterozygosity: $HET"
echo "------------------------------------------------"
echo "[$(date)] Pipeline completed for $SAMPLE_NAME."
```
