
#!/bin/bash -l
#SBATCH --job-name=EldsDeer_Het_Combined
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=90
#SBATCH --mem=100G
#SBATCH --partition=shared        
#SBATCH --array=0-26
#SBATCH --output=/anvil/scratch/x-bbist/Eld_Deer/Zheng_et_al_2024/All_27_dedeuplicated_BAMS/27_samples_Genome_Wide_Heterozygosity/logs/elds_deer_het_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# ----------------------------
# 1. ENVIRONMENT SETUP
# ----------------------------
module load anaconda/2024.02-py311

# Direct binary paths
ANGSD="$HOME/.conda/envs/angsd_env/bin/angsd"
REALSFS="$HOME/.conda/envs/angsd_env/bin/realSFS"

# ----------------------------
# 2. DEFINE MASTER PATHS
# ----------------------------
REFERENCE="/anvil/scratch/x-bbist/Eld_Deer/Reference_genome_asssembly/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
BAM_DIR="/anvil/scratch/x-bbist/Eld_Deer/Zheng_et_al_2024/All_27_dedeuplicated_BAMS"
OUTPUT_DIR="/anvil/scratch/x-bbist/Eld_Deer/Zheng_et_al_2024/All_27_dedeuplicated_BAMS/27_samples_Genome_Wide_Heterozygosity"
MASTER_FILE="$OUTPUT_DIR/All_27_Samples_Master_Heterozygosity.txt"

# Create output and cluster logging directories
mkdir -p "$OUTPUT_DIR/logs"
cd "$OUTPUT_DIR"

# ----------------------------
# 3. DYNAMIC BAM SELECTION VIA ARRAY ID
# ----------------------------
BAM_FILES=($(ls "$BAM_DIR"/*.bam))
BAM_INPUT="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$BAM_INPUT" ]; then
    echo "ERROR: No BAM file found mapped to Slurm Array ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

SAMPLE_NAME=$(basename "$BAM_INPUT" | cut -d'_' -f1)
OUT_PREFIX="Elds_Deer_${SAMPLE_NAME}"

# ----------------------------
# 4. INITIALIZE MASTER FILE HEADER (Only Task 0 does this)
# ----------------------------
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo -e "Species\tSample\tHeterozygosity" > "$MASTER_FILE"
fi

# ----------------------------
# 5. STEP 1: GENERATE SAF FILE
# ----------------------------
echo "[$(date)] Starting Step 1: Generating SAF file for Sample: $SAMPLE_NAME..."

$ANGSD -P 24 \
     -i "$BAM_INPUT" \
     -anc "$REFERENCE" \
     -ref "$REFERENCE" \
     -out "$OUT_PREFIX" \
     -dosaf 1 \
     -gl 1 \
     -C 50 \
     -minQ 20 \
     -minmapq 30 \
     -baq 1 \
     -uniqueOnly 1 \
     -remove_bads 1

# ----------------------------
# 6. STEP 2: CALCULATE FOLDED SFS
# ----------------------------
echo "[$(date)] Starting Step 2: Running realSFS for Sample: $SAMPLE_NAME..."

if [ -f "${OUT_PREFIX}.saf.idx" ]; then
    $REALSFS "${OUT_PREFIX}.saf.idx" \
         -P 24 \
         -fold 1 \
         -maxIter 100 \
         > "${OUT_PREFIX}.est.ml"
else
    echo "ERROR: SAF index file (${OUT_PREFIX}.saf.idx) not found."
    exit 1
fi

# ----------------------------
# 7. CALCULATE AND SAVE DUAL OUTPUTS
# ----------------------------
echo "------------------------------------------------"
echo "[$(date)] Saving heterozygosity metrics..."

# Run the division calculation
HET=$(awk '{print $2/($1+$2)}' "${OUT_PREFIX}.est.ml")

# OUTPUT A: Save to individual file per sample
INDIVIDUAL_FILE="${OUT_PREFIX}_final_heterozygosity.txt"
echo -e "Species\tSample\tHeterozygosity" > "$INDIVIDUAL_FILE"
echo -e "Elds_Deer\t${SAMPLE_NAME}\t$HET" >> "$INDIVIDUAL_FILE"
echo "Saved individual file to: $INDIVIDUAL_FILE"

# OUTPUT B: Append to the unified Master file safely using file locking (flock)
# This prevents jobs that finish at the exact same split-second from breaking the file.
exec 9>>"$MASTER_FILE"
flock -x 9
echo -e "Elds_Deer\t${SAMPLE_NAME}\t$HET" >> "$MASTER_FILE"
exec 9>&-

echo "Saved data row to master list: $MASTER_FILE"
echo "Sample ${SAMPLE_NAME} Heterozygosity: $HET"
echo "------------------------------------------------"
echo "[$(date)] Pipeline completed for $SAMPLE_NAME."
