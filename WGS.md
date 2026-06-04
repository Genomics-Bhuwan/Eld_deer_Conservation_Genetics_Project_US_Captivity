### Eld's deer pre-processing
```bash
#!/bin/bash
#SBATCH --job-name=EldDeer_pipeline
#SBATCH --output=logs/elddeer_%A_%a.out
#SBATCH --error=logs/elddeer_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=72:00:00
#SBATCH --partition=shared
#SBATCH --account=bio260092
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-26

# ── Modules ───────────────────────────────────────────────────────────────────
module purge
module load gcc/11.2.0
module load bwa/0.7.17
module load samtools/1.12
module load picard/2.25.7
module load cutadapt/2.10

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR="/anvil/scratch/x-bbist/Eld_Deer/Zheng_et_al_2024"
REF="/anvil/scratch/x-bbist/Eld_Deer/Reference_genome_asssembly/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
OUT_DIR="${PROJECT_DIR}/processed_results"
TRIMGALORE="/anvil/projects/x-bio260092/Guam_rail/Short_reads_Guam_rail_Project/software/TrimGalore-0.6.7/trim_galore"
JAVA=$(which java)

mkdir -p "${OUT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

# ── All 27 samples ────────────────────────────────────────────────────────────
SAMPLES=(
    CRR773090 CRR773091 CRR773093 CRR773094
    CRR773095 CRR773096 CRR773097 CRR773098
    CRR773099 CRR773100 CRR773101 CRR773102
    CRR773103 CRR773104 CRR773105 CRR773106
    CRR773107 CRR773108 CRR773109 CRR773110
    CRR773111 CRR773112 CRR773113 CRR773114
    CRR773115 CRR773116 CRR773117
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R1="${PROJECT_DIR}/${SAMPLE}/${SAMPLE}_f1.fq.gz"
R2="${PROJECT_DIR}/${SAMPLE}/${SAMPLE}_r2.fq.gz"

echo "=========================================="
echo "Array Job:  ${SLURM_ARRAY_JOB_ID}"
echo "Task Index: ${SLURM_ARRAY_TASK_ID}"
echo "Sample:     ${SAMPLE}"
echo "Node:       $(hostname)"
echo "Start:      $(date)"
echo "=========================================="

# ── 1. Adapter trimming ───────────────────────────────────────────────────────
echo ">>> Step 1: Adapter trimming for ${SAMPLE} ..."
TRIM_DIR="${OUT_DIR}/${SAMPLE}_trim"
mkdir -p "${TRIM_DIR}"

$TRIMGALORE \
    --paired \
    --cores $SLURM_CPUS_PER_TASK \
    --output_dir "${TRIM_DIR}" \
    "${R1}" "${R2}"

TRIM_R1="${TRIM_DIR}/${SAMPLE}_f1_val_1.fq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE}_r2_val_2.fq.gz"

echo ">>> Trimming done for ${SAMPLE}"

# ── 2. Index reference (only task 0 does this, others wait) ──────────────────
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    if [ ! -f "${REF}.bwt" ]; then
        echo ">>> Step 2: Building BWA index ..."
        bwa index "${REF}"
        samtools faidx "${REF}"
        echo ">>> Indexing done"
    else
        echo ">>> Index already exists, skipping"
    fi
else
    # Wait for index to be ready
    echo ">>> Waiting for BWA index to be ready ..."
    while [ ! -f "${REF}.bwt" ]; do
        sleep 60
    done
    echo ">>> Index ready, proceeding"
fi

# ── 3. Alignment ──────────────────────────────────────────────────────────────
echo ">>> Step 3: Aligning ${SAMPLE} ..."
bwa mem -t $SLURM_CPUS_PER_TASK \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${TRIM_R1}" "${TRIM_R2}" | \
    samtools view -bhS - | \
    samtools sort -@ $SLURM_CPUS_PER_TASK \
    -o "${OUT_DIR}/${SAMPLE}_temp_sorted.bam"

echo ">>> Alignment done for ${SAMPLE}"

# ── 4. Mark duplicates ────────────────────────────────────────────────────────
echo ">>> Step 4: Marking duplicates for ${SAMPLE} ..."
$JAVA -Xmx50g -jar $PICARD MarkDuplicates \
    I="${OUT_DIR}/${SAMPLE}_temp_sorted.bam" \
    O="${OUT_DIR}/${SAMPLE}_final.bam" \
    M="${OUT_DIR}/${SAMPLE}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

# ── 5. Index final BAM ────────────────────────────────────────────────────────
echo ">>> Step 5: Indexing final BAM for ${SAMPLE} ..."
samtools index "${OUT_DIR}/${SAMPLE}_final.bam"

# ── 6. Cleanup ────────────────────────────────────────────────────────────────
echo ">>> Step 6: Cleaning up temp files ..."
rm -f "${OUT_DIR}/${SAMPLE}_temp_sorted.bam"

echo ""
echo "=========================================="
echo "Pipeline complete for ${SAMPLE}"
echo "Output: ${OUT_DIR}/${SAMPLE}_final.bam"
echo "End: $(date)"
echo "=========================================="
```
