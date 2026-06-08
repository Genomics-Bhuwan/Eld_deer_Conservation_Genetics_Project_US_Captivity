#### Deepvariant
```bash
#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=batch
#SBATCH --job-name=DeepVariant_Array
#SBATCH --array=1-35
#SBATCH --output=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant/logs/dv_%A_%a.out
#SBATCH --error=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant/logs/dv_%A_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu

# ----------------- Load Anvil Modules ----------------- #
module purge
module load biocontainers
module load deepvariant

# ------------------- Directory Paths ------------------ #
BASE_DIR=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant
INPUT_DIR=${BASE_DIR}/input_bams
OUTPUT_DIR=${BASE_DIR}/gvcfs
REF_GENOME=${INPUT_DIR}/ref_genome.fna

# ---------------- Map Array ID to File ---------------- #
# Dynamically grabs the matching file row based on task index
BAM_FILE=$(ls ${INPUT_DIR}/*.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
BAM_BASENAME=$(basename "$BAM_FILE")

# Generate clean sample ID (strips both naming conventions)
SAMPLE_NAME=${BAM_BASENAME%_final.bam}
SAMPLE_NAME=${SAMPLE_NAME%_dedup.bam}

echo "Processing Sample Index ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME}"
echo "Using BAM: ${BAM_FILE}"

# ------------------- Run DeepVariant ------------------ #
# BioContainers natively wrappers the binary path. We just call it directly!
run_deepvariant \
    --model_type=WGS \
    --ref=${REF_GENOME} \
    --reads=${BAM_FILE} \
    --output_vcf=${OUTPUT_DIR}/${SAMPLE_NAME}.dv.vcf \
    --output_gvcf=${OUTPUT_DIR}/${SAMPLE_NAME}.dv.g.vcf \
    --num_shards=${SLURM_NTASKS_PER_NODE} \
    --vcf_stats_report=true \
    --disable_small_model=true \
    --logging_dir=${OUTPUT_DIR}/logs

echo "Completed processing ${SAMPLE_NAME} at $(date)"
```

#### Step 3: Joint Consolidation & Population Genotyping
For population analysis across 35 individuals, we should perform joint-genotyping on the generated .g.vcf data rather than isolating single-sample VCF vectors. This protects your population stats from massive artificial missing data patterns.

```bash
#!/bin/bash
#SBATCH --job-name=EldDeer_GLnexus
#SBATCH --output=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant/logs/glnexus_%j.out
#SBATCH --error=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant/logs/glnexus_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128       # Scale up to utilize the full 128 cores of the Anvil node
#SBATCH --mem=140G                # Allocates the entire system memory footprint available
#SBATCH -p wholenode              # Explicitly targets Anvil's node-exclusive whole node partition
#SBATCH --time=54:00:00           # Cleaned time configuration boundary (54 hours)

# ---------------- Load Environments ----------------- #
module load biocontainers/default
module load glnexus/1.4.1
module load bcftools/1.17

# ---------------- Set Variables --------------------- #
WORKING_DIR=/anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant
OUTPUT_DIR=${WORKING_DIR}/cohort_vcf

cd ${WORKING_DIR}

echo "=================================================="
echo "Starting GLnexus joint consolidation at $(date)"
echo "=================================================="

# Safety reset: Clear any broken temporary databases from previous attempts
rm -rf ${WORKING_DIR}/GLnexus.DB

# Execute the core joint calling
# GLnexus uses all *.dv.g.vcf.gz files in the folder and creates a unified binary BCF file
glnexus_cli \
  --config DeepVariant \
  --threads ${SLURM_CPUS_PER_TASK} \
  --mem-gbytes 115 \
  *.dv.g.vcf.gz > ${OUTPUT_DIR}/Eld_Deer_35_samples_raw.bcf

echo "GLnexus finished. Starting bcftools conversion to compressed VCF..."

# ---------------- Format Conversion ------------------ #
# Convert the binary BCF file into a standard, compressed multi-sample VCF
bcftools view \
  -Oz \
  -o ${OUTPUT_DIR}/Eld_Deer_35_samples_joint.vcf.gz \
  ${OUTPUT_DIR}/Eld_Deer_35_samples_raw.bcf

# Index the multi-sample population VCF
bcftools index ${OUTPUT_DIR}/Eld_Deer_35_samples_joint.vcf.gz

# Clean up the massive raw BCF to save your scratch space
rm ${OUTPUT_DIR}/Eld_Deer_35_samples_raw.bcf

echo "=================================================="
echo "Pipeline completely completed at $(date)!"
echo "=================================================="
```

#### Step 4: Downstream Filtering & Evaluation
Run your downstream splits to slice out your variant subgroups:

```bash
cd /anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant/cohort_vcf/VCF_filtering_for_DeepVariant_Data

# Load vcftools if it isn't already active
module load biocontainers/default
module load vcftools/0.1.16

# Run the test filtration command
vcftools \
  --gzvcf ../Eld_Deer_35_samples_joint.vcf.gz \
  --remove-indels \
  --maf 0.05 \
  --minGQ 20 \
  --max-missing 0.20 \
  --recode \
  --out Eld_Deer_35_samples_GQ20_test
```

