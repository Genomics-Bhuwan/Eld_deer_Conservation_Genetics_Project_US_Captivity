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
cd /anvil/scratch/x-bbist/Eld_Deer/Variant_Calling/DeepVariant
module load bcftools

# 1. Compress and Index gVCF profiles
for f in gvcfs/*.dv.g.vcf; do
    bcftools view -Oz -o "${f}.gz" "$f" &
done
wait

for f in gvcfs/*.dv.g.vcf.gz; do
    bcftools index "$f" &
done
wait

# 2. Jointly merge gVCF structures into a unified matrix
bcftools merge -g gvcfs/*.dv.g.vcf.gz -Oz -o vcf_results/Eld_Deer_35_samples_merged.g.vcf.gz

# 3. Extract the comprehensive Multi-Sample variant layer
bcftools convert -O z -o vcf_results/Eld_Deer_35_samples_raw.vcf.gz vcf_results/Eld_Deer_35_samples_merged.g.vcf.gz
bcftools index vcf_results/Eld_Deer_35_samples_raw.vcf.gz
```

Step 4: Downstream Filtering & Evaluation
Run your downstream splits to slice out your variant subgroups:

```bash
mkdir -p vcf_results/With_Indels vcf_results/Without_Indels
RAW_VCF=vcf_results/Eld_Deer_35_samples_raw.vcf.gz

# A. Structural/VEP Workflow (High Quality, Includes Indels)
vcftools --gzvcf $RAW_VCF \
         --minQ 30 \
         --recode --recode-INFO-all \
         --out vcf_results/With_Indels/Eld_Deer_merged_dv_with_indels

# B. Population Genomics Workflow (High Quality, Strict Biallelic SNPs)
vcftools --gzvcf $RAW_VCF \
         --minQ 30 \
         --min-alleles 2 \
         --max-alleles 2 \
         --remove-indels \
         --recode --recode-INFO-all \
         --out vcf_results/Without_Indels/Eld_Deer_merged_dv_biallelic_snps
````

Step 5: Quantify Outputs
To output your final metrics and see how your pipeline stacks up against GATK:

```bash
echo "Total Bi-allelic SNPs generated via DeepVariant:"
bcftools view -H vcf_results/Without_Indels/Eld_Deer_merged_dv_biallelic_snps.recode.vcf |
```
