#### Variant Calling using DeepVariant  
- https://github.com/google/deepvariant

- DeepVariant is a deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file (Poplin et al., 2018).

- My rmdup_.bam file holds five samples in it but DeepVariant runs for one bam for each sample at a time. Therefore, I am firstly seperating all the five samples and making 5 .bam files. Then will be running DeepVariant.

```bash
#!/bin/bash

# ----------------Modules / Paths----------------- #
APPTAINER=/usr/bin/apptainer
DEEPVARIANT_IMAGE=docker://google/deepvariant:1.9.0

# ----------------User Variables------------------ #
INPUT_DIR=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates
OUTPUT_DIR=/localscratch/bistbs/DeepVariant
REF_GENOME=Dama_gazelle_hifiasm-ULONT_primary.fasta
BAM_FILE=all_samples_merged_rmdup.bam

# ----------------Create output directories------------------ #
mkdir -p ${OUTPUT_DIR}/logs

# ----------------Run DeepVariant---------------- #
echo "Starting DeepVariant at $(date)"

$APPTAINER run \
  -B ${INPUT_DIR}:/input \
  -B ${OUTPUT_DIR}:/output \
  ${DEEPVARIANT_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/${REF_GENOME} \
    --reads=/input/${BAM_FILE} \
    --output_vcf=/output/${BAM_FILE%.bam}.dv.vcf \
    --output_gvcf=/output/${BAM_FILE%.bam}.dv.g.vcf \
    --num_shards=$(nproc) \
    --vcf_stats_report=true \
    --disable_small_model=true \
    --logging_dir=/output/logs

echo "DeepVariant finished at $(date)"

```
---
#### Since, I have one merged deduplicated .bam contatining five samples. But deepvariant is asking for individual .bam for each sample. Therefore, I am splitting the bam samplewise with below code.
---
##### Split the files
```bash
cd /scratch/bistbs/DeepVariant

samtools split -f "%!.bam" all_samples_merged_rmdup.bam
```
#### Step 1. Extracting the BAMs.
```bash
# Load samtools
module load samtools

# Move to your working directory
cd /scratch/bistbs/DeepVariant

# Create output folder
mkdir -p per_sample_BAMs logs

# Extract BAMs per sample in parallel (16 threads each)
samtools view -@ 16 -b -r SRR17129394 all_samples_merged_rmdup.bam > per_sample_BAMs/SRR17129394.bam &
samtools view -@ 16 -b -r SRR17134085 all_samples_merged_rmdup.bam > per_sample_BAMs/SRR17134085.bam &
samtools view -@ 16 -b -r SRR17134086 all_samples_merged_rmdup.bam > per_sample_BAMs/SRR17134086.bam &
samtools view -@ 16 -b -r SRR17134087 all_samples_merged_rmdup.bam > per_sample_BAMs/SRR17134087.bam &
samtools view -@ 16 -b -r SRR17134088 all_samples_merged_rmdup.bam > per_sample_BAMs/SRR17134088.bam &

# Wait for all background jobs to finish
wait

```

#### Step 2. Sort and Index them.
```bash
# Index all BAMs (recommended before DeepVariant)
for bam in per_sample_BAMs/*.bam; do
    samtools index -@ 16 "$bam" &
done
wait

echo "✅ All sample BAMs extracted and indexed successfully!"
```

#### Step 3. Run DeepVariant for each .bam or each sample
```bash
#!/bin/bash -l
#SBATCH --time=250:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=DeepVariant

APPTAINER=/usr/bin/apptainer
DEEPVARIANT_IMAGE=docker://google/deepvariant:1.9.0

# ----------------User Variables------------------ #
INPUT_DIR=/scratch/bistbs/DeepVariant
OUTPUT_DIR=/scratch/bistbs/DeepVariant/DeepVariant_Results
REF_GENOME=Dama_gazelle_primary.fasta
BAM_FILE=all_samples_merged_rmdup.bam

# ----------------Create output directories------------------ #
mkdir -p ${OUTPUT_DIR}/logs

# ----------------Run DeepVariant---------------- #
echo "Starting DeepVariant at $(date)"

$APPTAINER run \
  -B ${INPUT_DIR}:/input \
  -B ${OUTPUT_DIR}:/output \
  ${DEEPVARIANT_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/${REF_GENOME} \
    --reads=/input/${BAM_FILE} \
    --output_vcf=/output/${BAM_FILE%.bam}.dv.vcf \
    --output_gvcf=/output/${BAM_FILE%.bam}.dv.g.vcf \
    --num_shards=$(nproc) \
    --vcf_stats_report=true \
    --disable_small_model=true \
    --logging_dir=/output/logs

echo "DeepVariant finished at $(date)"
```

#### Step 4. After running five such slurm script, we get the two output files for each individual sample. For example: 
a. /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/SRR17129394.bam.dv.g.vcf
b. /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/SRR17129394.bam.dv.vcf
- For our purpose of population genomic analysis we use files with .bam.dv.vcf.
- I am compressing all the .vcf file to adjust it for merging and make it a one giant vcf file.

```bash
# 1️. Compress each DeepVariant VCF using bcftools 
for f in SRR*.bam.dv.vcf; do
  bcftools view -Oz -o "${f}.gz" "$f"
done

# 2️. Index the compressed VCFs
for f in SRR*.bam.dv.vcf.gz; do
  bcftools index "$f"
done

# 3. Merge all samples into one combined VCF
bcftools merge -Oz -o merged_deepvariant.vcf.gz SRR*.bam.dv.vcf.gz

# 4️. Index the merged file
bcftools index merged_deepvariant.vcf.gz
```

#### Step 5. Filter the high-quality variants
- Keep only variants with quality ≥30
 ---
 a. Filtering but keeping the indels for Structural annotation analysis as well as VEP and SnpEff analysis.
 ---
```bash
# -------------------------------
# Filter with VCFtools: Keep all variants (SNPs + INDELs) with high quality
# -------------------------------
vcftools --gzvcf $RAW_VCF \
         --minQ 30 \
         --recode --recode-INFO-all \
         --out /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/With_Indels/Dama_gazelle_merged_dv_with_indels

# -------------------------------
# Check missingness per individual
# -------------------------------
vcftools --vcf /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/With_Indels/Dama_gazelle_merged_dv_with_indels.recode.vcf \
         --missing-indv \
         --out /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/With_Indels/Dama_gazelle_merged_dv_with_indels
```

 ---
 b. Filtering but removing the indels for population genomic analysis.
 ---
```bash
# -------------------------------
# Filter with VCFtools: Keep only SNPs with high quality (remove INDELs)
# -------------------------------
vcftools --gzvcf $RAW_VCF \
         --minQ 30 \
         --remove-indels \
         --recode --recode-INFO-all \
         --out /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/Without_Indels/Dama_gazelle_merged_dv_without_indels

# -------------------------------
# Check missingness per individual
# -------------------------------
vcftools --vcf /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/Without_Indels/Dama_gazelle_merged_dv_without_indels.recode.vcf \
         --missing-indv \
         --out /scratch/bistbs/DeepVariant/per_sample_BAMs/DeepVariant_Results/Without_Indels/Dama_gazelle_merged_dv_without_indels
```

---
Count the number of bi-allelic SNPs only. Compare the bi-allelic SNPs between both methods: GATK and DeepVariant.
---
References
---
- Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, Newburger D, Dijamco J, Nguyen N, Afshar PT, Gross SS. A universal SNP and small-indel variant caller using deep neural networks. Nature biotechnology. 2018 Nov;36(10):983-7
