#### Running PSMC-FAC for high coverage and low coverage sampling for Eld's deer three subspecies

```
#!/bin/bash -l
#SBATCH --account=bio260092
#SBATCH --job-name=psmc_mpileup
#SBATCH --partition=shared
#SBATCH --array=1-35%10
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=42:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/anvil/projects/x-bio260092/Eld_der/PSMC_FAC/logs/mpileup_%A_%a.out
#SBATCH --error=/anvil/projects/x-bio260092/Eld_der/PSMC_FAC/logs/mpileup_%A_%a.err

module load biocontainers
module load bcftools/1.17

WORKDIR="/anvil/projects/x-bio260092/Eld_der/PSMC_FAC"
cd $WORKDIR

REF_GENOME="GCA_054824845.1_mRucEld1.hap2_genomic.fna"

mkdir -p vcf logs

# Extract the specific BAM filename for this array index task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

# Clean prefix string for outputs
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

echo "=== Processing Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME} ==="

bcftools mpileup --threads 4 -C50 -f "$REF_GENOME" -Ou "$BAM_FILE" | \
bcftools call --threads 4 -c -Oz -o "vcf/${SAMPLE_NAME}.vcf.gz"

bcftools index -t "vcf/${SAMPLE_NAME}.vcf.gz"

echo "=== Finished Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME} ==="
```

#### Step 2. Convert vcf to fq file format
```
#!/bin/bash -l
#SBATCH --account=bio260092
#SBATCH --job-name=vcf2psmcfa
#SBATCH --partition=shared
#SBATCH --array=1-35
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/anvil/projects/x-bio260092/Eld_der/PSMC_FAC/logs/vcf2psmcfa_%A_%a.out
#SBATCH --error=/anvil/projects/x-bio260092/Eld_der/PSMC_FAC/logs/vcf2psmcfa_%A_%a.err

module load biocontainers
module load bcftools/1.17

WORKDIR="/anvil/projects/x-bio260092/Eld_der/PSMC_FAC"
cd $WORKDIR

PSMC_DIR="/anvil/projects/x-bio260092/Eld_der/PSMC_FAC/psmc"
FQ2PSMCFA="$PSMC_DIR/utils/fq2psmcfa"

mkdir -p fq psmcfa logs

# Get BAM filename and clean sample name
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

VCF_IN="vcf/${SAMPLE_NAME}.vcf.gz"
FQ_OUT="fq/${SAMPLE_NAME}.fq.gz"
PSMCFA_OUT="psmcfa/${SAMPLE_NAME}.psmcfa"

# Dictionary storing exact average coverage for CRR samples
declare -A CRR_CVG=(
    ["CRR773090_final"]="2.92884" ["CRR773091_final"]="2.71407" ["CRR773093_final"]="2.82422"
    ["CRR773094_final"]="3.07432" ["CRR773095_final"]="2.96207" ["CRR773096_final"]="2.89966"
    ["CRR773097_final"]="2.82668" ["CRR773098_final"]="2.34077" ["CRR773099_final"]="2.62670"
    ["CRR773100_final"]="2.35612" ["CRR773101_final"]="2.31001" ["CRR773102_final"]="6.72462"
    ["CRR773103_final"]="7.82503" ["CRR773104_final"]="4.27568" ["CRR773105_final"]="5.70785"
    ["CRR773106_final"]="2.84364" ["CRR773107_final"]="2.76451" ["CRR773108_final"]="3.10650"
    ["CRR773109_final"]="3.00712" ["CRR773110_final"]="2.74377" ["CRR773111_final"]="3.51709"
    ["CRR773112_final"]="3.51194" ["CRR773113_final"]="3.24620" ["CRR773114_final"]="1.94030"
    ["CRR773115_final"]="3.79449" ["CRR773116_final"]="3.33746" ["CRR773117_final"]="2.65192"
)

# Parse parameters dynamically based on sample type
if [[ "$SAMPLE_NAME" =~ RE.*_10x ]]; then
    # Downsampled RE samples: Coverage is 10x
    cvg=10
    d_param=5
    D_param=$(( cvg * 2 ))
else
    # Original CRR samples: Retrieve exact depth from dictionary
    cvg_float=${CRR_CVG["$SAMPLE_NAME"]}
    
    # Calculate D = 2 * AvgCvg (rounded to integer)
    D_param=$(python3 -c "print(int(round(2 * ${cvg_float})))")
    
    # Set lower bound d = 2 (or 1 for <2x depth)
    d_param=$(python3 -c "print(max(1, int(round(${cvg_float} - 1))))")
fi

echo "=== Processing Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME} ==="
echo "Coverage: ${cvg:-$cvg_float}X | Parameters: -d ${d_param} -D ${D_param}"

# 1. Convert VCF to consensus FASTQ (.fq.gz) using system vcfutils.pl
bcftools view "$VCF_IN" | vcfutils.pl vcf2fq -d "$d_param" -D "$D_param" | gzip > "$FQ_OUT"

# 2. Convert FASTQ to PSMC input format (.psmcfa)
$FQ2PSMCFA -q20 "$FQ_OUT" > "$PSMCFA_OUT"

echo "=== Finished Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_NAME} ==="
```
