#!/bin/bash
#SBATCH --job-name=eldii_annotation_snpeff
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --partition=batch
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ==============================================================================
# PIPELINE: End-to-End Genome Structural Annotation Transfer and Variant Annotation
# Target Species: Rucervus eldii (Eld's deer) | Assembly: mRucEld1.hap2 (GCA_054824845.1)
# Reference Species: Cervus elaphus (Red deer) | Assembly: mCerEla1.1 (GCF_910594005.1)
# Workflow Steps:
#   1. Liftoff v1.6.3 (Apptainer) - Homology-based GFF3 Annotation Transfer
#   2. Local OpenJDK 21 - Runtime environment configuration
#   3. SnpEff - Custom Database Construction (Rucervus_eldii)
#   4. SnpEff - Variant Functional Effect Prediction on 35-sample VCF
# ==============================================================================

set -euo pipefail

# ==============================================================================
# 1. GLOBAL PATH DEFINITIONS AND DIRECTORY SETUP
# ==============================================================================
BASE_DIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis"
GENOME_DIR="${BASE_DIR}/Genome_Annotation"
LIFTOFF_DIR="${GENOME_DIR}/liftoff"
SNPEFF_WORKDIR="${BASE_DIR}/SNPeff_Elds_deer"
SNPEFF_DIR="${SNPEFF_WORKDIR}/snpEff"

# Tool Executables & Containers
APPTAINER_BIN="/usr/bin/apptainer"
LIFTOFF_SIF="${LIFTOFF_DIR}/liftoff_v1.6.3.sif"
JAVA_BIN="${SNPEFF_WORKDIR}/jdk-21.0.2+13/bin/java"

# Reference Genome Inputs (Red Deer)
REF_FASTA="${GENOME_DIR}/GCF_910594005.1_mCerEla1.1_genomic.fna"
REF_GFF="${GENOME_DIR}/GCF_910594005.1_mCerEla1.1_genomic.gff"

# Target Genome Inputs & Core Inputs (Eld's Deer)
TARGET_FASTA="${GENOME_DIR}/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
INPUT_VCF="${SNPEFF_WORKDIR}/Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf"

# Output Definitions
LIFTOFF_GFF="${SNPEFF_WORKDIR}/Rucervus_eldii_mCerEla1.1_liftoff.gff3"
UNMAPPED_LOG="${LIFTOFF_DIR}/unmapped_features_reddeer.txt"
ANNOTATED_VCF="${SNPEFF_WORKDIR}/Eld_Deer_35_samples_GQ20_annotated.vcf"
SUMMARY_HTML="${SNPEFF_WORKDIR}/snpeff_summary.html"
DB_NAME="Rucervus_eldii"

# Create Required Execution Directories
mkdir -p "${BASE_DIR}/logs"
mkdir -p "${LIFTOFF_DIR}"
mkdir -p "${SNPEFF_WORKDIR}/logs"
mkdir -p "${SNPEFF_DIR}/data/${DB_NAME}"

echo "========================================================================"
echo "PIPELINE STARTED: $(date)"
echo "Host Node: $(hostname)"
echo "Threads Allocated: ${SLURM_CPUS_PER_TASK}"
echo "Target Genome: ${TARGET_FASTA}"
echo "Reference Genome: ${REF_FASTA}"
echo "Reference GFF3: ${REF_GFF}"
echo "Input VCF: ${INPUT_VCF}"
echo "========================================================================"

# ==============================================================================
# STEP 1: ANNOTATION TRANSFER VIA LIFTOFF
# ==============================================================================
echo ""
echo "[STEP 1/4] Running Liftoff to map GFF3 annotation from Red Deer to Eld's Deer..."

cd "${LIFTOFF_DIR}"

${APPTAINER_BIN} exec "${LIFTOFF_SIF}" liftoff \
  -g "${REF_GFF}" \
  -o "${LIFTOFF_GFF}" \
  -u "${UNMAPPED_LOG}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  -flank 0.1 \
  -cds \
  -copies \
  "${TARGET_FASTA}" \
  "${REF_FASTA}"

echo "Liftoff complete."
echo "Target GFF3 output created: ${LIFTOFF_GFF}"
echo "Unmapped features saved to: ${UNMAPPED_LOG}"

# ==============================================================================
# STEP 2: VERIFY AND SET UP JAVA 21 ENVIRONMENT
# ==============================================================================
echo ""
echo "[STEP 2/4] Verifying Java 21 environment for SnpEff..."

cd "${SNPEFF_WORKDIR}"

# Check if JDK 21 exists; download via wget if not present
if [ ! -f "${JAVA_BIN}" ]; then
    echo "Java 21 binary not found. Downloading OpenJDK 21 tarball..."
    wget -q https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.2%2B13/OpenJDK21U-jdk_x64_linux_hotspot_21.0.2_13.tar.gz
    tar -xzf OpenJDK21U-jdk_x64_linux_hotspot_21.0.2_13.tar.gz
fi

echo "Java Executable: $(${JAVA_BIN} -version 2>&1 | head -n 1)"

# ==============================================================================
# STEP 3: CONSTRUCT CUSTOM SNPEFF DATABASE
# ==============================================================================
echo ""
echo "[STEP 3/4] Preparing data directory and building SnpEff database (${DB_NAME})..."

# Copy FASTA and Liftoff GFF3 files into standard SnpEff structure
cp "${TARGET_FASTA}" "${SNPEFF_DIR}/data/${DB_NAME}/sequences.fa"
cp "${LIFTOFF_GFF}" "${SNPEFF_DIR}/data/${DB_NAME}/genes.gff"

# Append genome entry to snpEff.config if not already defined
if ! grep -q "${DB_NAME}.genome" "${SNPEFF_DIR}/snpEff.config"; then
    echo "Registering ${DB_NAME} in snpEff.config..."
    echo -e "\n# Rucervus eldii custom genome annotation\n${DB_NAME}.genome : Rucervus eldii" >> "${SNPEFF_DIR}/snpEff.config"
fi

# Execute SnpEff Database Build
cd "${SNPEFF_DIR}"

"${JAVA_BIN}" -Xmx16g -jar snpEff.jar build \
  -gff3 \
  -v \
  -noCheckProtein \
  -noCheckCds \
  "${DB_NAME}"

echo "SnpEff database build complete for ${DB_NAME}."

# ==============================================================================
# STEP 4: VARIANT EFFECT ANNOTATION ON 35-SAMPLE VCF
# ==============================================================================
echo ""
echo "[STEP 4/4] Annotating 35-sample VCF variants using ${DB_NAME} database..."

cd "${SNPEFF_WORKDIR}"

"${JAVA_BIN}" -Xmx24g -jar "${SNPEFF_DIR}/snpEff.jar" \
  -v "${DB_NAME}" \
  -stats "${SUMMARY_HTML}" \
  "${INPUT_VCF}" \
  > "${ANNOTATED_VCF}"

echo ""
echo "========================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY: $(date)"
echo "Transformed GFF3: ${LIFTOFF_GFF}"
echo "Annotated VCF: ${ANNOTATED_VCF}"
echo "SnpEff HTML Summary Report: ${SUMMARY_HTML}"
echo "========================================================================"
