#### Variant Annotation for Eld's deer using SnpEff
```bash
#!/bin/bash
# SnpEff Database Build & Variant Annotation (SNPs + Indels)
# Species: Rucervus eldii
# ==============================================================================

set -euo pipefail

# 1. Define Directories and Paths
WORKDIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis/SNPeff_Elds_deer"
SNPEFF_DIR="${WORKDIR}/snpEff"
JAVA_BIN="${WORKDIR}/jdk-21.0.2+13/bin/java"
DB_NAME="Rucervus_eldii"

# Input Files
TARGET_FASTA="/home/bistbs/Elds_deer_Population_Genetic_Analysis/Genome_Annotation/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
LIFTOFF_GFF="${WORKDIR}/Rucervus_eldii_mCerEla1.1_liftoff.gff3"
INPUT_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_with_indels.recode.vcf"

# Output Files
OUTPUT_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_with_indels_annotated.vcf"
SUMMARY_HTML="${WORKDIR}/snpeff_summary_with_indels.html"

# 2. Prepare Database Folder Structure
mkdir -p "${SNPEFF_DIR}/data/${DB_NAME}"

cp "${TARGET_FASTA}" "${SNPEFF_DIR}/data/${DB_NAME}/sequences.fa"
cp "${LIFTOFF_GFF}" "${SNPEFF_DIR}/data/${DB_NAME}/genes.gff"

# Register genome in config if missing
if ! grep -q "${DB_NAME}.genome" "${SNPEFF_DIR}/snpEff.config"; then
    echo -e "\n# Rucervus eldii custom genome annotation\n${DB_NAME}.genome : Rucervus eldii" >> "${SNPEFF_DIR}/snpEff.config"
fi

# 3. Build SnpEff Database
cd "${SNPEFF_DIR}"

"${JAVA_BIN}" -Xmx16g -jar snpEff.jar build \
  -gff3 \
  -v \
  -noCheckProtein \
  -noCheckCds \
  "${DB_NAME}"

# 4. Annotate VCF (SNPs + Indels)
cd "${WORKDIR}"

"${JAVA_BIN}" -Xmx24g -jar "${SNPEFF_DIR}/snpEff.jar" \
  -v "${DB_NAME}" \
  -stats "${SUMMARY_HTML}" \
  "${INPUT_VCF}" \
  > "${OUTPUT_VCF}"
```
