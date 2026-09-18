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

#### Step 2. Consequences or impact class categorization based on synonymous variant, missense variant and loss of function.
```bash

#!/bin/bash
# ==============================================================================
# SnpSift Variant Extraction by Functional Consequence
# Target Species: Rucervus eldii (Eld's deer)
# Input VCF: Eld_Deer_35_samples_GQ20_with_indels_annotated.vcf
# ==============================================================================

set -euo pipefail

# 1. Define Paths and Executables
WORKING_DIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis/SNPeff_Elds_deer"
SNPSIFT_JAR="${WORKING_DIR}/snpEff/SnpSift.jar"
JAVA_EXEC="${WORKING_DIR}/jdk-21.0.2+13/bin/java"

INPUT_VCF="${WORKING_DIR}/Eld_Deer_35_samples_GQ20_with_indels_annotated.vcf"
OUTPUT_DIR="${WORKING_DIR}/Variants_by_Consequences"

# 2. Create the Output Directory
echo "Creating output directory: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# 3. Extract Missense Variants
echo "Extracting Missense variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'missense_variant'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Eld_deer_missense_sites.vcf

# 4. Extract Synonymous Variants
echo "Extracting Synonymous variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'synonymous_variant'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Eld_deer_synonymous_sites.vcf

# 5. Extract Loss of Function (LoF), Inframe Indels, & Splicing Variants
echo "Extracting Loss of Function (LoF) & Splicing variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "(ANN[*].EFFECT has 'transcript_ablation') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'stop_lost') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'inframe_insertion') | (ANN[*].EFFECT has 'inframe_deletion') | (ANN[*].EFFECT has 'splice_region_variant')" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Eld_deer_lof_sites.vcf

# 6. Extract Intergenic Variants
echo "Extracting Intergenic variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'intergenic_region'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Eld_deer_intergenic_sites.vcf

# -----------------------------------------------------------------
# Verification & Summary Block
# -----------------------------------------------------------------
echo "---------------------------------------------------"
echo "Filtering Complete! Total variant counts extracted:"
echo "---------------------------------------------------"
echo -n "Missense Sites:   " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_missense_sites.vcf | wc -l
echo -n "Synonymous Sites: " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_synonymous_sites.vcf | wc -l
echo -n "LoF/Splicing:     " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_lof_sites.vcf | wc -l
echo -n "Intergenic Sites: " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_intergenic_sites.vcf | wc -l
echo "---------------------------------------------------"
```
