#### Variant Annotation for Eld's deer using SnpEff
```bash
#!/bin/bash
# ==============================================================================
# Full Population Genomics Pipeline: SnpEff Annotation & Consequence Extraction
# Species: Rucervus eldii (Eld's deer)
# Samples: 35
# Primary Target: Autosomes CM139461.1 - CM139488.1 (Biallelic)
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Path & Environment Setup
# ------------------------------------------------------------------------------
WORKDIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis/SNPeff_Elds_deer"
JAVA_EXEC="${WORKDIR}/jdk-21.0.2+13/bin/java"
SNPEFF_JAR="${WORKDIR}/snpEff/snpEff.jar"
SNPSIFT_JAR="${WORKDIR}/snpEff/SnpSift.jar"

RAW_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_with_indels.recode.vcf"
AUTOSOME_LIST="${WORKDIR}/chr28_list.txt"
CHR28_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_chr1-28_with_indels.vcf"
BIALLELIC_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_chr1-28_biallelic_with_indels.vcf"
ANNOTATED_VCF="${WORKDIR}/Eld_Deer_35_samples_GQ20_chr1-28_biallelic_annotated.vcf"

OUT_DIR="${WORKDIR}/Variants_by_Consequences"

mkdir -p "${OUT_DIR}"
cd "${WORKDIR}"

echo "=== STAGE 1: Creating Primary Autosome List (CM139461.1 - CM139488.1) ==="
cat << 'EOF' > "${AUTOSOME_LIST}"
CM139461.1
CM139462.1
CM139463.1
CM139464.1
CM139465.1
CM139466.1
CM139467.1
CM139468.1
CM139469.1
CM139470.1
CM139471.1
CM139472.1
CM139473.1
CM139474.1
CM139475.1
CM139476.1
CM139477.1
CM139478.1
CM139479.1
CM139480.1
CM139481.1
CM139482.1
CM139483.1
CM139484.1
CM139485.1
CM139486.1
CM139487.1
CM139488.1
EOF

echo "=== STAGE 2: Filtering VCF to 28 Primary Autosomes (Removing Scaffolds) ==="
awk 'NR==FNR {chr[$1]; next} /^#/ || ($1 in chr)' "${AUTOSOME_LIST}" "${RAW_VCF}" > "${CHR28_VCF}"

echo "=== STAGE 3: Restricting to Biallelic Sites Only ==="
bcftools view -m2 -M2 "${CHR28_VCF}" -O v -o "${BIALLELIC_VCF}"

echo "=== STAGE 4: Annotating Biallelic VCF with SnpEff ==="
${JAVA_EXEC} -Xmx24g -jar "${SNPEFF_JAR}" \
  -v Rucervus_eldii \
  -stats "${WORKDIR}/snpeff_summary_chr1-28_biallelic.html" \
  "${BIALLELIC_VCF}" > "${ANNOTATED_VCF}"

echo "=== STAGE 5: Partitioning VCF into Functional Impact Classes ==="

# 5a. Loss of Function (LoF / HIGH Impact)
echo "Extracting LoF (HIGH Impact)..."
${JAVA_EXEC} -jar "${SNPSIFT_JAR}" filter "( ANN[*].IMPACT = 'HIGH' )" \
  "${ANNOTATED_VCF}" > "${OUT_DIR}/Eld_Deer_chr1-28_biallelic_LoF.vcf"

# 5b. Missense (MODERATE Impact)
echo "Extracting Missense (MODERATE Impact)..."
${JAVA_EXEC} -jar "${SNPSIFT_JAR}" filter "( ANN[*].EFFECT HAS 'missense_variant' )" \
  "${ANNOTATED_VCF}" > "${OUT_DIR}/Eld_Deer_chr1-28_biallelic_missense.vcf"

# 5c. Synonymous (LOW Impact)
echo "Extracting Synonymous (LOW Impact)..."
${JAVA_EXEC} -jar "${SNPSIFT_JAR}" filter "( ANN[*].EFFECT HAS 'synonymous_variant' )" \
  "${ANNOTATED_VCF}" > "${OUT_DIR}/Eld_Deer_chr1-28_biallelic_synonymous.vcf"

# 5d. Intergenic (MODIFIER Impact)
echo "Extracting Intergenic..."
${JAVA_EXEC} -jar "${SNPSIFT_JAR}" filter "( ANN[*].EFFECT HAS 'intergenic_region' )" \
  "${ANNOTATED_VCF}" > "${OUT_DIR}/Eld_Deer_chr1-28_biallelic_intergenic.vcf"

echo "=== STAGE 6: Exporting Genotype Matrices (.txt) and Missingness Filtering ==="
for TYPE in LoF missense synonymous intergenic; do
  VCF_FILE="${OUT_DIR}/Eld_Deer_chr1-28_biallelic_${TYPE}.vcf"
  TXT_FILE="${OUT_DIR}/Eld_Deer_chr1-28_biallelic_${TYPE}_genotypes.txt"
  
  MISS20_VCF="${OUT_DIR}/Eld_Deer_chr1-28_biallelic_${TYPE}_max20miss.vcf"
  MISS20_TXT="${OUT_DIR}/Eld_Deer_chr1-28_biallelic_${TYPE}_max20miss_genotypes.txt"
  
  echo "Processing ${TYPE} text matrix..."
  (bcftools query -l "${VCF_FILE}" | tr '\n' '\t' | sed 's/\t$/\n/' | awk '{print "CHROM\tPOS\t" $0}'; \
   bcftools query -f '%CHROM\t%POS[\t%GT]\n' "${VCF_FILE}") > "${TXT_FILE}"
   
  echo "Applying <=20% missingness filter to ${TYPE}..."
  bcftools view -i 'F_MISSING <= 0.2' "${VCF_FILE}" -O v -o "${MISS20_VCF}"
  
  (bcftools query -l "${MISS20_VCF}" | tr '\n' '\t' | sed 's/\t$/\n/' | awk '{print "CHROM\tPOS\t" $0}'; \
   bcftools query -f '%CHROM\t%POS[\t%GT]\n' "${MISS20_VCF}") > "${MISS20_TXT}"
done

echo "=============================================================================="
echo "PIPELINE COMPLETE!"
echo "Outputs generated in: ${OUT_DIR}"
echo "=============================================================================="
echo -n "Missense Sites:   " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_missense_sites.vcf | wc -l
echo -n "Synonymous Sites: " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_synonymous_sites.vcf | wc -l
echo -n "LoF/Splicing:     " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_lof_sites.vcf | wc -l
echo -n "Intergenic Sites: " && grep -v "^#" ${OUTPUT_DIR}/Eld_deer_intergenic_sites.vcf | wc -l
echo "---------------------------------------------------"
```
