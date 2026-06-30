```bash
#!/bin/bash
#SBATCH --job-name=EldsDeer_SnpEff
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00

set -euo pipefail

# --- Paths & Project Setup ---
PROJECT=$HOME/Elds_deer_Population_Genetic_Analysis/SNPeff_Elds_deer
SNPEFF=${PROJECT}/snpEff

# Custom genome identifier and name
GENOME_ID="mRucEld1"
GENOME_NAME="Elds_deer"

# Inputs based on your file listing
GENOME_FASTA="${PROJECT}/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
GENOME_GFF="${PROJECT}/GCF_910594005.1_mCerEla1.1_genomic.gff"
INPUT_VCF="${PROJECT}/Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf"

# Outputs
OUTPUT_VCF="${PROJECT}/Eld_Deer_35_samples_annotated.vcf.gz"
REPORT_HTML="${PROJECT}/snpeff_report.html"

# --- 1. Create Database Directory ---
echo "Creating database directory..."
mkdir -p ${SNPEFF}/data/${GENOME_ID}

# --- 2. Prepare Genome Files ---
echo "Preparing genome annotation files..."

# Note: Your GFF wasn't gzipped in the 'ls' output, 
# so we copy it directly instead of using gunzip.
cp ${GENOME_GFF} ${SNPEFF}/data/${GENOME_ID}/genes.gff
cp ${GENOME_FASTA} ${SNPEFF}/data/${GENOME_ID}/sequences.fa

# --- 3. Update the SnpEff Config ---
echo "Checking snpEff.config..."

if ! grep -q "^${GENOME_ID}\.genome" ${SNPEFF}/snpEff.config; then
    echo "" >> ${SNPEFF}/snpEff.config
    echo "${GENOME_ID}.genome : ${GENOME_NAME}" >> ${SNPEFF}/snpEff.config
fi

# --- 4. Build Database ---
echo "Building SnpEff database..."
cd ${SNPEFF}

# Added absolute path to snpEff.jar
java -Xmx70g -jar ${SNPEFF}/snpEff.jar build \
    -gff3 \
    -v \
    ${GENOME_ID}

# --- 5. Annotate the Variants ---
echo "Annotating variants..."

# Pipe the output through bgzip since your input was unzipped VCF 
# but we want a compressed output VCF.
# Added absolute path to snpEff.jar
java -Xmx32g -jar ${SNPEFF}/snpEff.jar ann \
    -v \
    -stats ${REPORT_HTML} \
    ${GENOME_ID} \
    ${INPUT_VCF} \
    | bgzip > ${OUTPUT_VCF}

# --- 6. Index the Output VCF ---
echo "Indexing annotated VCF..."
tabix -p vcf ${OUTPUT_VCF}

# --- 7. Check the Output ---
echo "Annotation complete."
echo ""
echo "Generated files:"
echo "  ${OUTPUT_VCF}"
echo "  ${OUTPUT_VCF}.tbi"
echo "  ${REPORT_HTML}"
echo ""

echo "Quick verification:"
bcftools view ${OUTPUT_VCF} | head -20

```
