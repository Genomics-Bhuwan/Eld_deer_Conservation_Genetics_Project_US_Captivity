#### Genome annotation for Burmese-brow antelered deer(Rucervus eldii thamin).
- We used Red deer reference genome assembly and annotation file of Red deer for annotation using LiftOff pipeline.
- Link: https://github.com/agshumate/Liftoff

#### Step 1. Slurm script for annotation
```bash
#!/bin/bash
#SBATCH --job-name=liftoff_reddeer
#SBATCH --output=logs/liftoff_%j.out
#SBATCH --error=logs/liftoff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ==============================================================================
# Pipeline: Structural Annotation Transfer via Liftoff
# Target Species: Rucervus eldii (Eld's deer)
# Reference Species: Cervus elaphus (Red deer)
# Tool: Liftoff v1.6.3 (Apptainer Container)
# ==============================================================================

set -euo pipefail

# 1. Environment and Executive Paths
APPTAINER_BIN="/usr/bin/apptainer"
WORKDIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis/Genome_Annotation/liftoff"
GENOME_DIR="/home/bistbs/Elds_deer_Population_Genetic_Analysis/Genome_Annotation"

# 2. Input and Output File Definitions
TARGET_FASTA="${GENOME_DIR}/GCA_054824845.1_mRucEld1.hap2_genomic.fna"
REF_FASTA="${GENOME_DIR}/GCF_910594005.1_mCerEla1.1_genomic.fna"
REF_GFF="${GENOME_DIR}/GCF_910594005.1_mCerEla1.1_genomic.gff"
CONTAINER="${WORKDIR}/liftoff_v1.6.3.sif"

OUT_GFF="${WORKDIR}/Rucervus_eldii_mCerEla1.1_liftoff.gff3"
UNMAPPED_TXT="${WORKDIR}/unmapped_features_reddeer.txt"

# 3. Setup Directory Structure
mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}" || exit 1

# 4. Logging Environment Info
echo "========================================================================"
echo "Job Started: $(date)"
echo "Host Node: $(hostname)"
echo "Allocated Threads: ${SLURM_CPUS_PER_TASK}"
echo "Target Assembly: ${TARGET_FASTA}"
echo "Reference Assembly: ${REF_FASTA}"
echo "Reference Annotation: ${REF_GFF}"
echo "========================================================================"

# 5. Execute Liftoff
${APPTAINER_BIN} exec "${CONTAINER}" liftoff \
  -g "${REF_GFF}" \
  -o "${OUT_GFF}" \
  -u "${UNMAPPED_TXT}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  -flank 0.1 \
  -cds \
  -copies \
  "${TARGET_FASTA}" \
  "${REF_FASTA}"

echo "========================================================================"
echo "Job Finished: $(date)"
echo "Output Annotation: ${OUT_GFF}"
echo "Unmapped Log: ${UNMAPPED_TXT}"
echo "========================================================================"
```
