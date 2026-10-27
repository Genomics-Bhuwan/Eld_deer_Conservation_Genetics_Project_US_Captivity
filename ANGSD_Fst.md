#### Step 1. 
```bash
#!/bin/bash
#SBATCH --job-name=1_SAF_parallel
#SBATCH --output=1_SAF_%A_%a.out
#SBATCH --error=1_SAF_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --partition=shared
#SBATCH --array=0-2                   # all three run in parallel; use 0-2%1 to serialize
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load biocontainers
module load angsd
module load samtools

THREADS=16
BASE="/anvil/projects/x-bio260092/Eld_der/ANGSD/Fst_across_subspecies"
REF="${BASE}/mRucEld1.hap1_genomic.fna"
AUTOSOMES_FILE="${BASE}/autosomes.txt"
OUT_DIR="${BASE}/angsd_out/FST_subspecies"

mkdir -p ${OUT_DIR}

POPS=("thamin" "siamensis" "hainanus")
MININDS=(5 5 15)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}
MININD=${MININDS[$SLURM_ARRAY_TASK_ID]}

echo "Task ${SLURM_ARRAY_TASK_ID}: ${POP}, minInd=${MININD}, start $(date)"

angsd -P ${THREADS} \
    -out ${OUT_DIR}/${POP} \
    -bam ${BASE}/${POP}_bamlist.txt \
    -rf ${AUTOSOMES_FILE} \
    -anc ${REF} -ref ${REF} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -minMapQ 30 -minQ 30 -C 50 -baq 1 \
    -minInd ${MININD} \
    -doCounts 1 -GL 1 -doSaf 1

echo "Completed ${POP} at $(date)"
```

#### Step 2. 1-SFS
```bash
#!/bin/bash
#SBATCH --job-name=2_SFS_1D_parallel
#SBATCH --output=2_SFS_1D_%A_%a.out
#SBATCH --error=2_SFS_1D_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G
#SBATCH --time=70:00:00
#SBATCH --partition=shared
#SBATCH --array=0-2
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=END,FAIL

# Load required modules
module load biocontainers
module load angsd

# Environment setup
THREADS=32
OUT_DIR="/anvil/projects/x-bio260092/Eld_der/ANGSD/Fst_across_subspecies/angsd_out/FST_subspecies"

# Map array index to subspecies name
POPS=("thamin" "siamensis" "hainanus")
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

echo "======================================================"
echo "Job Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing 1D SFS for Subspecies: ${POP}"
echo "Start Time: $(date)"
echo "======================================================"

realSFS ${OUT_DIR}/${POP}.saf.idx \
    -P ${THREADS} \
    -fold 1 > ${OUT_DIR}/${POP}.sfs

echo "======================================================"
echo "Completed 1D SFS for ${POP} at $(date)"
echo "======================================================"
```
