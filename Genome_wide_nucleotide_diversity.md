```bash
#!/bin/bash -l
#SBATCH --job-name=Elds_Diversity
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=120G
#SBATCH --partition=shared
#SBATCH --output=/anvil/scratch/x-bbist/ANGSD_new/logs/diversity_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

module load anaconda/2024.02-py311

ANGSD="$HOME/.conda/envs/angsd_env/bin/angsd"
REALSFS="$HOME/.conda/envs/angsd_env/bin/realSFS"
THETASTAT="$HOME/.conda/envs/angsd_env/bin/thetaStat"

WORKING_DIR="/anvil/scratch/x-bbist/ANGSD_new"
REFERENCE="${WORKING_DIR}/mRucEld1.hap1_genomic.fna"
BAM_LIST="${WORKING_DIR}/bam_list_35.txt"
OUT_DIR="${WORKING_DIR}/Genome_wide_nucleotide_diversity"

mkdir -p "$OUT_DIR"

# 1. Generate Joint SAF for all 35 samples
echo "[$(date)] Generating joint SAF for all 35 samples..."
$ANGSD -P 36 \
    -bam "$BAM_LIST" \
    -ref "$REFERENCE" \
    -anc "$REFERENCE" \
    -dosaf 1 \
    -gl 1 \
    -C 50 \
    -minQ 20 \
    -minMapQ 30 \
    -baq 1 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -out "${OUT_DIR}/all35_pop"

# 2. Estimate Population SFS
echo "[$(date)] Estimating population SFS..."
$REALSFS "${OUT_DIR}/all35_pop.saf.idx" -P 36 -fold 1 > "${OUT_DIR}/all35_pop.sfs"

# 3. Estimate Per-Site Thetas
echo "[$(date)] Calculating per-site theta values..."
$ANGSD -P 36 \
    -bam "$BAM_LIST" \
    -ref "$REFERENCE" \
    -anc "$REFERENCE" \
    -gl 1 \
    -doThetas 1 \
    -pest "${OUT_DIR}/all35_pop.sfs" \
    -dosaf 1 \
    -out "${OUT_DIR}/all35_thetas"

# 4. Calculate Global Summary Statistics (Pi, Watterson's Theta, Tajima's D)
echo "[$(date)] Calculating global summary statistics..."
$THETASTAT do_stat "${OUT_DIR}/all35_thetas.thetas.idx" -outnames "${OUT_DIR}/all35_diversity_global.stat"

# 5. Calculate Sliding Window Diversity (50kb windows, 10kb steps)
$THETASTAT do_stat "${OUT_DIR}/all35_thetas.thetas.idx" -win 50000 -step 10000 -outnames "${OUT_DIR}/all35_diversity_50kb.win.stat"

echo "[$(date)] Nucleotide diversity estimation complete."

```
