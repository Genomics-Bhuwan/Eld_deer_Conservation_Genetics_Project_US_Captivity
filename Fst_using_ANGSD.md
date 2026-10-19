```bash
#!/bin/bash
#SBATCH --job-name=Fst_ROH
#SBATCH --output=Fst_%j.out
#SBATCH --error=Fst_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=65
#SBATCH --mem=100G
#SBATCH --time=80:00:00
#SBATCH --partition=shared
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required software modules on Anvil
module load biocontainers
module load bioinfo
module load angsd
module load samtools

# -----------------------------------------------------------------------------
# DIRECTORY & REFERENCE SETUP
# -----------------------------------------------------------------------------
BASE_DIR="/anvil/projects/x-bio260092/Eld_der/ANGSD/Fst_across_subspecies"
OUT_DIR="${BASE_DIR}/FST"
mkdir -p "${OUT_DIR}/out"
cd "${BASE_DIR}"

# Reference sequence in your current directory
REF="${BASE_DIR}/mRucEld1.hap1_genomic.fna"

# Target chromosomes file
CHRS="${BASE_DIR}/autosomes.txt"

# -----------------------------------------------------------------------------
# CREATE SUBSPECIES BAM LISTS
# -----------------------------------------------------------------------------
# 1. Thamin (8 samples)
cat << 'EOF' > thamin_bamlist.txt
/anvil/projects/x-bio260092/Eld_der/RE113838_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE114379_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE114623_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE115125_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE115304_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE115445_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE115604_10X.bam
/anvil/projects/x-bio260092/Eld_der/RE116077_10X.bam
EOF

# 2. Siamensis (7 samples)
cat << 'EOF' > siamensis_bamlist.txt
/anvil/projects/x-bio260092/Eld_der/CRR773090_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773091_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773093_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773094_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773095_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773096_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773097_final_dedup.bam
EOF

# 3. Hainanus (20 samples)
cat << 'EOF' > hainanus_bamlist.txt
/anvil/projects/x-bio260092/Eld_der/CRR773098_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773099_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773100_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773101_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773102_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773103_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773104_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773105_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773106_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773107_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773108_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773109_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773110_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773111_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773112_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773113_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773114_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773115_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773116_final_dedup.bam
/anvil/projects/x-bio260092/Eld_der/CRR773117_final_dedup.bam
EOF

# -----------------------------------------------------------------------------
# STEP 1: CALCULATE SAF (Site Allele Frequency likelihoods) PER SUBSPECIES
# -----------------------------------------------------------------------------
echo "=== Step 1: Calculating SAF for each subspecies ==="

# Thamin (8 samples -> minInd 4)
angsd -P 64 -out ${OUT_DIR}/thamin \
  -minInd 4 -minMapQ 30 -minQ 30 \
  -bam ./thamin_bamlist.txt -rf ${CHRS} \
  -doCounts 1 -GL 1 -doSaf 1 -anc ${REF} -ref ${REF}

# Siamensis (7 samples -> minInd 4)
angsd -P 64 -out ${OUT_DIR}/siamensis \
  -minInd 4 -minMapQ 30 -minQ 30 \
  -bam ./siamensis_bamlist.txt -rf ${CHRS} \
  -doCounts 1 -GL 1 -doSaf 1 -anc ${REF} -ref ${REF}

# Hainanus (20 samples -> minInd 10)
angsd -P 64 -out ${OUT_DIR}/hainanus \
  -minInd 10 -minMapQ 30 -minQ 30 \
  -bam ./hainanus_bamlist.txt -rf ${CHRS} \
  -doCounts 1 -GL 1 -doSaf 1 -anc ${REF} -ref ${REF}

# -----------------------------------------------------------------------------
# STEP 2: CALCULATE 1D SFS PER SUBSPECIES
# -----------------------------------------------------------------------------
echo "=== Step 2: Calculating 1D SFS ==="
SUBS=("thamin" "siamensis" "hainanus")

for sub in "${SUBS[@]}"; do
    realSFS ${OUT_DIR}/${sub}.saf.idx -P 64 -fold 1 > ${OUT_DIR}/${sub}.sfs
done

# -----------------------------------------------------------------------------
# STEP 3: CALCULATE PAIRWISE 2D SFS BETWEEN SUBSPECIES
# -----------------------------------------------------------------------------
echo "=== Step 3: Calculating 2D SFS for pairs ==="
PAIRS=("thamin siamensis" "thamin hainanus" "siamensis hainanus")

for pair in "${PAIRS[@]}"; do
    pop1=$(echo $pair | cut -d' ' -f1)
    pop2=$(echo $pair | cut -d' ' -f2)
    
    echo "Processing 2D SFS for ${pop1} vs ${pop2}..."
    realSFS ${OUT_DIR}/${pop1}.saf.idx ${OUT_DIR}/${pop2}.saf.idx -P 64 > ${OUT_DIR}/${pop1}.${pop2}.ml
done

# -----------------------------------------------------------------------------
# STEP 4: INDEX PAIRWISE Fst COMPARISONS
# -----------------------------------------------------------------------------
echo "=== Step 4: Indexing pairwise Fst ==="

for pair in "${PAIRS[@]}"; do
    pop1=$(echo $pair | cut -d' ' -f1)
    pop2=$(echo $pair | cut -d' ' -f2)
    
    realSFS fst index ${OUT_DIR}/${pop1}.saf.idx ${OUT_DIR}/${pop2}.saf.idx \
      -sfs ${OUT_DIR}/${pop1}.${pop2}.ml \
      -fstout ${OUT_DIR}/out/${pop1}.${pop2} -P 64
done

# -----------------------------------------------------------------------------
# STEP 5: CALCULATE GLOBAL PAIRWISE Fst STATS
# -----------------------------------------------------------------------------
echo "=== Step 5: Calculating global Fst values ==="
rm -f ${OUT_DIR}/out/global_Fst_summary.txt

for pair in "${PAIRS[@]}"; do
    pop1=$(echo $pair | cut -d' ' -f1)
    pop2=$(echo $pair | cut -d' ' -f2)
    
    echo "--- Global Fst: ${pop1} vs ${pop2} ---" >> ${OUT_DIR}/out/global_Fst_summary.txt
    realSFS fst stats ${OUT_DIR}/out/${pop1}.${pop2}.fst.idx -P 64 >> ${OUT_DIR}/out/global_Fst_summary.txt
done

# -----------------------------------------------------------------------------
# STEP 6: SLIDING WINDOW Fst ANALYSIS (50 kb window, 25 kb step)
# -----------------------------------------------------------------------------
echo "=== Step 6: Performing sliding window analysis ==="

# Joint indexing across all 3 subspecies
realSFS fst index ${OUT_DIR}/thamin.saf.idx ${OUT_DIR}/siamensis.saf.idx ${OUT_DIR}/hainanus.saf.idx \
  -sfs ${OUT_DIR}/*.ml \
  -fstout ${OUT_DIR}/out/all_subspecies -P 64

# Sliding window calculation
realSFS fst stats2 ${OUT_DIR}/out/all_subspecies.fst.idx \
  -win 50000 -step 25000 -P 64 > ${OUT_DIR}/out/slidingwindow_all_subspecies.txt

echo "Analysis complete! Final outputs are in ${OUT_DIR}/out/"
```
