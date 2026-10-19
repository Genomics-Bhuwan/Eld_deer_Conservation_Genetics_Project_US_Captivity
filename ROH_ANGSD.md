```bash
#!/bin/bash
#SBATCH --job-name=Inbreeding_ROH
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=65
#SBATCH --mem=100G
#SBATCH --time=80:00:00
#SBATCH --partition=shared
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -e

# ---------------- PATHS ----------------
BAMDIR=/anvil/projects/x-bio260092/Eld_der
WORKDIR=/anvil/projects/x-bio260092/Eld_der/ANGSD/Inbreeding
REFGZ=$WORKDIR/mRucEld1.hap1_genomic.fna.gz
REF=$WORKDIR/mRucEld1.hap1_genomic.fna
OUT=$WORKDIR/eld_roh

cd $WORKDIR

# ---------------- AUTOSOMES (written here) ----------------
cat > $WORKDIR/autosomes_roh.txt << 'EOF'
CM139490.1
CM139491.1
CM139492.1
CM139493.1
CM139494.1
CM139495.1
CM139496.1
CM139497.1
CM139498.1
CM139499.1
CM139500.1
CM139501.1
CM139502.1
CM139503.1
CM139504.1
CM139505.1
CM139506.1
CM139507.1
CM139508.1
CM139509.1
CM139510.1
CM139511.1
CM139512.1
CM139513.1
CM139514.1
CM139515.1
CM139516.1
CM139517.1
EOF

# ---------------- BAM LIST (full paths, written here) ----------------
cat > $WORKDIR/bam_list_35_full.txt << EOF
$BAMDIR/CRR773090_final_dedup.bam
$BAMDIR/CRR773091_final_dedup.bam
$BAMDIR/CRR773093_final_dedup.bam
$BAMDIR/CRR773094_final_dedup.bam
$BAMDIR/CRR773095_final_dedup.bam
$BAMDIR/CRR773096_final_dedup.bam
$BAMDIR/CRR773097_final_dedup.bam
$BAMDIR/CRR773098_final_dedup.bam
$BAMDIR/CRR773099_final_dedup.bam
$BAMDIR/CRR773100_final_dedup.bam
$BAMDIR/CRR773101_final_dedup.bam
$BAMDIR/CRR773102_final_dedup.bam
$BAMDIR/CRR773103_final_dedup.bam
$BAMDIR/CRR773104_final_dedup.bam
$BAMDIR/CRR773105_final_dedup.bam
$BAMDIR/CRR773106_final_dedup.bam
$BAMDIR/CRR773107_final_dedup.bam
$BAMDIR/CRR773108_final_dedup.bam
$BAMDIR/CRR773109_final_dedup.bam
$BAMDIR/CRR773110_final_dedup.bam
$BAMDIR/CRR773111_final_dedup.bam
$BAMDIR/CRR773112_final_dedup.bam
$BAMDIR/CRR773113_final_dedup.bam
$BAMDIR/CRR773114_final_dedup.bam
$BAMDIR/CRR773115_final_dedup.bam
$BAMDIR/CRR773116_final_dedup.bam
$BAMDIR/CRR773117_final_dedup.bam
$BAMDIR/RE113838_10X.bam
$BAMDIR/RE114379_10X.bam
$BAMDIR/RE114623_10X.bam
$BAMDIR/RE115125_10X.bam
$BAMDIR/RE115304_10X.bam
$BAMDIR/RE115445_10X.bam
$BAMDIR/RE115604_10X.bam
$BAMDIR/RE116077_10X.bam
EOF

# sanity checks
echo "BAMs listed: $(wc -l < bam_list_35_full.txt) (expect 35)"
echo "Autosomes:   $(wc -l < autosomes_roh.txt) (expect 28)"
for b in $(cat bam_list_35_full.txt); do
    [ -s "$b" ] || { echo "MISSING BAM: $b"; exit 1; }
done

# ---------------- REFERENCE (unzip + index if needed) ----------------
module load biocontainers
module load samtools
if [ ! -s "$REF" ]; then
    zcat $REFGZ > $REF
fi
if [ ! -s "$REF.fai" ] || [ "$REF" -nt "$REF.fai" ]; then
    samtools faidx $REF
fi

# ---------------- 1. ANGSD ----------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate angsd_env

angsd -GL 1 -dobcf 1 -dopost 1 -domajorminor 1 -domaf 1 \
      -minQ 30 -minMapQ 30 -minInd 28 \
      -SNP_pval 1e-6 \
      -rf $WORKDIR/autosomes_roh.txt \
      -P 64 \
      -bam $WORKDIR/bam_list_35_full.txt \
      -ref $REF \
      -out $OUT

conda deactivate

# ---------------- 2. Allele frequencies ----------------
module load biocontainers
module load bcftools
module load htslib
module load samtools

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' ${OUT}.bcf | bgzip -c > ${OUT}.freqs.tab.gz
tabix -s1 -b2 -e2 ${OUT}.freqs.tab.gz

# ---------------- 3.# HMM-based ROH detection using genotype likelihoods (PL) and the AF prior
bcftools roh --AF-file ${OUT}.freqs.tab.gz \
             --output ${OUT}_ROH_raw.txt \
             --threads 8 \
             ${OUT}.bcf

# ---------------- 4. Keep only "RG" (ROH region) records from the bcftools roh output
awk '$1=="RG"' ${OUT}_ROH_raw.txt > ${OUT}_ROH_RG_all.txt

# ---------------- 5. Per-sample split ----------------
bcftools query -l ${OUT}.bcf > $WORKDIR/sample.names
mkdir -p $WORKDIR/per_sample
for i in $(cat $WORKDIR/sample.names); do
    awk -v s="$i" '$2==s' ${OUT}_ROH_RG_all.txt > $WORKDIR/per_sample/$(basename $i).ROH.txt
done

# ---------------- 6. Summary per individual ----------------
awk 'BEGIN{OFS="\t"; print "sample","n_ROH","total_ROH_bp"}
     {n[$2]++; len[$2]+=$6} END{for(s in n) print s, n[s], len[s]}' \
     ${OUT}_ROH_RG_all.txt > ${OUT}_ROH_summary.txt
```
