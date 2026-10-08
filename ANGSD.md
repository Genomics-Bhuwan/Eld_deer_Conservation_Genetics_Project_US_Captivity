#### Testing the bam and listing it
```bash
cd /anvil/scratch/x-bbist/Eld_Deer/ANGSD_new
module load samtools

# BAM list
ls $PWD/CRR*_final_dedup.bam $PWD/RE*_10X.bam > bam_list_35.txt
wc -l bam_list_35.txt                       # must be 35

# no output = all indexes present
for b in $(cat bam_list_35.txt); do [ -f ${b}.bai ] || echo "missing index: $b"; done

# expect one line: "35 SN:CM139490.1"
for b in $(cat bam_list_35.txt); do samtools view -H $b | grep -m1 '^@SQ' | cut -f2; done | sort | uniq -c

# check the BAMs are real files, not broken links
ls -lL *.bam | wc -l                        # expect 35

# scaffolds >= 1 Mb (names only, for -rf)
awk '$2>=1000000{print $1}' mRucEld1.hap1_genomic.fna.fai > scaffolds_1Mb.txt
wc -l scaffolds_1Mb.txt
awk '$2>=1000000{s+=$2} END{print s/1e9 " Gb in scaffolds >= 1 Mb"}' mRucEld1.hap1_genomic.fna.fai
```

#### Step 2. Batch script for analysis of ANGSD

```bash
cd /anvil/scratch/x-bbist/Eld_Deer/ANGSD_new
mkdir -p logs out

cat > angsd_all35.sbatch <<'EOF'
#!/bin/bash
#SBATCH -J angsd_all35
#SBATCH -A YOUR_ALLOCATION
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -o logs/angsd_all35_%j.out
#SBATCH -e logs/angsd_all35_%j.err

source ~/.bashrc
conda activate angsd_env

cd /anvil/scratch/x-bbist/Eld_Deer/ANGSD_new

angsd -b bam_list_35.txt -ref mRucEld1.hap1_genomic.fna -rf scaffolds_autosomes.txt \
  -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -skipTriallelic 1 \
  -minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 \
  -minInd 25 -setMinDepthInd 1 -setMaxDepth 280 \
  -doCounts 1 -doDepth 1 -maxDepth 1000 \
  -doIBS 1 -makeMatrix 1 -doCov 1 \
  -nThreads 16 \
  -out out/all35
EOF

grep -c "YOUR_ALLOCATION" angsd_all35.sbatch    # must print 0 after you edit it
grep -c "scaffolds_autosomes.txt" angsd_all35.sbatch   # expect 1
sbatch angsd_all35.sbatch
```
