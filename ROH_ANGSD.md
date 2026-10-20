```bash
Step 0. Go to your folder and load the modules
bash
# go to the working folder
cd /home/x-bbist/Inbreeding

# make a folder for the SLURM log files
mkdir -p /home/x-bbist/Inbreeding/logs

# load the compiler and GSL library that ngsF-HMM needs
module purge
module load gcc/11.2.0 gsl/2.4
Step 1. Convert the beagle file to the ngsF-HMM format
# The beagle file was generated after running the ANGSD
# remove the first 3 columns (marker, allele1, allele2) and the header line
zcat all35.beagle.gz | cut -f 4- | awk 'NR>1' | gzip > all35.glf.gz

# take chromosome and position of each SNP from the mafs file (header removed)
zcat all35.mafs.gz | cut -f 1,2 | tail -n +2 > all35.pos
Step 2. Check the two files match

# number of SNPs in the position file
wc -l all35.pos

# number of SNPs in the likelihood file (must be the same number)
zcat all35.glf.gz | wc -l

# number of columns in the first line (must be 105 = 35 samples x 3)
zcat all35.glf.gz | head -1 | awk '{print NF}'
Step 3. Thin the SNPs (keep one SNP per 2 kb)
bash
# list the line numbers of the SNPs to keep
awk -F'\t' -v d=2000 '$1!=c || $2-l>=d {c=$1; l=$2; print NR}' all35.pos > keep.idx

# keep only those lines in the position file
awk 'NR==FNR{k[$1];next} FNR in k' keep.idx all35.pos > thin.pos

# keep only those lines in the likelihood file
awk 'NR==FNR{k[$1];next} FNR in k' keep.idx <(zcat all35.glf.gz) | gzip > thin.glf.gz
Step 4. Check the thinned files

# all three numbers must be identical (yours was 1005269)
wc -l keep.idx
wc -l thin.pos
zcat thin.glf.gz | wc -l

# SNPs per chromosome (you should see 28)
cut -f1 thin.pos | sort | uniq -c | sort -k1,1nr
Step 5. Check ngsF-HMM is built
bash
# the program should exist (you built it already)
ls -lh /home/x-bbist/Inbreeding/ngsF-HMM/ngsF-HMM
Step 6. Check sample order
# the sample order in sample.names must match the BAM order used in ANGSD
paste sample.names bam_list_35.txt | head -40
Step 7. Write the SLURM script
cat > /home/x-bbist/Inbreeding/run_ngsFHMM.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=Inbreeding_ROH
#SBATCH --account=bio260092
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=80:00:00
#SBATCH --array=1-5
#SBATCH --output=/home/x-bbist/Inbreeding/logs/ROH_%A_%a.out
#SBATCH --error=/home/x-bbist/Inbreeding/logs/ROH_%A_%a.err
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=END,FAIL

module purge
module load gcc/11.2.0 gsl/2.4

cd /home/x-bbist/Inbreeding

# count the SNPs
N=$(wc -l < /home/x-bbist/Inbreeding/thin.pos)

# run ngsF-HMM; SLURM sets $SLURM_ARRAY_TASK_ID to 1,2,3,4,5 (one seed per task)
/home/x-bbist/Inbreeding/ngsF-HMM/ngsF-HMM --geno /home/x-bbist/Inbreeding/thin.glf.gz --pos /home/x-bbist/Inbreeding/thin.pos --n_ind 35 --n_sites $N --lkl --freq r --seed $SLURM_ARRAY_TASK_ID --min_epsilon 1e-7 --min_iters 10 --max_iters 500 --n_threads 8 --out /home/x-bbist/Inbreeding/thin_roh_s$SLURM_ARRAY_TASK_ID
EOF
Step 8. Submit and monitor
# submit the 5 seeds
sbatch /home/x-bbist/Inbreeding/run_ngsFHMM.sh

Step 9. Compare the 5 seeds (after all finish)
# first line of each .indF file is the log-likelihood; the highest (closest to 0) is the best seed
head -1 thin_roh_s1.indF
head -1 thin_roh_s2.indF
head -1 thin_roh_s3.indF
head -1 thin_roh_s4.indF
head -1 thin_roh_s5.indF

# F (inbreeding) per sample for each seed: lines 2-36 of the .indF file, first column
sed -n '2,36p' thin_roh_s1.indF | cut -f1 > F_s1.txt
sed -n '2,36p' thin_roh_s2.indF | cut -f1 > F_s2.txt
sed -n '2,36p' thin_roh_s3.indF | cut -f1 > F_s3.txt
sed -n '2,36p' thin_roh_s4.indF | cut -f1 > F_s4.txt
sed -n '2,36p' thin_roh_s5.indF | cut -f1 > F_s5.txt

# put the sample names and all 5 seeds side by side
paste sample.names F_s1.txt F_s2.txt F_s3.txt F_s4.txt F_s5.txt | column -t

If the five columns are close, the run converged. If they differ a lot, rerun with --min_epsilon 1e-9.

Step 10. Check the .ibd file layout
# replace s1 with your best seed
# line count, then the number of fields on the first 5 lines
wc -l thin_roh_s1.ibd
awk '{print NR, NF}' thin_roh_s1.ibd | head -5

```
