```bash
#!/bin/bash
#SBATCH --job-name=NGS_admix
#SBATCH --output=ngsadmix_%j.out
#SBATCH --error=ngsadmix_%j.err
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
WORKDIR=/anvil/projects/x-bio260092/Eld_der/ANGSD/Admixture
BEAGLE=$WORKDIR/all35_thin10kb.beagle.gz
OUTDIR=$WORKDIR/ADX

# ---------------- SETTINGS ----------------
KMIN=1
KMAX=10
NREPS=99          # replicates per K
NJOBS=16          # NGSadmix runs in parallel
THREADS=4         # threads per run (NJOBS x THREADS = 32 = -n above)

cd $WORKDIR

# ---------------- CHECKS ----------------
[ -x ~/.conda/envs/angsd_env/bin/NGSadmix ] || { echo "NGSadmix not found"; exit 1; }
[ -s "$BEAGLE" ] || { echo "Missing $BEAGLE"; exit 1; }
echo "Sites in thinned file: $(( $(zcat $BEAGLE | wc -l) - 1 ))"

# ---------------- RUN NGSadmix ----------------
mkdir -p $OUTDIR
for j in $(seq 1 $NREPS); do mkdir -p $OUTDIR/R_${j}; done

for j in $(seq 1 $NREPS); do
  for k in $(seq $KMIN $KMAX); do
    echo "$k $j"
  done
done | xargs -P $NJOBS -L 1 bash -c '
  k=$0; j=$1
  ~/.conda/envs/angsd_env/bin/NGSadmix -likes '"$BEAGLE"' -K $k -P '"$THREADS"' \
           -minMaf 0.05 -maxiter 50000 -tol 1e-9 -tolLike50 1e-9 \
           -seed $((k*1000 + j)) \
           -outfiles '"$OUTDIR"'/R_${j}/K${k}_R${j} > /dev/null 2>&1
  echo "done K=$k rep=$j"
'

# ---------------- LIKELIHOOD SUMMARY ----------------
echo -e "K\trep\tlikelihood" > $OUTDIR/likelihoods.txt
for j in $(seq 1 $NREPS); do
  for k in $(seq $KMIN $KMAX); do
    LOG=$OUTDIR/R_${j}/K${k}_R${j}.log
    L=$(grep -o 'best like=[-0-9.]*' $LOG | head -1 | cut -d= -f2)
    echo -e "${k}\t${j}\t${L}" >> $OUTDIR/likelihoods.txt
  done
done

echo "Finished. Results in $OUTDIR"
```
