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


#####Step 2. Use of Clumpp
```
#!/bin/bash
set -e

# Define Paths
WORKDIR=/anvil/projects/x-bio260092/Eld_der/ANGSD/Admixture
CLUMPP_DIR=$WORKDIR/CLUMPP_Linux64.1.1.2
CLUMPP_WORK=$WORKDIR/CLUMPP_run

mkdir -p $CLUMPP_WORK
cd $CLUMPP_WORK

# Step 1: Create 'pre' and 'post' helper files for 35 individuals
# DATATYPE 1 requires:
# Column 1: Individual ID (1..35)
# Column 2: Population ID (1)
# Last Column: End marker (1)
rm -f pre post
for i in $(seq 1 35); do
    echo -e "${i}\t1" >> pre
    echo -e "1" >> post
done

# Step 2: Loop through K values 1 to 10
for K in $(seq 1 10); do
    echo "=========================================="
    echo "Processing K = ${K}"
    echo "=========================================="

    # Remove temporary files from previous runs
    rm -f tmp_K${K}_R*.Q pop_K${K}_combined.txt

    # Format .qopt files from all 99 replicate runs into individual .Q files
    for r in $(seq 1 99); do
        QFILE=$WORKDIR/ADX/R_${r}/K${K}_R${r}.qopt
        if [ -s "$QFILE" ]; then
            paste pre "$QFILE" post > tmp_K${K}_R${r}.Q
        else
            echo "Warning: File $QFILE is missing or empty!"
        fi
    done

    # Concatenate all 99 formatted replicate matrices into ONE combined popfile for this K
    cat tmp_K${K}_R*.Q > pop_K${K}_combined.txt

    # Verify that the combined popfile is not empty before invoking CLUMPP
    if [ ! -s pop_K${K}_combined.txt ]; then
        echo "Error: Combined popfile pop_K${K}_combined.txt is empty. Skipping K=${K}."
        continue
    fi

    # Write the paramfile specifically for this K
    cat << EOF > paramfile_K${K}
DATATYPE 1
INDFILE NOTNEEDED.indfile 
POPFILE pop_K${K}_combined.txt 
OUTFILE pop_K${K}-combined-merged.txt 
MISCFILE pop_K${K}-combined-miscfile.txt 
K ${K} 
C 35 
R 99 
M 2 
W 0 
S 2 
GREEDY_OPTION 2 
REPEATS 20 
PERMUTATIONFILE NOTNEEDED.permutationfile 
PRINT_PERMUTED_DATA 1 
PERMUTED_DATAFILE pop_K${K}-combined-aligned.txt 
PRINT_EVERY_PERM 0 
EVERY_PERMFILE pop_K${K}-combined.every_permfile 
PRINT_RANDOM_INPUTORDER 0 
RANDOM_INPUTORDERFILE pop_K${K}-combined.random_inputorderfile 
OVERRIDE_WARNINGS 0 
ORDER_BY_RUN 0 
EOF

    # Step 3: Run CLUMPP for the current K value
    $CLUMPP_DIR/CLUMPP paramfile_K${K}

    # Clean up temporary files for this K
    rm -f tmp_K${K}_R*.Q

    echo "Completed K = ${K}"
done

echo "All K values successfully processed!"
```
