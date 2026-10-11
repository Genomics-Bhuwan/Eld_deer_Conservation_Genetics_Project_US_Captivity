```bash
#!/bin/bash
#SBATCH --job-name=pcangsd_hwe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --output=/anvil/scratch/x-bbist/Eld_Deer/ANGSD_new/out/PCANGSD_for_PCA/pcangsd_%j.out
#SBATCH --error=/anvil/scratch/x-bbist/Eld_Deer/ANGSD_new/out/PCANGSD_for_PCA/pcangsd_%j.err

# Load Conda and activate environment
module load conda/2026.03
conda activate /home/x-bbist/.conda/envs/pcangsd_env

# Define directories
WORKDIR=/anvil/scratch/x-bbist/Eld_Deer/ANGSD_new/out
OUTDIR=${WORKDIR}/PCANGSD_for_PCA
BEAGLE_FILE=${WORKDIR}/all35.beagle.gz
PCANGSD_BIN=/home/x-bbist/.conda/envs/pcangsd_env/bin/pcangsd

mkdir -p "$OUTDIR"
cd "$WORKDIR" || exit 1

# Execution parameters
THREADS=5
E_VARS=3

# Check input file
if [[ ! -f "$BEAGLE_FILE" ]]; then
    echo "ERROR: Beagle file not found: $BEAGLE_FILE"
    exit 1
fi

# Step 1: Inbreeding estimation and site-level analysis
"$PCANGSD_BIN" \
    -b "$BEAGLE_FILE" \
    --inbreed-sites \
    -e "$E_VARS" \
    -t "$THREADS" \
    --sites-save \
    -o "$OUTDIR/hwe_global_loci_e${E_VARS}"

# Stop if Step 1 fails
if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 1 failed."
    exit 1
fi

# Step 2: HWE filtering and covariance matrix analysis
"$PCANGSD_BIN" \
    -b "$BEAGLE_FILE" \
    -e "$E_VARS" \
    -t "$THREADS" \
    --hwe "$OUTDIR/hwe_global_loci_e${E_VARS}.lrt.sites.npy" \
    --hwe-tole 1e-6 \
    --sites-save \
    -o "$OUTDIR/pcangsd.HWE.global.loci.results"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 2 failed."
    exit 1
fi

echo "PCAngsd analysis completed successfully."
```
