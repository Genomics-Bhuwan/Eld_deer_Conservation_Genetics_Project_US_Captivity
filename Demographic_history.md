VCF file
  │
  ▼
Step 1: Split VCF by subspecies (3 groups)
  │
  ▼
Step 2: Extract SFS per subspecies using easySFS or vcftools
  │
  ▼
Step 3: Create blueprint (.blueprint) files for Stairway Plot 2
  │
  ▼
Step 4: Build and run Stairway Plot 2
  │
  ▼
Step 5: Visualize demographic history plots


#### Step 1. Create the text file for the sub-species
```bash
# Make your working directory
~/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2

# Create a population map file (sample → population) ──
# This is needed by easySFS and vcftools to know which samples belong where

cat > popmap.txt << 'EOF'
RE116077    thamin
RE115445    thamin
RE115304    thamin
RE115125    thamin
RE114623    thamin
RE114379    thamin
RE113838    thamin
RE115604    thamin
CRR773117   hainanus
CRR773116   hainanus
CRR773115   hainanus
CRR773114   hainanus
CRR773113   hainanus
CRR773112   hainanus
CRR773111   hainanus
CRR773110   hainanus
CRR773109   hainanus
CRR773108   hainanus
CRR773107   hainanus
CRR773106   hainanus
CRR773105   hainanus
CRR773104   hainanus
CRR773103   hainanus
CRR773102   hainanus
CRR773101   hainanus
CRR773100   hainanus
CRR773099   hainanus
CRR773098   hainanus
CRR773097   siamensis
CRR773096   siamensis
CRR773095   siamensis
CRR773094   siamensis
CRR773093   siamensis
CRR773091   siamensis
CRR773090   siamensis
EOF

# Also create per-subspecies sample lists (needed for vcftools --keep)
grep "thamin"   popmap.txt | awk '{print $1}' > samples_thamin.txt
grep "hainanus" popmap.txt | awk '{print $1}' > samples_hainanus.txt
grep "siamensis" popmap.txt | awk '{print $1}' > samples_siamensis.txt
```

#### Step 2: Generate the SFS using easySFS (recommended)
- easySFS is the most reliable VCF→SFS tool for Stairway Plot 2.
- Install easySFS if not installed

```bash
cd ~
git clone https://github.com/isaacovercast/easySFS.git
cd easySFS
pip install -r requirements.txt   # or conda install

# ── Set your VCF path ──
VCF=/home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf

POPMAP=~/Demographic_history/Stairway_Plot2/per_subspecies/popmap.txt

# ── FIRST: Preview projection sizes ──
# This tells you how many segregating sites you keep at each projection level.
# Always run this BEFORE the real run to choose the best projection.
python3 easySFS.py \
    -i $VCF \
    -p $POPMAP \
    --preview \
    -a          # -a = use all sites (not just variable); remove if you want SNPs only

# ── THEN: Run with chosen projections ──
# Rule of thumb: project to ~80% of max haplotypes to balance sites vs. sample size
python3 easySFS/easySFS.py \
     -i Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf \
     -p popmap.txt \
     --proj 12,4,4 \
     --ploidy 2 \
     -a
` The output SFS files will be in:

~/Demographic_history/Stairway_Plot2/easySFS_output/fastsimcoal2/

as thamin_MSFS.obs, hainanus_MSFS.obs, siamensis_MSFS.obs

Step 3: Get the total number of callable sites (L)
Stairway Plot 2 requires L = total number of sites sequenced (not just SNPs). This is critical.
bash# Option A: If you have a callable sites BED file from your variant calling pipeline
awk '{sum += $3 - $2} END {print sum}' callable_sites.bed

# Option B: Count sites from the VCF (less accurate, use only if no BED available)
# Total SNPs is NOT L — L is the whole surveyed genome length

# Option C: Use the genome size minus unsequenced regions
# If you filtered to only certain chromosomes/scaffolds, sum their lengths from the .fai index:
VCF=/home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf

# Get reference contigs from VCF header
grep "^##contig" $VCF | sed 's/.*length=\([0-9]*\).*/\1/' | awk '{sum+=$1} END{print "Total L =", sum}'

# ⚠️  Write down your L value — you need it in Step 4

Step 4: Read the SFS values from easySFS output
bash# Look at the SFS output for each population
# The folded SFS from easySFS looks like:
# 1 entry per frequency class, space-separated

cat ~/Demographic_history/Stairway_Plot2/easySFS_output/fastsimcoal2/thamin/thamin_MSFS.obs
# e.g.: 1 observation
#       0  452  310  198  ... (counts at freq 0, 1, 2, ...)
# The first value (freq=0, monomorphic) is usually dropped for Stairway Plot

# For Stairway Plot 2, you need the 1D SFS for each population separately.
# Extract just that population's marginal SFS.

Step 5: Install Stairway Plot 2

cd ~/Demographic_history/Stairway_Plot2/

# Download (already have zip from GitHub)
# If not downloaded yet:
wget https://github.com/xiaoming-liu/stairway-plot-v2/raw/master/stairway_plot_v2.1.3.zip
unzip stairway_plot_v2.1.3.zip
cd stairway_plot_v2.1.3

# Verify Java is installed (required!)
java -version   # needs Java 8 or higher

# Compile (if .class files not included)
javac *.java

Step 6: Create blueprint files for each subspecies
This is the core configuration file for Stairway Plot 2. Create one per subspecies.
cd /home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/

# 1. Update thamin.blueprint (6 numbers)
cat << 'EOF' > thamin.blueprint
# Main parameters
popid: thamin
n_samples: 12
sequence_length: 2500000000
mutation_rate: 2.11e-9
generation_time: 1.5

# SFS configurations
n_snps: 354967
sfs: 78515 70801 64285 57042 55902 28422

# Run settings
n_est: 200
n_rand: 15 30 45 60
project_dir: thamin_SP2_run
stairway_plot_dir: /home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/stairway_plot_v2.1.3
EOF

# 2. Update hainanus.blueprint (2 numbers)
cat << 'EOF' > hainanus.blueprint
# Main parameters
popid: hainanus
n_samples: 4
sequence_length: 2500000000
mutation_rate: 2.11e-9
generation_time: 1.5

# SFS configurations
n_snps: 90648
sfs: 34755 55893

# Run settings
n_est: 200
n_rand: 2 4 6 8
project_dir: hainanus_SP2_run
stairway_plot_dir: /home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/stairway_plot_v2.1.3
EOF

# 3. Update siamensis.blueprint (2 numbers)
cat << 'EOF' > siamensis.blueprint
# Main parameters
popid: siamensis
n_samples: 4
sequence_length: 2500000000
mutation_rate: 2.11e-9
generation_time: 1.5

# SFS configurations
n_snps: 25911
sfs: 3254 22658

# Run settings
n_est: 200
n_rand: 2 4 6 8
project_dir: siamensis_SP2_run
stairway_plot_dir: /home/bistbs/Elds_deer_Population_Genetic_Analysis/Demographic_history/Stairway_Plot2/stairway_plot_v2.1.3
EOF

# Verify that all 3 files were written successfully
ls -l *.blueprint
⚠️ The SFS values above are placeholders. You MUST replace them with the actual counts from your easySFS output in Step 2.

### Step 7: Build the run scripts with Stairbuilder
bash
cd ~/Demographic_history/Stairway_Plot2/stairway_plot_v2.1.3/

# This generates a shell script for each blueprint
java -cp stairway_plot_v2.1.3.jar Stairbuilder thamin.blueprint
java -cp stairway_plot_v2.1.3.jar Stairbuilder hainanus.blueprint
java -cp stairway_plot_v2.1.3.jar Stairbuilder siamensis.blueprint

# Each command produces a .sh file, e.g. thamin.blueprint.sh
ls *.blueprint.sh

Step 8: Run Stairway Plot 2
bash# ── Run each subspecies (can run in parallel in separate terminals) ──

# Option A: Sequential (safe on single machine)
bash thamin.blueprint.sh   2>&1 | tee thamin_run.log
bash hainanus.blueprint.sh 2>&1 | tee hainanus_run.log
bash siamensis.blueprint.sh 2>&1 | tee siamensis_run.log

# Option B: Parallel background jobs (if you have enough cores)
bash thamin.blueprint.sh   > thamin_run.log 2>&1 &
bash hainanus.blueprint.sh > hainanus_run.log 2>&1 &
bash siamensis.blueprint.sh > siamensis_run.log 2>&1 &
wait  # waits for all background jobs to finish

# ── Monitor progress ──
tail -f thamin_run.log      # Ctrl+C to stop watching
Running time: ~30 min to a few hours depending on genome size and server speed.

Step 9: Generate the final plots
bash# After Step 8 completes, generate the demographic history plots
# Stairway Plot 2 auto-generates plots from the summary files

# Check output files exist
ls thamin_output/
# Should contain: thamin.final.summary  thamin.*.addTheta  etc.

# Generate plots (this uses Stairpainter internally via the blueprint.sh)
# If plots weren't auto-generated, run:
java -cp stairway_plot_v2.1.3.jar Stairpainter thamin.blueprint
java -cp stairway_plot_v2.1.3.jar Stairpainter hainanus.blueprint
java -cp stairway_plot_v2.1.3.jar Stairpainter siamensis.blueprint

# Output plots will be in each project_dir as .png and .pdf files
ls thamin_output/*.png
ls hainanus_output/*.png
ls siamensis_output/*.png
