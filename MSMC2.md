#### DOWNLOAD MSMC2
- git clone https://github.com/stschiff/msmc-tools.git 

#### Step 1. Phase the multi-sample vcf using BEAGLE.
- This step seperates the maternal and paternal alleles in the vcf.
```bash
java -Xmx24g -jar beagle5.jar gt=Eld_Deer_35_samples_GQ20_test_biallelic.recode.vcf out=Eld_Deer_phased

## Index the multi-sample phased vcf
bcftools index Eld_Deer_phased.vcf.gz

#### Step 2. Generate the indiviudal genomic mask files. 
- This loop generates the individual callability masks for the samples so MSMC2 knows which regions have reliable data.
```bash
for sample in $(bcftools query -l Eld_Deer_phased.vcf.gz); do
    echo "Processing mask for sample: $sample"
    bcftools view -s $sample Eld_Deer_phased.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\n' | \
    awk 'BEGIN{c="";s=0;e=0} {if($1!=c || $2!=e+1){if(c!="")print c"\t"s"\t"e; c=$1;s=$2;e=$2}else{e=$2}} END{print c"\t"s"\t"e}' | \
    bgzip -c > ${sample}.mask.bed.gz
    tabix -p bed ${sample}.mask.bed.gz
done
```
#### Step 3: Convert Sample Subspecies to Multihetsep Formats
- Run the python utility script to compile the specialized input text files for each of the three subspecies group.
# 1. Generate Thamin input
# 1. Generate Thamin Population Input
python3 msmc-tools/generate_multihetsep.py \
    --mask RE116077.mask.bed.gz --mask RE115445.mask.bed.gz \
    --mask RE115304.mask.bed.gz --mask RE115125.mask.bed.gz \
    Eld_Deer_phased.vcf.gz > thamin.multihetsep.txt

# 2. Generate Hainanus Population Input
python3 msmc-tools/generate_multihetsep.py \
    --mask CRR773117.mask.bed.gz --mask CRR773116.mask.bed.gz \
    --mask CRR773115.mask.bed.gz --mask CRR773114.mask.bed.gz \
    Eld_Deer_phased.vcf.gz > hainanus.multihetsep.txt

# 3. Generate Siamensis Population Input
python3 msmc-tools/generate_multihetsep.py \
    --mask CRR773097.mask.bed.gz --mask CRR773096.mask.bed.gz \
    --mask CRR773095.mask.bed.gz --mask CRR773094.mask.bed.gz \
    Eld_Deer_phased.vcf.gz > siamensis.multihetsep.txt

- Finally, run the calculation engine sequentially for each subspecies. The -t 4 flag tells it to utilize 4 CPU cores on your login node:

Bash
# Run Thamin calculation
msmc2 -t 4 -p 1*2+15*1+1*2 -o thamin_msmc2 -I 0,1,2,3,4,5,6,7 thamin.multihetsep.txt

# Run Hainanus calculation
msmc2 -t 4 -p 1*2+15*1+1*2 -o hainanus_msmc2 -I 0,1,2,3,4,5,6,7 hainanus.multihetsep.txt

# Run Siamensis calculation
msmc2 -t 4 -p 1*2+15*1+1*2 -o siamensis_msmc2 -I 0,1,2,3,4,5,6,7 siamensis.multihetsep.txt
