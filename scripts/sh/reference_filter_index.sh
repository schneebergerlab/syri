cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data
# Get reference genome
wget  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz

# index genome
samtools faidx hg38.fa.gz

# Filter alt and fix sequences
grep -P  'alt|fix' hg38.fa.gz.fai | cut -f1 > chr_to_filt.txt
hometools getchr -v -F chr_to_filt.txt -o hg38.filt.fa.gz hg38.fa.gz

# Get autosomal chromosomes
rm -f chromosomes.txt
for i in {1..22}; do
  echo chr$i >> chromosomes.txt
done
hometools getchr -F chromosomes.txt -o hg38.only_chromosome.fa.gz hg38.fa.gz &

# Get indices
samtools faidx hg38.filt.fa.gz
samtools faidx hg38.only_chromosome.fa.gz

bowtie2-build --threads 40 hg38.filt.fa.gz hg38.filt.fa
picard CreateSequenceDictionary -R hg38.filt.fa.gz

# Liftover GIAB_SVs from GRCh37 to GRCh38 using liftover
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh38
rm -f chr_name_maps.txt
for i in {1..22} X Y; do
  echo -e $i'\t'chr$i >> chr_name_maps.txt
done
bcftools annotate --rename-chrs chr_name_maps.txt HG002_SVs_Tier1_v0.6.vcf.gz | bgzip > HG002_SVs_Tier1_v0.6.renamed.vcf.gz
gatk LiftoverVcf \
-C hg19ToHg38.over.chain.gz \
-I HG002_SVs_Tier1_v0.6.renamed.vcf.gz \
-O hg38.HG002_SVs_Tier1_v0.6.vcf.gz \
-R ../../../hg38.fa.gz \
--REJECT rejected_variants.vcf

#-----------------------------------------------------------------------------------------------------------------------
## The liftover based VCF was very different from variant calles suggesting issues in lifting over. Using GRCh37/hg19 instead
# Get genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
samtools faidx hg19.fa.gz

# The alternate haplotypes have 'hap' string for hg19
grep -P  'hap' hg19.fa.gz.fai | cut -f1 > hg19_chr_to_filt.txt
hometools getchr -v -F hg19_chr_to_filt.txt -o hg19.filt.fa.gz hg19.fa.gz

# Get autosomal chromosomes
hometools getchr -F chromosomes.txt -o hg19.only_chromosome.fa hg19.fa &

# Get indices
samtools faidx hg19.filt.fa.gz
samtools faidx hg19.only_chromosome.fa.gz

bowtie2-build --threads 40 hg19.filt.fa.gz hg19.filt.fa &
picard CreateSequenceDictionary -R hg19.filt.fa.gz &


