########################################################################################################################
## Pre-processing
########################################################################################################################
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa

#-----------------------------------------------------------------------------------------------------------------------
# GIAB Benchmark

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh37/

# BED file defines regions in which HG002_SVs_Tier1_v0.6.vcf.gz should contain close to all of true insertions and deletions >=50bp. This is more conservative than HG002_SVs_Tier1_v0.6.bed in that it excludes the VDJ and X and Y. (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/README_SV_v0.6.txt)
#conbed=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh37/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed
conbed=/home/ra98jam/d16/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh37/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed

sed -i 's/^/chr/g' HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed

# Rename chromosomes in VCF file
for i in {1..22} X Y; do
  echo -e $i'\t'chr$i >> chr_name_maps.txt
done
bcftools annotate --rename-chrs chr_name_maps.txt HG002_SVs_Tier1_v0.6.vcf.gz > HG002_SVs_Tier1_v0.6.renamed.vcf

# Get only Tier1
vcftools --vcf HG002_SVs_Tier1_v0.6.renamed.vcf --recode --recode-INFO-all --bed $conbed --out HG002_SVs_Tier1_v0.6.onlyTier1

# Get passed variants
vcftools --vcf HG002_SVs_Tier1_v0.6.onlyTier1.recode.vcf --remove-filtered-all --recode --recode-INFO-all --bed $conbed --out HG002_SVs_Tier1_v0.6.passed

# Filter SVs smaller than 50 bps
SURVIVOR filter HG002_SVs_Tier1_v0.6.passed.recode.vcf NA 50 -1 0 -1 HG002_SVs_Tier1_v0.6.passed.large.vcf


#-----------------------------------------------------------------------------------------------------------------------
# SYRI

# Hg38 -----------------------------------------------------------------------------------------------------------------

# Merge Syri calls
cd /home/ra98jam/d16/projects/sv_benchmark/results/variant_calls/syri # linked to dss-0016
/home/ra98jam/d16/projects/syri/scripts/python/syri/syri/scripts/vcfasm hapmerge hg002.patsyri.vcf hg002.matsyri.vcf tmp.vcf hg002
/home/ra98jam/d16/projects/syri/scripts/python/syri/syri/scripts/vcfasm hapmerge hg002.nofilt.patsyri.vcf hg002.nofilt.matsyri.vcf tmp1.vcf hg002
# remove 'pat' sample from the merged file
bcftools view -s ^pat tmp.vcf >  hg002.syri.vcf
bcftools view -s ^pat tmp1.vcf >  hg002.nofilt.syri.vcf

sed -i 's/VCFv4.3/VCFv4.2/g' hg002.syri.vcf
sed -i 's/VCFv4.3/VCFv4.2/g' hg002.nofilt.syri.vcf
# SNPs
vcftools --vcf hg002.syri.vcf --remove-indels --recode --recode-INFO-all --bed $conbed --out hg002.syri.SNPs
vcftools --vcf hg002.nofilt.syri.vcf --remove-indels --recode --recode-INFO-all --bed $conbed --out hg002.nofilt.syri.SNPs
# Indels
vcftools --vcf hg002.syri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg002.syri.indels
vcftools --vcf hg002.nofilt.syri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg002.nofilt.syri.indels

# Hg19 -----------------------------------------------------------------------------------------------------------------
# Merged/diploid
# Merge VCFs using survivor
SURVIVOR merge survivor.in 50 1 1 1 0 50 sample_merged.vcf # distance 50
SURVIVOR merge survivor.nofilt.in 50 1 1 1 0 50 sample_merged.nofilt.vcf

# Create a single sample file with diploid GT values based on the SUPP_VEC info
/home/ra98jam/d16/projects/sv_benchmark/scripts/py/add_genotype_survivor_merged.py sample_merged.vcf hg19.hg002.syri.vcf
/home/ra98jam/d16/projects/sv_benchmark/scripts/py/add_genotype_survivor_merged.py sample_merged.nofilt.vcf hg19.hg002.nofilt.syri.vcf

# Indels
vcftools --vcf hg19.hg002.syri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.syri.indels
vcftools --vcf hg19.hg002.nofilt.syri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.nofilt.syri.indels
# select large indels (filter out smaller indels)
SURVIVOR filter hg19.hg002.syri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.syri.indels.large.vcf
SURVIVOR filter hg19.hg002.nofilt.syri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.nofilt.syri.indels.large.vcf

# Haplotype-paternal
sed -i 's/VCFv4.3/VCFv4.2/g' hg19.hg002.patsyri.vcf
vcftools --vcf hg19.hg002.patsyri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.patsyri.indels
SURVIVOR filter hg19.hg002.patsyri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.patsyri.indels.large.vcf
awk '/^chr/ {print $0"|1"; next} 1' hg19.hg002.patsyri.indels.large.vcf > hg19.hg002.patsyri.indels.large.dip_gt.vcf

sed -i 's/VCFv4.3/VCFv4.2/g' hg19.hg002.nofilt.patsyri.vcf
vcftools --vcf hg19.hg002.nofilt.patsyri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.nofilt.patsyri.indels
SURVIVOR filter hg19.hg002.nofilt.patsyri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.nofilt.patsyri.indels.large.vcf
awk '/^chr/ {print $0"|1"; next} 1' hg19.hg002.nofilt.patsyri.indels.large.vcf > hg19.hg002.nofilt.patsyri.indels.large.dip_gt.vcf


# Haplotype-maternal
sed -i 's/VCFv4.3/VCFv4.2/g' hg19.hg002.matsyri.vcf
vcftools --vcf hg19.hg002.matsyri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.matsyri.indels
SURVIVOR filter hg19.hg002.matsyri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.matsyri.indels.large.vcf
awk '/^chr/ {print $0"|1"; next} 1' hg19.hg002.matsyri.indels.large.vcf > hg19.hg002.matsyri.indels.large.dip_gt.vcf

sed -i 's/VCFv4.3/VCFv4.2/g' hg19.hg002.nofilt.matsyri.vcf
vcftools --vcf hg19.hg002.nofilt.matsyri.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.nofilt.matsyri.indels
SURVIVOR filter hg19.hg002.nofilt.matsyri.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.nofilt.matsyri.indels.large.vcf
awk '/^chr/ {print $0"|1"; next} 1' hg19.hg002.nofilt.matsyri.indels.large.vcf > hg19.hg002.nofilt.matsyri.indels.large.dip_gt.vcf


#-----------------------------------------------------------------------------------------------------------------------
# svim-asm

# Hg19 -----------------------------------------------------------------------------------------------------------------
# merged
vcftools --vcf hg19.hg002_dip/variants.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.dip.indels
SURVIVOR filter hg19.hg002.dip.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.dip.indels.large.vcf

# Haplotype-paternal
vcftools --vcf hg19.hg002_pat/variants.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.pat.indels
SURVIVOR filter hg19.hg002.pat.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.pat.indels.large.vcf

# Haplotype-maternal
vcftools --vcf hg19.hg002_mat/variants.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.mat.indels
SURVIVOR filter hg19.hg002.mat.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.mat.indels.large.vcf



#-----------------------------------------------------------------------------------------------------------------------
# Sniffles2

# Hg19 -----------------------------------------------------------------------------------------------------------------
vcftools --gzvcf hg19.hg002.vcf.gz --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.indels
SURVIVOR filter hg19.hg002.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.indels.large.vcf



#-----------------------------------------------------------------------------------------------------------------------
# cutesv
# Hg19 -----------------------------------------------------------------------------------------------------------------
# Hifi
vcftools --gzvcf hg19.hg002.vcf.gz --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.indels
SURVIVOR filter hg19.hg002.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.indels.large.vcf
# diploid assembly
vcftools --vcf hg19.hg002.dip_asm.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.diploid.indels
SURVIVOR filter hg19.hg002.diploid.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.diploid.indels.large.vcf



#-----------------------------------------------------------------------------------------------------------------------
# pbsv

# Hg19 -----------------------------------------------------------------------------------------------------------------
vcftools --vcf hg19.hg002.vcf --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg19.hg002.indels
SURVIVOR filter hg19.hg002.indels.recode.vcf NA 50 -1 0 -1 hg19.hg002.indels.large.vcf



#-----------------------------------------------------------------------------------------------------------------------
# DeepVariant
# Get passed SNPs
vcftools --vcf SRR12898346.vcf --remove-filtered-all --remove-indels --recode --recode-INFO-all --bed $conbed --out SRR12898346.pass_SNPs
# Get passed indels
vcftools --vcf SRR12898346.vcf --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --bed $conbed --out SRR12898346.pass_indels

#-----------------------------------------------------------------------------------------------------------------------
# GATK
# SNPs
vcftools --gzvcf hg002.vcf.gz --remove-indels --recode --recode-INFO-all --bed $conbed --out hg002.SNPs
# Indels
vcftools --gzvcf hg002.vcf.gz --keep-only-indels --recode --recode-INFO-all --bed $conbed --out hg002.indels



########################################################################################################################
## Truvari
########################################################################################################################


########################################################################################################################
## vcfeval
########################################################################################################################

########################################################################################################################
## vcfdist
########################################################################################################################
conbed=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Evaluate SNPs
giabsnps=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.recode.vcf

#-----------------------------------------------------------------------------------------------------------------------
# Syri
vcfdist hg002.syri.SNPs.recode.vcf $giabsnps $ref -b $conbed -p vcfdist/
vcfdist hg002.nofilt.syri.SNPs.recode.vcf $giabsnps $ref -b $conbed -p vcfdist_nofilt/

#-----------------------------------------------------------------------------------------------------------------------
# Deepvariant
vcfdist SRR12898346.pass_SNPs.recode.vcf $giabsnps $ref -b $conbed -p vcfdist/

#-----------------------------------------------------------------------------------------------------------------------
# GATK
vcfdist hg002.SNPs.recode.vcf $giabsnps $ref -b $conbed -p vcfdist/ -v 0

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Evaluate indels
giabindels=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/HG002_GRCh38_1_22_v4.2.1_benchmark.indels.recode.vcf

#-----------------------------------------------------------------------------------------------------------------------
# Syri
vcfdist hg002.syri.indels.recode.vcf $giabindels $ref -b $conbed -p vcfdist_indels/
vcfdist hg002.nofilt.syri.indels.recode.vcf $giabindels $ref -b $conbed -p vcfdist_nofilt_indels/

#-----------------------------------------------------------------------------------------------------------------------
# Deepvariant
vcfdist SRR12898346.pass_indels.recode.vcf $giabindels $ref -b $conbed -p vcfdist_indels/

#-----------------------------------------------------------------------------------------------------------------------
# GATK
vcfdist hg002.indels.recode.vcf $giabindels $ref -b $conbed -p vcfdist_indels/ -v 0


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Evaluate SVs
giabsvs=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh37/HG002_SVs_Tier1_v0.6.passed.large.vcf
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa
conbed=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/SV_GRCh37/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed
#-----------------------------------------------------------------------------------------------------------------------
# Syri
vcfdist hg19.hg002.syri.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ -t 4 -v 0 &
vcfdist hg19.hg002.nofilt.syri.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_nofilt_svs/ -t 4 -v 0 &

# Paternal haplotype
vcfdist hg19.hg002.patsyri.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_pat_svs/
vcfdist hg19.hg002.patsyri.indels.large.dip_gt.vcf $giabsvs $ref -b $conbed -p vcfdist_pat_dipgt_svs/ &
vcfdist hg19.hg002.nofilt.patsyri.indels.large.dip_gt.vcf $giabsvs $ref -b $conbed -p vcfdist_nofilt_pat_dipgt_svs/ &

# Maternal haplotype
vcfdist hg19.hg002.matsyri.indels.large.dip_gt.vcf $giabsvs $ref -b $conbed -p vcfdist_mat_dipgt_svs/ -t 4 -v 0 &
vcfdist hg19.hg002.nofilt.matsyri.indels.large.dip_gt.vcf $giabsvs $ref -b $conbed -p vcfdist_nofilt_mat_dipgt_svs/ -t 4 -v 0 &

#-----------------------------------------------------------------------------------------------------------------------
# svim-asm
# merged
vcfdist hg19.hg002.dip.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ -t 4 -v 0 &
# Paternal haplotype
vcfdist hg19.hg002.pat.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_pat_svs/ -t 4 -v 0 &
# Maternal haplotype
vcfdist hg19.hg002.mat.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_mat_svs/ -t 4 -v 0 &

#-----------------------------------------------------------------------------------------------------------------------
# sniffles2
vcfdist hg19.hg002.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ &

#-----------------------------------------------------------------------------------------------------------------------
# cutesv
# HiFi reads
vcfdist hg19.hg002.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ -t 10 -v 0
# Diploid assembly
vcfdist hg19.hg002.diploid.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_diploid_svs/ -t 10 -v 0 

#-----------------------------------------------------------------------------------------------------------------------
# pbsv

vcfdist hg19.hg002.indels.large.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ -t 10 -v 0
#vcfdist hg19.hg002.vcf $giabsvs $ref -b $conbed -p vcfdist_svs/ -t 4 -v 0 &



