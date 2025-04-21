########################################################################################################################
## Pre-processing
########################################################################################################################
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa

#-----------------------------------------------------------------------------------------------------------------------
# GIAB Benchmark

cd /home/ra98jam/d16/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/

# Get SNPs only file
vcftools --gzvcf HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --remove-indels --recode --recode-INFO-all --out HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs

# Get indels only file
vcftools --gzvcf HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --keep-only-indels --recode --recode-INFO-all --out HG002_GRCh38_1_22_v4.2.1_benchmark.indels

# Consistent BED regions: Regions containing benchmark variants from GIAB
#conbed=/home/ra98jam/d16/projects/sv_benchmark/data/benchmark_vcf/hg002_giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.with_header.bed


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







