########################################################################################################################
## syri
########################################################################################################################
pat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852605.3_Q100_hg002v1.0.1.pat_genomic_autosomal_chr.fna.gz
mat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852615.3_Q100_hg002v1.0.1.mat_genomic_autosomal_chr.fna.gz

# Hg38 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.only_chromosome.fa.gz
## align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
{
minimap2 -ax asm5 -t 20 --eqx $ref $pat \
| samtools sort -@ 2 -O BAM - \
> hg002.pat.sorted.bam
samtools index hg002.pat.sorted.bam
minimap2 -ax asm5 -t 20 --eqx $ref $mat \
| samtools sort -@ 2 -O BAM - \
> hg002.mat.sorted.bam
samtools index hg002.mat.sorted.bam
} &

## Run syri
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/syri
syri -c ../../alignments/hg002.pat.sorted.bam -r $ref -q $pat -F B --nc 5 --prefix hg002.pat --samplename pat &
syri -c ../../alignments/hg002.mat.sorted.bam -r $ref -q $mat -F B --nc 5 --prefix hg002.mat --samplename mat &
syri -c ../../alignments/hg002.pat.sorted.bam -r $ref -q $pat -F B -f --nc 5 --prefix hg002.nofilt.pat --samplename pat &
syri -c ../../alignments/hg002.mat.sorted.bam -r $ref -q $mat -F B -f --nc 5 --prefix hg002.nofilt.mat --samplename mat &

# Hg19 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.only_chromosome.fa.gz
## align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
{
minimap2 -ax asm5 -t 20 --eqx $ref $pat \
| samtools sort -@ 2 -O BAM - \
> hg19.hg002.pat.sorted.bam
samtools index hg19.hg002.pat.sorted.bam
minimap2 -ax asm5 -t 20 --eqx $ref $mat \
| samtools sort -@ 2 -O BAM - \
> hg19.hg002.mat.sorted.bam
samtools index hg19.hg002.mat.sorted.bam
} &

## Run syri
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/syri
syri -c ../../alignments/hg19.hg002.pat.sorted.bam -r $ref -q $pat -F B --nc 5 --prefix hg19.hg002.pat --samplename pat &
syri -c ../../alignments/hg19.hg002.mat.sorted.bam -r $ref -q $mat -F B --nc 5 --prefix hg19.hg002.mat --samplename mat &
syri -c ../../alignments/hg19.hg002.pat.sorted.bam -r $ref -q $pat -F B -f --nc 5 --prefix hg19.hg002.nofilt.pat --samplename pat &
syri -c ../../alignments/hg19.hg002.mat.sorted.bam -r $ref -q $mat -F B -f --nc 5 --prefix hg19.hg002.nofilt.mat --samplename mat &


########################################################################################################################
## svim-asm
########################################################################################################################
pat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852605.3_Q100_hg002v1.0.1.pat_genomic_autosomal_chr.fna.gz
mat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852615.3_Q100_hg002v1.0.1.mat_genomic_autosomal_chr.fna.gz

# Hg38 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.only_chromosome.fa.gz
## Align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
{
minimap2 -a -x asm5 --cs -r2k -t 4 $ref $pat > hg002.pat.svim_asm.sam
minimap2 -a -x asm5 --cs -r2k -t 4 $ref $mat > hg002.mat.svim_asm.sam
samtools sort -m4G -@4 -o hg002.pat.svim_asm.bam hg002.pat.svim_asm.sam
samtools sort -m4G -@4 -o hg002.mat.svim_asm.bam hg002.mat.svim_asm.sam
samtools index hg002.pat.svim_asm.bam
samtools index hg002.mat.svim_asm.bam
} &
## Run svim-asm in haploid mode
svim-asm haploid hg002_pat ../../alignments/hg002.pat.svim_asm.bam $ref --sample pat &
svim-asm haploid hg002_mat ../../alignments/hg002.mat.svim_asm.bam $ref --sample mat &
## Run svim-asm in diploid mode
svim-asm diploid hg002_dip ../../alignments/hg002.pat.svim_asm.bam ../../alignments/hg002.mat.svim_asm.bam $ref --sample HG002 &

# Hg19 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.only_chromosome.fa.gz
## Align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
{
minimap2 -a -x asm5 --cs -r2k -t 4 $ref $pat > hg19.hg002.pat.svim_asm.sam
minimap2 -a -x asm5 --cs -r2k -t 4 $ref $mat > hg19.hg002.mat.svim_asm.sam
samtools sort -m4G -@4 -o hg19.hg002.pat.svim_asm.bam hg19.hg002.pat.svim_asm.sam
samtools sort -m4G -@4 -o hg19.hg002.mat.svim_asm.bam hg19.hg002.mat.svim_asm.sam
samtools index hg19.hg002.pat.svim_asm.bam
samtools index hg19.hg002.mat.svim_asm.bam
} &
## Run svim-asm in haploid mode
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/svim_asm
svim-asm haploid hg19.hg002_pat ../../alignments/hg19.hg002.pat.svim_asm.bam $ref --sample pat &
svim-asm haploid hg19.hg002_mat ../../alignments/hg19.hg002.mat.svim_asm.bam $ref --sample mat &
## Run svim-asm in diploid mode
svim-asm diploid hg19.hg002_dip ../../alignments/hg19.hg002.pat.svim_asm.bam ../../alignments/hg19.hg002.mat.svim_asm.bam $ref --sample HG002 &


########################################################################################################################
## Sniffles2
########################################################################################################################
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/sniffle2
# Hg38 -----------------------------------------------------------------------------------------------------------------
## Run sniffles
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/HG002.pbmm2.sorted.bam
sniffles --input $bam --vcf hg002.vcf.gz --reference $ref -t 10

# Hg19 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/hg19.HG002.pbmm2.sorted.bam
sniffles --input $bam --vcf hg19.hg002.vcf.gz --reference $ref -t 10


########################################################################################################################
## cutesv
########################################################################################################################
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/cutesv

# Hg38 -----------------------------------------------------------------------------------------------------------------
## Run cuteSV for hifi reads
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/HG002.pbmm2.sorted.bam
cuteSV $bam $ref hg002.vcf.gz hg002 --retain_work_dir -t 20 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -S HG002 --genotype &

## Run cuteSV for diploid assembly calling (https://github.com/tjiangHIT/cuteSV/wiki/Diploid-assembly-based-SV-detection-using-cuteSV)
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.only_chromosome.fa.gz
pat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852605.3_Q100_hg002v1.0.1.pat_genomic_autosomal_chr.fna.gz
mat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852615.3_Q100_hg002v1.0.1.mat_genomic_autosomal_chr.fna.gz
### Align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/cutesv
cat <(cat $pat ) <(cat $mat) > hg002.fa.gz
minimap2 --paf-no-hit -a -x asm5 --cs -r 2k -t 8 $ref hg002.fa.gz \
| samtools sort -@ 8 -o hg002.dip_asm.sorted.bam - &
samtools index hg002.dip_asm.sorted.bam
cuteSV hg002.dip_asm.sorted.bam $ref hg002.dip_asm.vcf hg002_dip_asm \
    -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5
python  /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/bmtools/lib/python3.8/site-packages/cuteSV/diploid_calling.py hg002_dip_asm/hg002.dip_asm.vcf hg002.dip_asm.merged.vcf

# Hg19 -----------------------------------------------------------------------------------------------------------------

## Run cuteSV for hifi reads
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/hg19.HG002.pbmm2.sorted.bam
cuteSV $bam $ref hg19.hg002.vcf.gz hg19.hg002 --retain_work_dir -t 20 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -S HG002 --genotype &

## Run cuteSV for diploid assembly calling (https://github.com/tjiangHIT/cuteSV/wiki/Diploid-assembly-based-SV-detection-using-cuteSV)
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.only_chromosome.fa.gz
pat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852605.3_Q100_hg002v1.0.1.pat_genomic_autosomal_chr.fna.gz
mat=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GCA_018852615.3_Q100_hg002v1.0.1.mat_genomic_autosomal_chr.fna.gz
### Align genomes
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/cutesv
cat <(cat $pat ) <(cat $mat) > hg002.fa.gz
minimap2 --paf-no-hit -a -x asm5 --cs -r 2k -t 20 $ref hg002.fa.gz \
| samtools sort -@ 20 -o hg19.hg002.dip_asm.sorted.bam - &
samtools index hg19.hg002.dip_asm.sorted.bam
cuteSV hg19.hg002.dip_asm.sorted.bam $ref hg19.hg002.dip_asm.vcf hg19.hg002_dip_asm \
    -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5
python  /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/bmtools/lib/python3.8/site-packages/cuteSV/diploid_calling.py hg19.hg002.dip_asm.vcf hg19.hg002.dip_asm.merged.vcf

# 28.04.2025 Running alignments using the new documentation https://github.com/tjiangHIT/cuteSV/wiki/Diploid-assembly-based-SV-detection-using-cuteSV

cat <(pigz -d -c -p 10 $pat | sed 's/>/>cutesvh1_/g') \
    <(pigz -d -c -p 10 $mat | sed 's/>/>cutesvh2_/g') \
    | bgzip -@10  > hg002.v2.fa.gz

minimap2 --paf-no-hit -a -x asm5 --cs -r 2k -t 25 $ref hg002.v2.fa.gz \
|  samtools sort -@ 8 -o hg002.v2.dip_asm.sorted.bam -
samtools index  hg002.v2.dip_asm.sorted.bam
cuteSV hg002.v2.dip_asm.sorted.bam $ref hg19.hg002.v2.dip_asm.vcf ./ \
    -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5
python3  /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/bmtools/lib/python3.8/site-packages/cuteSV/diploid_calling.py hg19.hg002.v2.dip_asm.vcf  hg19.hg002.v2.dip_asm.merged.vcf

########################################################################################################################
## PBSV
########################################################################################################################
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/pbsv

# Hg38 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/HG002.pbmm2.sorted.bam

# Discover signatures of structural variation
pbsv discover --hifi -s HG002 $bam HG002.svsig.gz

# index svsig.gz to allow random access via `pbsv call -r`
tabix -c '#' -s 3 -b 4 -e 4 HG002.svsig.gz

# Call SVs
pbsv call --hifi -j 10 $ref HG002.svsig.gz hg002.vcf

# Hg19 -----------------------------------------------------------------------------------------------------------------
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/hg19.HG002.pbmm2.sorted.bam

# Discover signatures of structural variation
pbsv discover --hifi -s HG002 $bam hg19.HG002.svsig.gz

# index svsig.gz to allow random access via `pbsv call -r`
tabix -c '#' -s 3 -b 4 -e 4 hg19.HG002.svsig.gz

# Call SVs
pbsv call --hifi -j 10 $ref hg19.HG002.svsig.gz hg19.hg002.vcf

########################################################################################################################
## deepvariant
########################################################################################################################
# Align reads
/home/ra98jam/d16/projects/sv_benchmark/scripts/sh/SBATCH/SBATCH_align_SRR12898346.sh
# Call variants on AI cluster
/home/ra98jam/d16/projects/sv_benchmark/scripts/sh/SBATCH/SBATCH_run_deepvariant.sh

########################################################################################################################
## GATK
########################################################################################################################
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/SRR12898346.sorted.dupmark.bam
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/gatk
gatk HaplotypeCaller \
-I $bam \
-O hg002.vcf.gz \
-R $ref \
-A VariantType

## Filter variant calls or not ??
