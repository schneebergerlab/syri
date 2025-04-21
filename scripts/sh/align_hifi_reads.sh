# Extract reads from BAM to a fastq file
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa.gz
bam=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/GRCh38.m84039_241001_220042_s2.hifi_reads.bc2018.bam
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/
samtools fastq -n --reference $ref $bam | bgzip -@40  > HG002.hifi.fq.gz

# Hg38 -----------------------------------------------------------------------------------------------------------------
# Align reads with PBMM2
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg38.filt.fa
hififq=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/HG002.hifi.fq.gz
pbmm2 align --sort -J 4 --preset HIFI --sample HG002 --rg "@RG\tID:HG002\tPL:Pacbio_hifi\tSM:1\tPU:1" -j 20 $ref $hififq HG002.pbmm2.sorted.bam
samtools index HG002.pbmm2.sorted.bam


# Hg19 -----------------------------------------------------------------------------------------------------------------
# Align reads with PBMM2
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments
ref=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa
hififq=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/HG002.hifi.fq.gz
pbmm2 align --sort -j 40 -J 4 --preset HIFI --sample HG002 --rg "@RG\tID:HG002\tPL:Pacbio_hifi\tSM:1\tPU:1" $ref $hififq hg19.HG002.pbmm2.sorted.bam
samtools index hg19.HG002.pbmm2.sorted.bam
