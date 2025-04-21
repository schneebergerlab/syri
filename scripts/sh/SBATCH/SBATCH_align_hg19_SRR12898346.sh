#!/bin/bash
###SBATCH --array=0-9
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=2-00:00:00
#SBATCH -J align_SRR12898346_hg19

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/
idx=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/hg19.filt.fa
R1=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/SRR12898346_1.fastq.gz
R2=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/reads_assemblies/hg002/SRR12898346_2.fastq.gz
samplename="SRR12898346"


/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/bowtie2 --local \
  --very-sensitive-local \
  --threads ${SLURM_CPUS_PER_TASK} \
  -x $idx \
  -1 $R1 \
  -2 $R2 \
  --rg "ID:$samplename" \
  --rg "PL:ILLUMINA" \
  --rg "SM:1" \
  --rg "PU:1" \
| /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/samtools sort -@ ${SLURM_CPUS_PER_TASK} -O BAM - \
> hg19.SRR12898346.sorted.bam
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/samtools index hg19.SRR12898346.sorted.bam


picard MarkDuplicates \
-I hg19.SRR12898346.sorted.bam \
-M hg19.mark_dup_stats.txt \
-O hg19.SRR12898346.sorted.dupmark.bam
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/samtools index hg19.SRR12898346.sorted.dupmark.bam
