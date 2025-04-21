#!/bin/bash
#SBATCH -p lrz-dgx-a100-80x8
#SBATCH --gres=gpu:2
#SBATCH -o deepvar.out
#SBATCH -e deepvar.err
###SBATCH --container-image="docker://nvcr.io/nvidia/clara/clara-parabricks-deepvariant:4.0.0-1"
###SBATCH --container-mounts="./tmpdata:/mnt/"
#SBATCH -J deepvar_test



#reffile=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#bamfile=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/tmp.bam
#varfile=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments/tmp.vcf

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/deepvariant/
#
#reffile=GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#bamfile=tmp.bam
#varfile=tmp.vcf
#
# Check that the enroot container is present on the GPU node
#if  enroot list | grep -q 'deepvariant' ; then
if  enroot list | grep -q 'parabricks' ; then
  :
else
#  enroot create --name deepvariant ~/nvidia+clara+clara-parabricks-deepvariant+4.0.0-1.sqsh
  enroot create --name parabricks ~/nvidia+clara+clara-parabricks+4.4.0-1.sqsh
fi
#cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/tests/chr22.fa .
#cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/tests/tmp.sorted.bam .
#

#     python script.py --epochs 55 --batch-size 512
#time enroot start --mount ./:/mnt/ parabricks pbrun deepvariant  \
time enroot start \
    --mount /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/data/:/mnt/genome \
    --mount /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/alignments:/mnt/alignment \
    --mount /dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/deepvariant:/mnt/outdir \
    parabricks \
    pbrun deepvariant  \
    --ref /mnt/genome/hg38.filt.fa \
    --in-bam /mnt/alignment/SRR12898346.sorted.dupmark.bam \
    --out-variants /mnt/outdir/SRR12898346.vcf \
    --num-gpus 2
#
#time enroot start --mount ./:/mnt/ parabricks pbrun deepvariant  \
#--ref /mnt/chr22.fa \
#--in-bam /mnt/tmp.sorted.bam \
#--out-variants /mnt/tmp.vcf \
#--num-gpus 1