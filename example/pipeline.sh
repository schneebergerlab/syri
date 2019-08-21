cwd="."     # Change to working directory
cd $cwd
PATH_TO_SYRI="../syri/bin/syri" #Change the path to point to syri
## Get Yeast Reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64/GCA_000146045.2_R64_genomic.fna.gz
gzip -df GCA_000146045.2_R64_genomic.fna.gz

## Get Query genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/977/955/GCA_000977955.2_Sc_YJM1447_v1/GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz
gzip -df GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz

## Remove mitochondria
head -151797 GCA_000977955.2_Sc_YJM1447_v1_genomic.fna > GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.filtered

## NOTE: IDEALLY, SYRI EXPECTS THAT THE HOMOLOGOUS CHROMOSOMES IN THE TWO GENOMES
## WOULD HAVE EXACTLY SAME CHROMOSOME ID. THEREFORE, IT IS HIGHLY RECOMMENDED THAT
## THE USER PRE-PROCESSES THE FASTA FILES TO ENSURE THAT HOMOLOGOUS CHROMOSOMES
## HAVE EXACTLY THE SAME ID IN BOTH FASTA FILES CORRESPONDING TO THE TWO GENOMES.
## IN CASE, THAT IS NOT THE CASE, SYRI WOULD TRY TO FIND HOMOLOGOUS GENOMES USING
## WHOLE GENOME ALIGNMENTS, BUT THAT METHOD IS HEURISTICAL AND CAN RESULT IN 
## SUBOPTIMAL RESULTS.
## ALSO, IT IS REQUIRED THAT THE TWO GENOMES (FASTA FILES) SHOULD HAVE SAME NUMBER
## OF CHROMOSOMES.

ln -sf GCA_000146045.2_R64_genomic.fna refgenome
ln -sf GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.filtered qrygenome

## Perform whole genome alignment
# Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
minimap2 -ax asm5 --eqx refgenome qrygenome > out.sam

## Note: It is recommended that the user tests different alignment settings to find what
## alignment resolution suits their biological problem. Some alignment tools find
## longer alignments (with lots of gaps) while other find smaller more fragmented
## alignments. The smaller alignments generally have higher alignment identity scores
## and are more helpful in identifying smaller genomic structural rearrangments.
## But they could also lead to significant increase in redundant alignments which
## leads to increase in runtime of the alignment tool and SyRI.

## Ryn SyRI with SAM or BAM file as input
python3 $PATH_TO_SYRI -c out.sam -r refgenome -q qrygenome -k -F S
# OR
samtools view -b out.sam > out.bam
python3 $PATH_TO_SYRI -c out.bam -r refgenome -q qrygenome -k -F B

## SyRI would report genomic structural differences in syri.out and syri.vcf.


## Using SyRI to identify genomic rearrangements from whole-genome alignments 
## generated using MUMmer. A .tsv (out.filtered.coords) file is used as the input.
nucmer --maxmatch -c 100 -b 500 -l 50 refgenome qrygenome       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta     # Remove small and lower quality alignments
show-coords -THrd out.filtered.delta > out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI
python3 $PATH_TO_SYRI -c out.filtered.coords -d out.filtered.delta -r refgenome -q qrygenome

