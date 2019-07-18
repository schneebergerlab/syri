## Get Yeast Reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64/GCA_000146045.2_R64_genomic.fna.gz
gzip -d GCA_000146045.2_R64_genomic.fna.gz

## Get Query genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/977/955/GCA_000977955.2_Sc_YJM1447_v1/GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz
gzip -d GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz

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

## Perform whole genome alignment
ln -s GCA_000146045.2_R64_genomic.fna refgenome
ln -s GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.filtered qrygenome
nucmer --maxmatch -c 100 -b 500 -l 50 refgenome qrygenome       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta > out.delta              # Remove small and lower quality alignments
show-coords -THrd out.delta > out.coords                        # Convert alignment information to a .TSV format as required by SyRI

## Ryn SyRI 
python3 /biodata/dep_coupland/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c out.coords -d out.delta -r refgenome -q qrygenome

#OUTPUT: SyRI would report genomic structural differences in syri.out and syri.vcf.
