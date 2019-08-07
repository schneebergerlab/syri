## An example pipeline for running SyRI
### Set variables
```
cwd="."     # Change to working directory
PATH_TO_SYRI="../syri/bin/syri" #Change the path to point to syri executable
```

### Got to the working directory and downloand the Yeast reference genome and a query assembly
```
cd $cwd
## Get Yeast Reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64/GCA_000146045.2_R64_genomic.fna.gz
## Get Query genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/977/955/GCA_000977955.2_Sc_YJM1447_v1/GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz

```

### Unzip the pre-process the genome
```
gzip -df GCA_000146045.2_R64_genomic.fna.gz
gzip -df GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.gz
## Remove mitochondrial DNA
head -151797 GCA_000977955.2_Sc_YJM1447_v1_genomic.fna > GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.filtered
```
#### Note
Ideally, syri expects that the homologous chromosomes in the two genomes would have exactly same chromosome id. Therefore, it is recommended that the user pre-processes the fasta files to ensure that homologous chromosomes have exactly the same id in both fasta files corresponding to the two genomes. In case, that is not the case, syri would try to find homologous genomes using whole genome alignments, but that method is heuristical and can result in suboptimal results. Also, it is recommended that the two genomes (fasta files) should have same number of chromosomes.

### Perform whole genome alignment
```
ln -sf GCA_000146045.2_R64_genomic.fna refgenome
ln -sf GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.filtered qrygenome
nucmer --maxmatch -c 100 -b 500 -l 50 refgenome qrygenome       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta              # Remove small and lower quality alignments
show-coords -THrd out.filtered.delta > out.filtered.coords                        # Convert alignment information to a .TSV format as required by SyRI
```

### Ryn SyRI 
```
python3 $PATH_TO_SYRI -c out.filtered.coords -d out.filtered.delta -r refgenome -q qrygenome
```
SyRI would report genomic structural differences in syri.out and syri.vcf.
