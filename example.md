### Running SyRI for identificaion of genomic differences
## Genome difference identification:

### Pre-requisite:
Chromosomal assemblies need to be aligned. This can be done using mummer's nucmer
tool. The nucmer specific parameter would depend on the complexity of the concerned
genomes (amount of repeat regions) and assembly quality (gaps).

```bash
## sample nucmer run, filtering of alignment (optional, but recommended), and transforming delta file into tab-delimited file format.
nucmer --maxmatch -c 500 -b 500 -l 100 refgenome qrygenome;
delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta; 
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords;
```

Here, `--maxmatch` (for nucmer), `-m` (for delta-filter), `-THrd` (for show-coords) are critical and should be used as such to ensure that the output format is correct.

### SR identification using `syri`:
This is the main method of this package. It takes `*.coords` file as input and process them to annotate structural rearrangements. 
syri can be run using the following command in working directory:
```bash
syri /path/to/coords/file [options]
```
<!---where, the accepted parameters are:
```
optional arguments:
  -h, --help          show this help message and exit
  -b BRUTERUNTIME     Cutoff to restrict brute force methods to take too much
                      time (in seconds). Smaller values would make algorithm
                      faster, but could have marginal effects on accuracy. In
                      general case, would not be required. (default: 60)
  -c TRANSUNICOUNT    Number of uniques bps for selecting translocation.
                      Smaller values would select smaller TLs better, but may
                      increase time and decrease accuracy. (default: 1000)
  -p TRANSUNIPERCENT  Percent of unique region requried to select
                      tranalocation. Value should be in range (0,1]. Smaller
                      values would selection of translocation which are more
                      overlapped with other regions. (default: 0.5)
  -nC NCORES          number of cores to use in parallel (max is number of
                      chromosomes) (default: 1)
  -d DIR              path to working directory (if not current directory)
                      (default: /biodata/dep_coupland/grp_schneeberger/project
                      s/SynSearch/scripts/python/syri/)
  -i INCREASEBY       Minimum score increase required to add another alingment
                      to translocation cluster solution (default: 1000)
  --prefix PREFIX     Prefix to add before the output file Names (default: )
  -s SEED             seed for generating random numbers (default: 1)

```-->

The output is stored in seven files corresponding to syntenic regions (synOut.txt) and six classes of SRs inversion (invOut.txt), translocation (TLOut.txt), inverted translocation (invTLOut.txt), duplication (dupOut.txt), inverted duplication (invDupOut.txt), and cross-chromosomal exchange (ctxOut.txt). The files use a two layer structure reporting annotated block and the alignments which constitute the block.

```
#	Chr1	8241	610363	-	Chr1	1	601274              
8241	550541	1	542302  
549844	587482	541241	578850  
588093	610363	578947	601274  
#	Chr1	610355	1160239	-	Chr1	602856	1153904  
610355	670785	602856	663228  
671022	768174	663228	760407  
768166	821005	761285	814172  
```
Here, the lines starting with '#' correspond to alignment block, and lines below it (till the beginning of next annotated block) are the alignments in this block. For blocks, the columns are:

- RefChr
- RefStart
- RefEnd
- "-"
- QryChr
- QryStart
- QryEnd

For alignments, the columns are:

- RefStart
- RefEnd
- QryStart
- QryEnd.

### SV identification using `getsv`:
This tool uses the output of syri and identifies structure variations between the two genomes, outputs divergent (not aligned) regions, and transform syri's output into a list of alignment.

```bash
#Usage:
getsv [-d /path/to/directory/if/not/current/directory]
```
It generates three output files.

notAligned.txt: Lists all divergenet (not aligned) regions in the two genomes. It is in .tsv format with the columns being:

- Genome identifier. R for reference and Q for query
- Start position
- End position
- Chromosome ID

mumSNPIn.txt: list of alignments that were aligned. Input for getShV. Can directly be used with show-snps (from mummer). In tsv format, each row is an alignment, and columns being:

- Start position in reference
- End position in reference
- Start position in query
- End position in query
- Chromosome ID in reference
- Chromosome ID in query

sv.txt: Structural variations in the two genomes. Reports, CNVs (CNV), highly different regions (HDR), indels (InDel), tandem repeat (CNV+Tandem), CNV with an indel between the two repeat sites, and Indel which are accompanied by a SNP. TSV format where each row is a sv, with columns being:

- Variation type
- Start position in reference
- End position in reference
- Start position in query
- End position in query
- Chromosome ID in reference
- Chromosome ID in query
- Secondary variation. Format-> type:genome containing the variation:position

```bash
# sample output
InDel+SNP   3900260 3901020 3909459 3909459 Chr1    Chr1    SNP:Q:3909459-3909459

## The region Ref-Chr1:3900260-3901020 should be at Qry-Chr1:3909459-3909459 but is deleted in the query genome. However, there is also SNP at Qry-Chr1:3909459.
```

### Short variation identification using `getshv`:
SyRI allows identification of short variations (SNPs and indels) along with information of their actual biological confirmation within the compared genomes. Annotated alignment are compared and processed to find ShVs. Alignment file generated by getsv, and the alignment delta file are parsed to show-snps (from mummer). Its output is then further processed to classify ShVs to incorporate confirmation information.

To use `getshv`:
```bash
getshv alignment_file delta_file [options]
```

Output files:

- snps.txt: show-snps output
- snps_no_indels.txt: show-snps output without indels
- snps_no_indels_buff0.txt: -b filtered snps
- snps_no_indels_buff0_syn: -b filtered snps in syntenic regions
- snps_no_indels_buff0_inv: -b filtered snps in inverted regions
- snps_no_indels_buff0_TL: -b filtered snps in translocated regions
- snps_no_indels_buff0_invTL: -b filtered snps in inverted translocated regions
- snps_no_indels_buff0_dup: -b filtered snps in duplicated regions
- snps_no_indels_buff0_invDup: -b filtered snps in inverted duplicated regions
- snps_no_indels_buff0_ctx: -b filtered snps in cross-chromosomal exchange
- snps_no_indels_buff0_strict_syn: -b filtered which are syntenic regions which are not overlapping with duplications

similarly, indels are divided into files corresponding to different regions (indels_syn, indels_inv, indels_TL, indels_invTL, indels_dup, indels_invDup, indels_ctx)

Outfile file format is same as that from show-snps with parameters -H, -l, -r, -T.
