## Identifying genomic differences using SyRI

SyRI requires assemblies to be at chromosome-level for accurate identification of SRs. If chromosome-level assemblies are not available, one can create pseudo-chromosome level assemblies using [scaforder](scaforder.md) utility. 

### Whole-genome alignment
SyRI uses whole-genome alignments as input. These can be generated using the [MUMmer3](http://mummer.sourceforge.net/) package. Firstly, the genomes (in multi-fasta format) are aligned using the NUCmer utility.
```bash
nucmer --maxmatch -c 500 -b 500 -l 100 refgenome qrygenome;
```

Here, `-c`,`-b`, and `-l` are parameters used to control the alignment resolution and need to be twerked based on the genome size and complexity. Mere details are available [here](http://mummer.sourceforge.net/manual/#nucmer).

NUCmer would generate a `out.delta` file as output. The identified alignments are filtered using and then converted into a tab-separated [format](fileformat.md) as required by SyRI.

```bash
delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta; 
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords;
```

Users can change values for `-i`, and `-l` input to suite their genomes and problem. More information is available [here](http://mummer.sourceforge.net/manual/#filter).

Here, `--maxmatch` (for nucmer), `-m` (for delta-filter), `-THrd` (for show-coords) are essential and should be used as such.

### SR identification using `syri`:
This is the main method of this package. It takes genome alignments coordinates as input in a tsv format. Additionally, fasta files for the two genomes will also be required if structure variations are also needed. Further, for short variation identification `delta` file (as generated from NUCmer) will also be requried.

The usage and parameters are:

```
usage: syri [-h] -c INFILE [-r REF] [-q QRY] [-d DELTA]
            [-log {DEBUG,INFO,WARN}] [-lf LOG_FIN] [-dir DIR]
            [--prefix PREFIX] [-seed SEED] [-nc NCORES] [-k] [-o FOUT]
            [-novcf] [-nosr] [-b BRUTERUNTIME] [-unic TRANSUNICOUNT]
            [-unip TRANSUNIPERCENT] [-inc INCREASEBY] [--no-chrmatch] [-nosv]
            [-nosnp] [--all] [--allow-offset OFFSET] [-ss SSPATH] [-buff BUFF]

Input Files:
  -c INFILE             File containing alignment coordinates in a tsv format
                        (default: None)
  -r REF                Genome A (which is considered as reference for the
                        alignments). Required for local variation (large
                        indels, CNVs) identification. (default: None)
  -q QRY                Genome B (which is considered as query for the
                        alignments). Required for local variation (large
                        indels, CNVs) identification. (default: None)
  -d DELTA              .delta file from mummer. Required for short variation
                        (SNPs/indels) identification (default: None)


optional arguments:
  -h, --help            show this help message and exit
  -log {DEBUG,INFO,WARN}
                        log level (default: INFO)
  -lf LOG_FIN           Name of log file (default: syri.log)
  -dir DIR              path to working directory (if not current directory)
                        (default: None)
  --prefix PREFIX       Prefix to add before the output file Names (default: )
  -seed SEED            seed for generating random numbers (default: 1)
  -nc NCORES            number of cores to use in parallel (max is number of
                        chromosomes) (default: 1)
  -k                    Keep internediate output files (default: False)
  -o FOUT               Output file name (default: syri)
  -novcf                Do not combine all files into one output file
                        (default: False)

SR identification:
  -nosr                 Set to skip structural rearrangement identification
                        (default: False)
  -b BRUTERUNTIME       Cutoff to restrict brute force methods to take too
                        much time (in seconds). Smaller values would make
                        algorithm faster, but could have marginal effects on
                        accuracy. In general case, would not be required.
                        (default: 60)
  -unic TRANSUNICOUNT   Number of uniques bps for selecting translocation.
                        Smaller values would select smaller TLs better, but
                        may increase time and decrease accuracy. (default:
                        1000)
  -unip TRANSUNIPERCENT
                        Percent of unique region requried to select
                        translocation. Value should be in range (0,1]. Smaller
                        values would selection of translocation which are more
                        overlapped with other regions. (default: 0.5)
  -inc INCREASEBY       Minimum score increase required to add another
                        alignment to translocation cluster solution (default:
                        1000)
  --no-chrmatch         Do not allow SyRI to automatically match chromosome
                        ids between the two genomes if they are not equal
                        (default: False)

ShV identification:
  -nosv                 Set to skip structural variation identification
                        (default: False)
  -nosnp                Set to skip SNP/Indel (within alignment)
                        identification (default: False)
  --all                 Use duplications too for variant identification
                        (default: False)
  --allow-offset OFFSET
                        BPs allowed to overlap (default: 0)
  -ss SSPATH            path to show-snps from mummer (default: show-snps)
  -buff BUFF            Remove SNPs which have other variants or alignment
                        break within buff size bps (default: 0)
```
#### Parameters related to SR identification
In case the chromosome IDs for the two assemblies are not identical, SyRI would try to find homologous chromosomes and then map their IDs to be identical. This behaviour can be turned off using the `--no-chrmatch` parameter.

Other parameters in this section regulate how translocation and duplications (TDs) are identified. For small networks of overlapping candidate TDs, SyRI uses a brute-force method to find the optimal set of TDs. The time allowed to this method can be restricted using the `-b` parameter. If for a network, brute-force method take more than the assigned time, then it will automatically switch to a randomized-greedy method. The `-unic` and `-unip` parameters state how unique a candidate TD need to be. Candidates which overlap highly with syntenic path and inversions and thus do not pass these thresholds will be filtered out. From a network of candidate TDs, it is possible to select different set of candidates. The `-inc` threshold is used decide whether a new set of candidates is better then the current candidate and thus can be selected as the solution or not. A new set will be considered as the better set if:

```

```
```math
a + v = 2
```

$a+b=2$
    h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x

The output is stored in seven files corresponding to syntenic regions (synOut.txt) and six classes of SRs inversion (invOut.txt), translocation (TLOut.txt), inverted translocation (invTLOut.txt), duplication (dupOut.txt), inverted duplication (invDupOut.txt), and cross-chromosomal exchange (ctxOut.txt). The files use a two layer structure reporting annotated block and the alignments which constitute the block.


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
