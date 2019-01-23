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

Other parameters in this section regulate how translocation and duplications (TDs) are identified. For small networks of overlapping candidate TDs, SyRI uses a brute-force method to find the optimal set of TDs. The time allowed to this method can be restricted using the `-b` parameter. If for a network, brute-force method take more than the assigned time, then it will automatically switch to a randomized-greedy method. The `-unic` and `-unip` parameters state how unique a candidate TD need to be. Candidates which overlap highly with syntenic path and inversions and thus do not pass these thresholds will be filtered out. From a network of candidate TDs, it is possible to select different set of candidates. The `-inc` threshold is used decide whether a new set of candidates is better then the current candidate and thus can be selected as the solution or not. A new set will be considered as the better set if one of the following conditions satisfy: <br />
1. <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)+inc\\" title="eq1" /> <br />
2. <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)" title="eq2" /> <br /> and <br /> <img src="https://latex.codecogs.com/svg.latex?number\_of\_candidate(new\_set)\leq{}number\_of\_candidate(current\_set)" title="eq2" /> <br />
3. <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)-inc" title="eq3" /> <br /> and <br /> <img src="https://latex.codecogs.com/svg.latex?number\_of\_candidate(new\_set)<number\_of\_candidate(current\_set)" title="eq3" />

#### Parameters for local variaition identificaiton
The `--allow-offset` parameter is used to define a threshold to decide whether two consecutive alignments within an annotated blocks are overlapping or not. Alignments, for which number of overlapping bases will be more than `--allow-offset` will result CNVs.
For short variation identification (SNPs, small indels), SyRI uses a wrapper around `show-snps` utility in MUMmer package which parses out short variaitons from annotated alignments. Consequently, short variations can only be identified using alignments generated by MUMmer. If the `show-snps` is not in enironment path, then `-ss` can be used to provide path to it. By default, SyRI do not report short variations within duplicated regions because there the lack one-to-one mapping between regions, which in turn renders short variations ambiguous. However, user can set `--all` which will return all short variations within all annotated alignments. Further, `-buff` parameter can be used to remove low-confidence (which have other variations in its neighbouring regions) SNPs.
