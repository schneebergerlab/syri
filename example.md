## Identifying genomic differences using SyRI

SyRI requires assemblies to be at chromosome-level for accurate identification of SRs. If chromosome-level assemblies are not available, one can create pseudo-chromosome level assemblies using [scaforder](scaforder.md) utility. 

### Whole-genome alignment
SyRI uses whole-genome alignments as input. These can be generated using whole-genome aligner of user's choice (check [installation guide](install.md). Here, we would use [MUMmer3](http://mummer.sourceforge.net/) package.

Firstly, the genomes (in multi-fasta format) are aligned using the NUCmer utility.
```bash
nucmer --maxmatch -c 500 -b 500 -l 100 refgenome qrygenome;
```

Here, `-c`,`-b`, and `-l` are parameters used to control the alignment resolution and need to be adjusted based on the genome size and complexity. Mere details are available [here](http://mummer.sourceforge.net/manual/#nucmer).

NUCmer would generate a `out.delta` file as output. The identified alignments are filtered using `delta-filter`  and then converted into a tab-separated [format](fileformat.md) using `show-coords`.

```bash
delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta; 
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords;
```

Users can change values for `-i`, and `-l` input to suite their genomes and specifc scientific problem. More information is available [here](http://mummer.sourceforge.net/manual/#filter).

For identificaiton of structural rearrangements (which include duplications), that overlapping alignments are not filtered out. In the example above, `--maxmatch` (for nucmer) results in identificaiton of all alignments. The `-m` (for delta-filter) parameter removes redundant alignments, though it is not necessary but is used as it helps in significantly reducing number of alignments which in turn reduces time and memory required by SyRI. Finally, `-THrd` (for show-coords) converts the alignments form `.delta` format to `.tsv` format consisting of alignment coordinates required by SyRI.

For alignments generated using MUMmer3, CIGAR strings are not required. For other, aligners it is necessary to have CIGAR string for identification of SNPs and short indels (structural rearrangements and structural variaition can be identified without CIGAR strings).

### SR identification using `syri`
SyRI takes genome alignments coordinates as input. Additionally, fasta files for the two genomes will also be required if structure variations are also needed. Further, for short variation identification, when CIGAR strings are not available, `.delta` file (as generated from NUCmer) will also be requried.

The usage and parameters are:

```
usage: syri [-h] -c INFILE [-r REF] [-q QRY] [-d DELTA] [-o FOUT] [-k]
            [--log {DEBUG,INFO,WARN}] [--lf LOG_FIN] [--dir DIR]
            [--prefix PREFIX] [--seed SEED] [--nc NCORES] [--novcf] [--nosr]
            [-b BRUTERUNTIME] [--unic TRANSUNICOUNT] [--unip TRANSUNIPERCENT]
            [--inc INCREASEBY] [--no-chrmatch] [--nosv] [--nosnp] [--all]
            [--allow-offset OFFSET] [--cigar] [-s SSPATH]

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
                        (SNPs/indels) identification when CIGAR string is not
                        available (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -o FOUT               Output file name (default: syri)
  -k                    Keep internediate output files (default: False)
  --log {DEBUG,INFO,WARN}
                        log level (default: INFO)
  --lf LOG_FIN          Name of log file (default: syri.log)
  --dir DIR             path to working directory (if not current directory)
                        (default: None)
  --prefix PREFIX       Prefix to add before the output file Names (default: )
  --seed SEED           seed for generating random numbers (default: 1)
  --nc NCORES           number of cores to use in parallel (max is number of
                        chromosomes) (default: 1)
  --novcf               Do not combine all files into one output file
                        (default: False)

SR identification:
  --nosr                Set to skip structural rearrangement identification
                        (default: False)
  -b BRUTERUNTIME       Cutoff to restrict brute force methods to take too
                        much time (in seconds). Smaller values would make
                        algorithm faster, but could have marginal effects on
                        accuracy. In general case, would not be required.
                        (default: 60)
  --unic TRANSUNICOUNT  Number of uniques bps for selecting translocation.
                        Smaller values would select smaller TLs better, but
                        may increase time and decrease accuracy. (default:
                        1000)
  --unip TRANSUNIPERCENT
                        Percent of unique region requried to select
                        translocation. Value should be in range (0,1]. Smaller
                        values would selection of translocation which are more
                        overlapped with other regions. (default: 0.5)
  --inc INCREASEBY      Minimum score increase required to add another
                        alignment to translocation cluster solution (default:
                        1000)
  --no-chrmatch         Do not allow SyRI to automatically match chromosome
                        ids between the two genomes if they are not equal
                        (default: False)

ShV identification:
  --nosv                Set to skip structural variation identification
                        (default: False)
  --nosnp               Set to skip SNP/Indel (within alignment)
                        identification (default: False)
  --all                 Use duplications too for variant identification
                        (default: False)
  --allow-offset OFFSET
                        BPs allowed to overlap (default: 0)
  --cigar               Find SNPs/indels using CIGAR string. Necessary for
                        alignment generated using aligners other than nucmers
                        (default: False)
  -s SSPATH             path to show-snps from mummer (default: show-snps)
```

#### SR identification
In case the chromosome IDs for the two assemblies are not identical, SyRI would try to find homologous chromosomes and then map their IDs to be identical. This behaviour can be turned off using the `--no-chrmatch` parameter.

Other parameters in this section regulate how translocation and duplications (TDs) are identified. For small networks of overlapping candidate TDs, SyRI uses a brute-force method to find the optimal set of TDs. The time allowed to this method can be restricted using the `-b` parameter. If for a network, brute-force method take more than the assigned time, then it will automatically switch to a randomized-greedy method. The `-unic` and `-unip` parameters state how unique a candidate TD need to be. Candidates which overlap highly with syntenic path and inversions and thus do not pass these thresholds will be filtered out. From a network of candidate TDs, it is possible to select different set of candidates. The `-inc` threshold is used decide whether a new set of candidates is better then the current candidate and thus can be selected as the solution or not. A new set will be considered as the better set if one of the following conditions satisfy: <br />
* <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)+inc\\" title="eq1" />
* <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)" title="eq2" /> <br /> and <br /> <img src="https://latex.codecogs.com/svg.latex?number\_of\_candidate(new\_set)\leq{number\_of\_candidate(current\_set)}" title="eq2" />
* <img src="https://latex.codecogs.com/svg.latex?score(new\_set)>score(current\_set)-inc" title="eq3" /> <br /> and <br /> <img src="https://latex.codecogs.com/svg.latex?number\_of\_candidate(new\_set)<number\_of\_candidate(current\_set)" title="eq3" />

#### Parameters for local variaition identificaiton
The `--allow-offset` parameter is used to define a threshold to decide whether two consecutive alignments within an annotated blocks are overlapping or not. Alignments, for which number of overlapping bases will be more than `--allow-offset` will result CNVs.
Short variations (SNPs, small indels) are identified by either using CIGAR strings for the alignments or using the `show-snps` utility in MUMmer package (would require `.delta` file). User can set `--cigar` when using CIGAR. If the `show-snps` is not in enironment path, then `-ss` can be used to provide path to it. By default, SyRI do not report short variations within duplicated regions because they lack one-to-one mapping between regions, which in turn renders short variations ambiguous. However, user can set `--all` which will return all short variations within all annotated alignments.
