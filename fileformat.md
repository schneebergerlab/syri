## File formats
### Input file format
SyRI takes whole genome alignment coordinates in a TSV format with the following columns:

|Column Number | Value   | Type |
| ------------ | ---------- | ----------- |
|1       | genome A start position (1-based, includes start position) | int        |
|2       | genome A end position (1-based, includes end position) | int |
|3       | genome B start position (1-based. Includes the start position.) | int |
|4       | genome B end position (1-based. Includes the end position.)   |  int |
|5       | alignment length on genome A  | int |
|6       | alignment length on genome B    |    int |
|7       | alingment identity (in percent, 0-100)  | float |
|8       | alignment direction in genome A  (always 1)  | int |
|9       | alignment direction in genome B (1 for directed alignments, -1 for inverted alignments)       | int |
|10      | chromosome ID in genome A     |         string |
|11      | chromosome ID in genome B         |         string |
|12      | CIGAR string corresponding to the alignment (Optional; '=' for match, 'X' for mismatch, 'D' for deletion, 'I' for insertion)    |         string |


Genomes are required to be provideed in multi-fasta format. Alternatively, nucmer generated `.delta` file can also be provided in place of CIGAR string for SNP identification.

### Output file format
SyRI outputs results in TSV format and VCF file format.

#### TSV format specifications
|Column Number | Value   | Type |
| ------------ | ---------- | ----------- |
|1       | chromosome ID in genome A     |         string |
|2       | genome A start position (1-based, includes start position) | int |
|3       | genome A end position (1-based, includes end position) | int |
|4       | Sequence in genome A (Only for SNPs and indels) | string |
|5       | Sequence in genome B (Only for SNPs and indels) | string |
|6       | chromosome ID in genome B     |         string |
|7       | genome B start position (1-based, includes start position) | int |
|8       | genome B end position (1-based, includes end position) | int |
|9       | Unique ID  (annotation type + number)  | string |
|10      | Parent ID  (annotation type + number)  | string |
|11      | Annotation type  | string |
|12      | Copy status (for duplications)| string |


#### VCF format
Above information is translated to VCF (v4.3) file format where genome A is considered as the reference genome. However, since VCF is based on reference genome position, we do not output un-aligned regions in genome B in VCF file, as it was not possible to write them in context of position in reference genome.


<!--
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
-->
