## Synteny and Rearrangement Identifier (SyRI)
SyRI is a comprehensive tool for identifying all genomic differences between  genomes of related genomes using whole-genome assemblies. Whole-genome alignments between the assemblies are used to identify syntenic path (longest set of co-linear regions), structural rearrangements (inversions, translocations, and duplications), local variations (SNPs, indels, CNVs etc) and un-aligned regions.

SyRI is based on the observance that every region which is not structurally conserved between two genomes correspond to structural rearrangments. Consequently, SyRI starts by identifying all structurally conserved regions between the two genomes. These regions are called syntenic regions and correspond to the longest set co-linear regions between the two genomes. After this step, all non-syntenic regions are structural rearrangements (SRs) by definition and are then classified as either inversion, translocation, or duplication. This approach transforms the challenging problem of SR identification to a comparatively easier problem of SR classificaiton.

Further, SyRI also identifies local variations within all syntenic and structurally rearranged regions. Local variations consists of short variations like SNPs, and small indels as well as structural variations like large indels, CNVs (copy-number variations), and HDRs.

# Contents
1. [Pre-requisite](prereq.md)
2. [Installation guide](install.md)
3. [Genomic difference identification](example.md)
4. [Pseudo-genome generation](scaforder.md)
5. [File-format](fileformat.md)
6. [Citation](#citation)
