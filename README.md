## Synteny and Rearrangement Identifier (SyRI)
SyRI is a comprehensive tool for predicting genomic differences between related genomes using whole-genome assemblies (WGA). The assemblies are aligned using whole-genome alignment tools, and these alignments are then used as input to SyRI. SyRI identifies syntenic path (longest set of co-linear regions), structural rearrangements (inversions, translocations, and duplications), local variations (SNPs, indels, CNVs etc) within syntenic and structural rearrangements, and un-aligned regions.

SyRI uses an unprecedented approach where it starts by identifying longest syntenic path (set of co-linear regions). Since, all non-syntenic regions corresponds to genomic regions which have rearranged between the two genomes, identification of syntenic simultaneously identifies all structural rearrangements as well. After this step, all aligned non-syntenic regions are then classified as either inversion, translocation, or duplication based on the conformation of the constituting alignments. This approach transforms the challenging problem of SR identification to a comparatively easier problem of SR classificaiton.

Further, SyRI also identifies local variations within all syntenic and structurally rearranged regions. Local variations consists of short variations like SNPs, and small indels as well as structural variations like large indels, CNVs (copy-number variations), and HDRs. Short variations are parsed out from the constituting alignments, where as structural variations are predicting by comparing the overlaps and gaps between consecutive alignments of a syntenic or rearranged region.

# Contents
1. [Installation guide](install.md)
2. [Genomic difference identification](example.md)
3. [Pseudo-genome generation](scaforder.md)
4. [File-format](fileformat.md)
5. [Citation](#citation)


fads asdf
