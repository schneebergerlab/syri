## Pre-requisite:
1. Python3.5 and the following packages: Cython, numpy, scipy, pandas (version: 0.23.4), python-igraph, biopython, psutil, pysam, and matplotlib
2. C/C++ compiler: g++

## Recent major updates:
* Added a new method `plotsr` in `./syri/bin/` for plotting genomic structures (syntenic and structurally rearranged regions)
* Optimised inversion identification improving worst-case performance by many folds.
 

## Installation:
Download/clone the repository, open the folder and run:

`python setup.py install`

This will install the cython modules.

## Running:
Executables will be in `syri/bin/` folder and can be run directly from the terminal.

Detailed information is available at: https://schneebergerlab.github.io/syri

## Citation:
Please cite:

`Goel, M., Sun, H., Jiao, W. et al. SyRI: finding genomic rearrangements and local sequence differences from whole-genome assemblies. Genome Biol 20, 277 (2019) doi:10.1186/s13059-019-1911-0`

## Current Limitations:
1. The homologous chromosomes in the two genomes need to represent the same strand. If the chromosomes are from different strands, then the alignments between these chromosomes would be inverted. As SyRI only checks directed alignments for syntenic region identification, it would not be able to find syntenic regions and can crash.  
Current solution to this problem is to manually check alignments. If the majority of alignments between homologous chromosomes are inverted, then the chromosome in the query genome needs to be reverse-complemented. Then the corrected query genome needs to be aligned against the reference genome. We are working on a method which can generate dot plots to automatically identify and reverse-complement such inverted-chromosomes.

2. Large translocations and duplications (consisting of multiple alignments) can result in high memory-usage and CPU runtime.
