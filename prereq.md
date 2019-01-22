## Pre-requisite for installing SyRI

### Python
SyRI needs Python3 (recommended Python3.5) and depends on the following python libraries: [Cython](https://cython.org/#download), [Numpy](https://www.numpy.org/), [Scipy](https://www.scipy.org/install.html), [Pandas](https://pandas.pydata.org/), [Igraph](https://igraph.org/python/), [Biopython](https://biopython.org/). 

### Whole genome aligner 
To generate whole-genome alignment we used and recommend [MUMmer3](http://mummer.sourceforge.net/) package. User can use any whole-genome aligner of their choice for indetification of structural rearrangments and structural variations as long as the input file is in correct [file-format](fileformat.md). However, currently, SyRI can identify parse short variations (SNPs and small indels) from files generated using MUMmer3 only.
