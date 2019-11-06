## Pre-requisite for installing SyRI

### C/C++ compiler
g++
### Python
SyRI needs Python3.5 and depends on the following python libraries: [Cython](https://cython.org/#download), [Numpy](https://www.numpy.org/), [Scipy](https://www.scipy.org/install.html), [Pandas](https://pandas.pydata.org/), [Igraph](https://igraph.org/python/), [Biopython](https://biopython.org/), [psutil](https://github.com/giampaolo/psutil), [pysam](https://pysam.readthedocs.io/en/latest/index.html)

All these packages are available through anaconda and can be installed using:

```bash
conda install cython numpy scipy pandas biopython psutil
conda install -c conda-forge python-igraph
conda install -c bioconda pysam
```

### Whole genome aligner 
SyRI uses whole genome alignments as input. Users can use any aligner of their choice. SyRI accepts alignment input in SAM/BAM format or in a `tab-separated value` format with [CIGAR](https://samtools.github.io/hts-specs/SAMv1.pdf) string information for each alignment. If the user wants to use [MUMmer](http://mummer.sourceforge.net/), then ```.delta``` file can be used in place of CIGAR strings. See [fileformat](fileformat.md) for more information.

## Installation guide
1. Download (or clone) the repository for SyRI
2. Open the directory where ```setup.py``` is (we will call it as  ```cwd```) and in the terminal run:

```bash
python3 setup.py install		            # Install syri
chmod +x syri/bin/syri syri/bin/chroder		# Make it executable
```

All executables would be in ```cwd/syri/bin/```.
