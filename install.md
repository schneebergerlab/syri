## Pre-requisite for installing SyRI
### Python
SyRI needs Python3 (recommended Python3.5 or above) and depends on the following python libraries: [Cython](https://cython.org/#download), [Numpy](https://www.numpy.org/), [Scipy](https://www.scipy.org/install.html), [Pandas](https://pandas.pydata.org/), [Igraph](https://igraph.org/python/), [Biopython](https://biopython.org/), [psutil](https://github.com/giampaolo/psutil)

All these packages are available through and be installed using:

```bash
conda install cython numpy scipy pandas biopython psutil
conda install -c conda-forge python-igraph
```

### Whole genome aligner 
SyRI uses whole genome alignments as input. Users can use any aligner of their choice as long as the alignments can be parsed in .tsv format and [CIGAR](https://samtools.github.io/hts-specs/SAMv1.pdf) string information is present for each alignment. If the user wants to use [MUMmer3](http://mummer.sourceforge.net/), then ```.delta``` file can be used in place of CIGAR strings. Also, see [fileformat](fileformat.md) for more information.

## Installation guide
1. Download (or clone) the repository for SyRI
2. Open the directory where ```setup.py``` is (we will call it as  ```cwd```) and in the terminal run:

```bash
python3 setup.py install
```

All executables would be in ```cwd/syri/bin/```.
