## Requires:
- Python3 (preferred python3.5)
- Python packages:
    - Cython
    - Numpy
    - Scipy
    - Pandas
    - Igraph
    - Biopython
- MUMmer3

## Installation guide:
1. Download the directory
2. Open it and run:

```bash
$ python3 setup.py install
```

This will compile cython code.

All executables would be in ```syri/bin/```.

## Genome difference identification:

### Pre-requisite:
Chromosomal assemblies need to be aligned. This can be done using mummer's nucmer
tool. The nucmer specific parameter would depend on the complexity of the concerned
genomes (amount of repeat regions) and assembly quality (gaps).

```bash
## sample nucmer run, filtering of alignment (optional, but recommended), and transforming delta file into tab-delimited file format.
nucmer --maxmatch -c 500 -b 500 -l 100 refgenome qrygenome;
delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta; 
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords;
```

Here, `--maxmatch` (for nucmer), `-m` (for delta-filter), `-THrd` (for show-coords) are critical and should be used as such to ensure that the output format is correct.

### SR identification using `syri`:
This is the main method of this package. It takes `*.coords` file as input and process them to annotate structural rearrangements. 
syri can be run using the following command in working directory:
```bash
$ syri /path/to/coords/file [options]
```
<!---where, the accepted parameters are:
```
optional arguments:
  -h, --help          show this help message and exit
  -b BRUTERUNTIME     Cutoff to restrict brute force methods to take too much
                      time (in seconds). Smaller values would make algorithm
                      faster, but could have marginal effects on accuracy. In
                      general case, would not be required. (default: 60)
  -c TRANSUNICOUNT    Number of uniques bps for selecting translocation.
                      Smaller values would select smaller TLs better, but may
                      increase time and decrease accuracy. (default: 1000)
  -p TRANSUNIPERCENT  Percent of unique region requried to select
                      tranalocation. Value should be in range (0,1]. Smaller
                      values would selection of translocation which are more
                      overlapped with other regions. (default: 0.5)
  -nC NCORES          number of cores to use in parallel (max is number of
                      chromosomes) (default: 1)
  -d DIR              path to working directory (if not current directory)
                      (default: /biodata/dep_coupland/grp_schneeberger/project
                      s/SynSearch/scripts/python/syri/)
  -i INCREASEBY       Minimum score increase required to add another alingment
                      to translocation cluster solution (default: 1000)
  --prefix PREFIX     Prefix to add before the output file Names (default: )
  -s SEED             seed for generating random numbers (default: 1)

```-->

The output is stored in seven files corresponding to syntenic regions (synOut.txt) and six classes of SRs inversion (invOut.txt), translocation (TLOut.txt), inverted translocation (invTLOut.txt), duplication (dupOut.txt), inverted duplication (invDupOut.txt), and cross-chromosomal exchange (ctxOut.txt). The files use a two layer structure reporting annotated block and the alignments which constitute the block.

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
Here, the lines starting with '#' correspond to alignment block, and lines below it (till the beginning of next annotated block) are the alignments in this block. For blocks, the columns are: RefChr | RefStart | RefEnd |  - | QryChr |  QryStart |  QryEnd , and for alignments: RefStart | RefEnd | QryStart |  QryEnd.

### SV identification using `getsv`:
This tool uses the output of syri and identifies structure variations between the two genomes, outputs divergent (not aligned) regions, and creates input alignment file which can be used with show-snps (from mummer) to get SNPs and indels.

```bash
#Usage:
$ getsv [-d /path/to/directory/(if)/not/current/directory]
```
getsv can simply be run in the working directory (containing output files from syri) without any parameter, or 
### **scaffoldOrder.py**
- Script to generate pseudo-chromosomes using mummer alignment (coords file) and mummerplot output (.gp file)
- Run from the folder containing the coords file from mummer:
    - python \<path to scaffoldOrder.py\>  \<path to thequery genome\> \<coords file name\>
- Automatically scans the working directory for mummerplot output files. Each chromosome must have a separate mummerplot output
  - Run mummerplot as: mummerplot -l -f -r \<chrID\> -p \<chrID\> --postscript out.delta
- Test files are at /netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_Ler_inhouse

### **synsearchFunctions.py**
- collection of functions used by other files

### **myUsefulFunctions.py**
- collection of generally useful python functions

### **multiSV.py**
- Script to analyse multiple. Still work in progress. I think the code is very intuitive, so use can try it out on their own.

### **jobStarter.py**
- Script with example how to run process from within python. Not useful as of now






#### **Linux**:
- Install python, pip, and related libraries
	- ```apt install python3 python3-dev python3-pip python3-tk libxml2-dev zlib1g-dev```
- Install python modules
	- ```pip3 install setuptools numpy scipy pandas matplotlib python-igraph biopython glob2```
- Clone the directory
	- ```git clone https://www.github.com/schneebergerlab/syri```
- Install tool
	- ```cd syri```
	- ```python3 setup.py install```
- If pip fail to install python-igraph, then try ```easy_install python-igraph``` or first install igraph C core manually by ```apt install libigraph0-dev``` and then use pip to install python-igraph

#### **Windows**
- Download and install python3 from https://www.python.org/downloads/windows/
- Add python and pip executables to Environment Variables Path
- Download igraph binary from https://www.lfd.uci.edu/~gohlke/pythonlibs/ sutiable for your system
- Install the .whl from command line
	- Open cmd with admin right
	- Go to directory where .whl is saved
	- Run ```pip install <python_igraphXXXXX.whl>```
- Install remaining python modules
	- ```pip install setuptools numpy scipy pandas matplotlib biopython glob2```
- Download and extract syri https://www.github.com/schneebergerlab.syri
- From command line, go to the extracted directory and install module
	- ```python setup.py install```

General requirements:
  - python3+, corredponding pip (or easy_install, conda), python-dev
=======
General requirements: lol
  - python3+, corredponding pip (or easy_install, conda)
  - python dependencies: os, sys, pandas, numpy, igraph, collections, matplotlib, multiprocessing, functools, time, scipy, Biopython, subprocess, glob, operator, math, itertools
  - libxml2,libxml2-dev, zlib, python3-tk

#### **Installing dependencies**:
Linux:
  - Install igraph C libraries: ```apt install libigraph0-dev```	
  - Install dependencies: ```pip install setuptools numpy scipy pandas matplotlib python-igraph biopython glob2 datetime multiprocess```
  - Depending on your setup, may also need to install python-tk and python-scipy libraries: ```apt install python-tk python-scipy```

Windows:
  - Install igraph C libraries:
      - Download igraph binary suitable for your python from https://www.lfd.uci.edu/~gohlke/pythonlibs/#python-igraph
      - install using ```pip install <.whl>```
  - Install dependencies:
  	```pip install setuptools numpy scipy pandas matplotlib python-igraph biopython glob2 datetime multiprocess```


#### **Installing syri**:
  - Download and extract the package
  - Open folder
  - Run ```python setup.py install```



<!---

Table can be generated like this:
# **Syntax_block**
# 
# RefChr | RefStart | RefEnd |  - | QryChr |  QryStart |  QryEnd 
# --- | --- | --- | --- | --- | --- | ---
#  Chr1 |	8241 | 610363 |	- |	Chr1 |	1	 | 601274 
# 
# **Syntax_alignment**
# 
# | RefStart | RefEnd | QryStart |  QryEnd |
# |--- | --- | --- | ---|
# | 8241	| 550541 |	1 |	542302 |
-->
