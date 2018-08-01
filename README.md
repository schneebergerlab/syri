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

## **Installation guide**:
1. Download the directory
2. Open it and run:

```bash
python3 setup.py install
```

This will compile cython code.

## **Genome difference identification**:

### Pre-requisite:
Chromosomal assemblies need to be aligned. This can be done using mummer's nucmer
tool. The nucmer specific parameter would depend on the complexity of the concerned
genomes (amount of repeat regions) and assembly quality (gaps).

```bash
## sample nucmer run
nucmer --maxmatch -c 500 -b 500 -l 100 refgenome qrygenome;
delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta; 
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords;
```

Here, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as suchHere, --maxmatch (for nucmer), -m (for delta-filter), -THrd(for show-coords) are critical and should be used as such
## nucmer

### SR identification:







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



  
### **SynSearch.py**:
- Main script to find structural rearrangement between two genomes.
- run from the folder containing the coords file from mummer:
    - python \<path to SynSearch.py\> \<coords file name\> \<number of CPU cores to use\>
- Each chromosomes runs on a separate CPU core, so using 1 CPU will analyse chromosomes sequentially, while using 5 cores (for A.thaliana) will analyse all chromosomes in parallel
- Output will be created in same directory (6 files per chromosome + 1 file for CTX)
- Sample input files are in /netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/ directory.
  - Ex: /netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_ler_Chr/out_m_i90_l100.coords

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

