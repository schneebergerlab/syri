<<<<<<< HEAD
General requirements:
  - python3+, corredponding pip (or easy_install, conda), python-dev
=======
General requirements: lol
  - python3+, corredponding pip (or easy_install, conda)
>>>>>>> 7fc4f6347ffb4347eaffcde6aee5c4f5993d94bd
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

