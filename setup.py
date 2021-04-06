from setuptools import setup, Extension
import sys
import numpy

try:
    from Cython.Build import cythonize
except ImportError:
    print("Cython is not installed. Please install it. Exiting")
    sys.exit(1)


    
setup(name="syri",
        version='1.4.1',
        ext_modules=cythonize([Extension('syri.pyxFiles.synsearchFunctions',['syri/pyxFiles/synsearchFunctions.pyx']),
        Extension('syri.inversions',['syri/pyxFiles/inversions.pyx']),
        Extension('syri.tdfunc',['syri/pyxFiles/tdfunc.pyx']),
        Extension('syri.findshv',['syri/pyxFiles/findshv.pyx']),
        Extension('syri.findsv',['syri/pyxFiles/findsv.pyx']),
        Extension('syri.writeout',['syri/pyxFiles/writeout.pyx'])]),
        packages=["syri","syri.bin", "syri.bin.func"],
        include_dirs=[numpy.get_include()])
