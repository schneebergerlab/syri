from distutils.core import setup
import sys
import numpy

try:
    from Cython.Build import cythonize
except ImportError:
    print("Cython is not installed. Please install it. Exiting")
    sys.exit(1)


    
setup(name="synfunc",
        ext_modules=cythonize("syri/pyxFiles/synsearchFunctions.pyx"),
        packages=["syri","syri.bin", "syri.bin.func"],
        include_dirs=[numpy.get_include()])
