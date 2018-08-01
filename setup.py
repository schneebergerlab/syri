from distutils.core import setup
#from Cython.Build import cythonize

try:
    from Cython.Build import cythonize
except:
    print("Cython is not installed. Please install it. Exiting")
    sys.exit(1)
    
setup(name="synfunc",
        ext_modules=cythonize("syri/pyxFiles/synsearchFunctions.pyx"),
        packages=["syri","syri.bin", "syri.bin.func"])
