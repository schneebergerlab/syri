from distutils.core import setup
from Cython.Build import cythonize

setup(name="synfunc",
        ext_modules=cythonize("pyxFiles/synsearchFunctions.pyx"))
