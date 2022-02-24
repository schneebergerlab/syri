from setuptools import setup, Extension
import numpy
from Cython.Build import cythonize
from syri import __version__

setup(name="syri",
      version='{}'.format(__version__),
      description='Synteny and rearrangement identifier between whole-genome assemblies',
      author='Manish Goel',
      author_email='goel@mpipz.mpg.de',
      url='https://github.com/schneebergerlab/syri',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('syri.synsearchFunctions', ['syri/pyxFiles/synsearchFunctions.pyx']),
                             Extension('syri.inversions', ['syri/pyxFiles/inversions.pyx']),
                             Extension('syri.tdfunc', ['syri/pyxFiles/tdfunc.pyx']),
                             Extension('syri.findshv', ['syri/pyxFiles/findshv.pyx']),
                             Extension('syri.findsv', ['syri/pyxFiles/findsv.pyx']),
                             Extension('syri.writeout', ['syri/pyxFiles/writeout.pyx'])]),
      packages=["syri", "syri.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['bin/syri', 'bin/chroder'],
      long_description=open('README.rst').read(),
)