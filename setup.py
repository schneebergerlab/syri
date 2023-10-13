#!/usr/bin/env python3

import glob
from setuptools import setup, Extension
import numpy
from Cython.Build import cythonize

setup(name="syri",
      description='Synteny and rearrangement identifier between whole-genome assemblies',
      author='Manish Goel',
      author_email='goel@mpipz.mpg.de',
      url='https://github.com/schneebergerlab/syri',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([
          Extension(f"syri.{name.split('/')[-1].split('.')[0]}", [name])
          for name in glob.iglob('syri/pyxFiles/*.pyx')
          ]),
      packages=["syri", "syri.scripts"],
      include_dirs=[numpy.get_include()],
      entry_points={"console_scripts": ["syri=syri.scripts.syri:main", "chroder=syri.scripts.chroder:main"]},
      long_description=open('README.rst').read(),
)
