from setuptools import setup

setup(name='syri',
        version='0.1',
        description='Identify synteny and structural rearrangements',
        author='Manish Goel',
        author_email='goel@mpipz.mpg.de',
        packages=['syri',
                  'syri.methods'],
        install_requires=['pandas','numpy','python-igraph','matplotlib','multiprocess','biopython','glob2','scipy','datetime'],
        scripts=['bin/SynSearch',
                'bin/scaffoldOrder',
                'bin/getSVs'],
        zip_safe=False)
        
        

