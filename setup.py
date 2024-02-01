#!/usr/bin/env python

from distutils.core import setup

setup(name='motifSearch',
      version='1.0',
      description='For counting motifs in genomic DNA',
      author='Callum T-C',
      author_email='email@example.com',
      url='https://github.com/callum-jpg/motifSearch',
      install_requires=[
      'biopython==1.77',
      'numpy>=1.18.1',
      'pandas==1.1.2',
      'scikit-learn==1.0.1',
      'matplotlib==3.1.3'
      ]
     )
