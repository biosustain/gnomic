# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import codecs
from setuptools import setup, find_packages

setup(
    name='gnomic',
    version='0.6.1',
    packages=find_packages(exclude=['*tests*']),
    license='Apache',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='A grammar for describing microbial genotypes and phenotypes',
    long_description=codecs.open('README.rst', encoding='utf-8').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=[
        'grako>=3.18.1',
        'six>=1.8.0'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Other Environment',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    extras_require={
        'docs': ['Sphinx', 'sphinx-rtd-theme'],
    }
)
