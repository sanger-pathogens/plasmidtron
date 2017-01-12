import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

setup(
    name='plasmidtron',
    version='0.0.1',
    description='plasmidtron: a tool to assemble parts of a genome responsible for a trait',
    packages = find_packages(),
    author='Andrew J. Page, Alexander Wailan',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/plasmidtron',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
