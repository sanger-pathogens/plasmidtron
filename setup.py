import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = 'x'
if os.path.exists('VERSION'):
  version = open('VERSION').read().strip()

setup(
    name='plasmidtron',
    version=version,
    description='plasmidtron: a tool to assemble parts of a genome responsible for a trait',
	long_description=read('README.md'),
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
		'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
