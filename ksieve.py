#!/usr/bin/env python3
import argparse
import sys
import os

# input data
parser = argparse.ArgumentParser(
	description = 'Plasmid kmers',
	usage = 'ksieve [options] output_directory file_of_trait_fastqs file_of_nontrait_fastqs')
parser.add_argument('output_directory', help='Output directory')
parser.add_argument('file_of_trait_fastqs', help='File of filenames of trait FASTQs')
parser.add_argument('file_of_nontrait_fastqs', help='File of filenames of nontrait FASTQs')
parser.add_argument('--verbose',  '-v', action='count', help='Turn on debugging', default = 0)
parser.add_argument('--threads',  '-t', help='Number of threads', type=int,  default = 1)
parser.add_argument('--kmer',     '-k', help='Kmer to use, depends on read length', type=int,  default = 81)
options = parser.parse_args()



    # Run kmc to generate kmers for each set of FASTQs

# using Complex, create a file describing merging all the traits into one set, non traits into another set, then subtract.

# Given kmer database of difference, extract raw reads from each sample, but only for trait set.

# SPAdes assembly of each trait set of reads

# cleanup at the end.