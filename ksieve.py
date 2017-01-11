#!/usr/bin/env python3
import argparse
import sys
import os
import csv
import tempfile
import subprocess

class SampleData:
    def __init__(self,forward_file, reverse_file):
            self.forward_file = forward_file
			self.reverse_file = reverse_file
			self.database_name = ''
			self.file_of_fastq_files = ''
			self.basename=''
	



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

trait_samples = []
nontrait_samples = []

# read in the metadata about the fastq files
with open(options.file_of_trait_fastqs) as csvfile:
    spreadsheetreader = csv.reader(csvfile, delimiter = ',')
    for row in spreadsheetreader:
        forward_fastq_file = row[0]
        reverse_fastq_file = row[1]
		
		if not os.path.exists(forward_fastq_file) or os.path.exists(reverse_fastq_file):
			 sys.exit( "The input files dont exist")
		trait_samples.append( SampleData(forward_fastq_file,reverse_fastq_file) )

# read in the metadata about the fastq files
with open(options.file_of_nontrait_fastqs) as csvfile:
    spreadsheetreader = csv.reader(csvfile, delimiter = ',')
    for row in spreadsheetreader:
        forward_fastq_file = row[0]
        reverse_fastq_file = row[1]
		
		if not os.path.exists(forward_fastq_file) or os.path.exists(reverse_fastq_file):
			 sys.exit( "The input files dont exist")
		nontrait_samples.append( SampleData(forward_fastq_file,reverse_fastq_file) )	

# Run kmc to generate kmers for each set of FASTQs
for set_of_samples in [trait_samples, nontrait_samples]:
	for sample in set_of_samples:
		temp_working_dir = tempfile.mkdtemp(dir=os.getcwd())
		# TODO make sure to strip off path
		basename = sample.forward_fastq_file.replace('_1.fastq.gz','')
		sample.basename(basename)
		
		with open(temp_working_dir+'/fofn', 'w') as file_of_sample_fastqs:
			file_of_sample_fastqs.write(sample.forward_fastq_file + "\n")
			file_of_sample_fastqs.write(sample.reverse_fastq_file + "\n")
		
		database_name =  temp_working_dir+"/kmc_"+basename
		file_of_fastq_files = temp_working_dir+"/fofn"
		
		sample.database_name(database_name)
		sample.file_of_fastq_files(file_of_fastq_files)
		
		kmc_command = "kmc -ci20 -k"+options.kmer+" @"+file_of_fastq_files+" "+ database_name
		subprocess.call(kmc_command,shell=True)
	
# using Complex, create a file describing merging all the traits into one set, non traits into another set, then subtract.
	
# create complex input file
temp_working_dir = tempfile.mkdtemp(dir=os.getcwd())
complex_config_filename = temp_working_dir+'/complex_config_file'
with open(complex_config_filename, 'w') as complex_config_file:
	complex_config_file.write("INPUT:\n")
	
	# write out each database name
	for set_of_samples in [trait_samples, nontrait_samples]:
		for sample in set_of_samples:
			complex_config_file.write(sample.basename+' = '+sample.database_name+"\n")
	
	complex_config_file.write("OUTPUT:\n")
	complex_config_file.write("result =")
	trait_basenames = []
	for sample in trait_samples:
		trait_basenames.append(sample.basename)
	nontrait_basenames = []
	for sample in nontrait_samples:
		nontrait_basenames.append(sample.basename)	
		
	complex_config_file.write('('+  trait_basenames.join('+') +')')
	complex_config_file.write('-')
	complex_config_file.write('('+  nontrait_basenames.join('+') +')')
	complex_config_file.write("\n")
	# set operation

kmc_complex_command = "kmc_tools complex " + complex_config_filename
subprocess.call(kmc_complex_command,shell=True)


# For traits only
# Filter each FASTQ file against kmer database of the differences between traits and non-traits

# Get the names of the filtered reads for each sample (forward and reverse). Remove the /1 or /2 from the end to get common read name. unique. 
# Open original FASTQ files and filter against the read names. This ensures that the reads are still paired, even if only 1 of the reads hit the kmer database.



# Given kmer database of difference, extract raw reads from each sample, but only for trait set.

# SPAdes assembly of each trait set of reads

# cleanup at the end.