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
		self.filtered_forward_file=''
		self.filtered_reverse_file=''

# input data
parser = argparse.ArgumentParser(
	description = 'Plasmid kmers',
	usage = 'ksieve [options] output_directory file_of_trait_fastqs file_of_nontrait_fastqs')
parser.add_argument('output_directory', help='Output directory')
parser.add_argument('file_of_trait_fastqs', help='File of filenames of trait FASTQs')
parser.add_argument('file_of_nontrait_fastqs', help='File of filenames of nontrait FASTQs')
parser.add_argument('--verbose',  '-v', action='count', help='Turn on debugging', default = 0)
parser.add_argument('--threads',  '-t', help='Number of threads', type=int,  default = 1)
parser.add_argument('--kmer',	 '-k', help='Kmer to use, depends on read length', type=int,  default = 81)
parser.add_argument('--min_kmers_threshold',	 '-m', help='Exclude k-mers occurring less than this', type=int,  default = 20)
options = parser.parse_args()

if not os.path.exists(options.output_directory):
    os.makedirs(options.output_directory)
else:
	sys.exit("The output directory already exists")

trait_samples = []
nontrait_samples = []
# read in the metadata about the fastq files
with open(options.file_of_trait_fastqs) as csvfile:
	spreadsheetreader = csv.reader(csvfile, delimiter = ',')
	for row in spreadsheetreader:
		forward_file = row[0]
		reverse_file = row[1]

		if not os.path.exists(forward_file) or not os.path.exists(reverse_file):
			 sys.exit( "The input files do not exist: "+forward_file+ " "+reverse_file)
		trait_samples.append( SampleData(forward_file,reverse_file) )

# read in the metadata about the fastq files
with open(options.file_of_nontrait_fastqs) as csvfile:
	spreadsheetreader = csv.reader(csvfile, delimiter = ',')
	for row in spreadsheetreader:
		forward_file = row[0]
		reverse_file = row[1]

		if not os.path.exists(forward_file) or not os.path.exists(reverse_file):
			 sys.exit( "The input files dont exist: ")
		nontrait_samples.append( SampleData(forward_file,reverse_file) )

# Run kmc to generate kmers for each set of FASTQs
for set_of_samples in [trait_samples, nontrait_samples]:
	for sample in set_of_samples:
		temp_working_dir = tempfile.mkdtemp(dir=options.output_directory)
		# TODO make sure to strip off path
		basename = sample.forward_file.replace('_1.fastq.gz','')
		sample.basename = basename

		with open(temp_working_dir+'/fofn', 'w') as file_of_sample_fastqs:
			file_of_sample_fastqs.write(sample.forward_file + "\n")
			file_of_sample_fastqs.write(sample.reverse_file + "\n")

		database_name =  temp_working_dir+"/kmc_"+basename
		file_of_fastq_files = temp_working_dir+"/fofn"

		sample.database_name = database_name
		sample.file_of_fastq_files =file_of_fastq_files

		kmc_command = "kmc -t"+str(options.threads)+" -ci"+str(options.min_kmers_threshold)+" -k"+ str(options.kmer) +" @"+file_of_fastq_files+" "+ temp_working_dir+"/kmc_"+basename +" "+ temp_working_dir
		print('DEBUG: '+ kmc_command)
		subprocess.call(kmc_command,shell=True)

# using Complex, create a file describing merging all the traits into one set, non traits into another set, then subtract.

# create complex input file
temp_working_dir = tempfile.mkdtemp(dir=options.output_directory)
complex_config_filename = temp_working_dir+'/complex_config_file'
with open(complex_config_filename, 'w') as complex_config_file:
	complex_config_file.write("INPUT:\n")

	# write out each database name
	for set_of_samples in [trait_samples, nontrait_samples]:
		for sample in set_of_samples:
			complex_config_file.write(sample.basename.replace('#','_')+' = '+sample.database_name+"\n")

	complex_config_file.write("OUTPUT:\n")
	complex_config_file.write("result =")
	trait_basenames = []
	for sample in trait_samples:
		trait_basenames.append(sample.basename.replace('#','_'))
	nontrait_basenames = []
	for sample in nontrait_samples:
		nontrait_basenames.append(sample.basename.replace('#','_'))

	complex_config_file.write('('+  '*'.join(trait_basenames) +')')
	complex_config_file.write('-')
	complex_config_file.write('('+  '+'.join(nontrait_basenames) +')')
	complex_config_file.write("\n")
	# set operation
	# shoudl there be an intersection between the trait set???
	# 

kmc_complex_command = "kmc_tools -t"+str(options.threads)+" complex " + complex_config_filename
print('DEBUG: '+ kmc_complex_command)
subprocess.call(kmc_complex_command,shell=True)


# For traits only
# Filter each FASTQ file against kmer database of the differences between traits and non-traits

# Get the names of the filtered reads for each sample (forward and reverse). Remove the /1 or /2 from the end to get common read name. unique.
# Open original FASTQ files and filter against the read names. This ensures that the reads are still paired, even if only 1 of the reads hit the kmer database.
for sample in trait_samples:
	temp_working_dir_filter = tempfile.mkdtemp(dir=options.output_directory)
	intermediate_filtered_fastq = temp_working_dir_filter+'/intermediate.fastq'
	read_names_file = temp_working_dir_filter+'/read_names_file'

	kmc_filter_command = 'kmc_tools -t'+str(options.threads)+' filter result @'+sample.file_of_fastq_files+' '+intermediate_filtered_fastq
	print('DEBUG: '+ kmc_filter_command)
	subprocess.call(kmc_filter_command,shell=True)

	read_names_cmd = 'awk \'NR%4==1\' '+intermediate_filtered_fastq+' | sed \'s!@!!\' > ' + read_names_file
	print('DEBUG: '+ read_names_cmd)
	subprocess.call(read_names_cmd, shell=True)
	
	filtered_forward_file = temp_working_dir_filter+'/sample_1.fastq.gz'
	sample.filtered_forward_file = filtered_forward_file

	filtered_reverse_file = temp_working_dir_filter+'/sample_2.fastq.gz'
	sample.filtered_reverse_file = filtered_reverse_file

	filtered_fastq_command = 'fastaq filter --ids_file '+read_names_file+' --mate_in '+sample.reverse_file+' --mate_out '+sample.filtered_reverse_file+' '+sample.forward_file+ ' '+sample.filtered_forward_file
	print('DEBUG: '+ filtered_fastq_command)
	subprocess.call(filtered_fastq_command, shell=True)

	# delete intermediate fastq file
	os.remove(intermediate_filtered_fastq)
	os.remove(read_names_file)

for sample in trait_samples:
	spades_output_directory = options.output_directory+'/spades_'+sample.basename
	spades_command = 'spades-3.9.0.py --careful --only-assembler -k '+ str(options.kmer) +' -1 '+sample.filtered_forward_file+' -2 '+sample.filtered_reverse_file+' -o '+spades_output_directory
	subprocess.call(spades_command, shell=True)
	print('DEBUG: '+ spades_command)
