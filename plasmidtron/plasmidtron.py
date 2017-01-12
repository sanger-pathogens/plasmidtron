#!/usr/bin/env python3
import argparse
import sys
import os
import tempfile
import subprocess
from  plasmidtron.SampleData import SampleData
from  plasmidtron.SpreadsheetParser import SpreadsheetParser

class PlasmidTron:
	def __init__(self,options):
		self.output_directory        = options.output_directory 
		self.file_of_trait_fastqs    = options.file_of_trait_fastqs
		self.file_of_nontrait_fastqs = options.file_of_nontrait_fastqs
		self.verbose                 = options.verbose
		self.threads                 = options.threads
		self.kmer                    = options.kmer
		self.min_kmers_threshold     = options.min_kmers_threshold

	def run():
		if not os.path.exists(self.output_directory):
		    os.makedirs(self.output_directory)
		else:
			sys.exit("The output directory already exists")
		
		trait_samples = SpreadsheetParser(self.file_of_trait_fastqs)
		nontrait_samples = SpreadsheetParser(self.file_of_nontrait_fastqs)
		
		# Run kmc to generate kmers for each set of FASTQs
		for set_of_samples in [trait_samples, nontrait_samples]:
			for sample in set_of_samples:
				temp_working_dir = tempfile.mkdtemp(dir=self.output_directory)
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
		
				kmc_command = "kmc -t"+str(self.threads)+" -ci"+str(self.min_kmers_threshold)+" -k"+ str(self.kmer) +" @"+file_of_fastq_files+" "+ temp_working_dir+"/kmc_"+basename +" "+ temp_working_dir
				print('DEBUG: '+ kmc_command)
				subprocess.call(kmc_command,shell=True)
		
		# using Complex, create a file describing merging all the traits into one set, non traits into another set, then subtract.
		
		# create complex input file
		temp_working_dir = tempfile.mkdtemp(dir=self.output_directory)
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
		
			complex_config_file.write('('+  '+'.join(trait_basenames) +')')
			complex_config_file.write('-')
			complex_config_file.write('('+  '+'.join(nontrait_basenames) +')')
			complex_config_file.write("\n")
			# set operation
			# 
		
		# increase the threshold for kmer counts
		kmc_complex_command = "kmc_tools -t"+str(self.threads)+" complex " + complex_config_filename
		print('DEBUG: '+ kmc_complex_command)
		subprocess.call(kmc_complex_command,shell=True)
		
		
		# For traits only
		# Filter each FASTQ file against kmer database of the differences between traits and non-traits
		
		# Get the names of the filtered reads for each sample (forward and reverse). Remove the /1 or /2 from the end to get common read name. unique.
		# Open original FASTQ files and filter against the read names. This ensures that the reads are still paired, even if only 1 of the reads hit the kmer database.
		for sample in trait_samples:
			temp_working_dir_filter = tempfile.mkdtemp(dir=self.output_directory)
			intermediate_filtered_fastq = temp_working_dir_filter+'/intermediate.fastq'
			read_names_file = temp_working_dir_filter+'/read_names_file'
		
			kmc_filter_command = 'kmc_tools -t'+str(self.threads)+' filter result @'+sample.file_of_fastq_files+' '+intermediate_filtered_fastq
			print('DEBUG: '+ kmc_filter_command)
			subprocess.call(kmc_filter_command,shell=True)
		
			read_names_cmd = 'awk \'NR%4==1\' '+intermediate_filtered_fastq+' | sed \'s!@!!\' > ' + read_names_file
			print('DEBUG: '+ read_names_cmd)
			subprocess.call(read_names_cmd, shell=True)
			
			sample.filtered_forward_file = temp_working_dir_filter+'/sample_1.fastq.gz'
			sample.filtered_reverse_file = temp_working_dir_filter+'/sample_2.fastq.gz'
		
			filtered_fastq_command = 'fastaq filter --ids_file '+read_names_file+' --mate_in '+sample.reverse_file+' --mate_out '+sample.filtered_reverse_file+' '+sample.forward_file+ ' '+sample.filtered_forward_file
			print('DEBUG: '+ filtered_fastq_command)
			subprocess.call(filtered_fastq_command, shell=True)
		
			# delete intermediate fastq file
			os.remove(intermediate_filtered_fastq)
			os.remove(read_names_file)
		
		for sample in trait_samples:
			spades_output_directory = self.output_directory+'/spades_'+sample.basename
			spades_command = 'spades-3.9.0.py --careful --only-assembler -k '+ str(self.kmer) +' -1 '+sample.filtered_forward_file+' -2 '+sample.filtered_reverse_file+' -o '+spades_output_directory
			subprocess.call(spades_command, shell=True)
			print('DEBUG: '+ spades_command)
		