'''Run KMC over a single sample (2 FASTQs or 1 FASTA) to get kmers database'''

import os
import tempfile
import subprocess
import logging
import shutil
from plasmidtron.SampleData import SampleData

class Kmc:
	def __init__(self,output_directory, sample, threads, kmer, min_kmers_threshold, max_kmers_threshold, verbose):
		self.logger = logging.getLogger(__name__)
		self.output_directory = output_directory
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.min_kmers_threshold = min_kmers_threshold
		self.max_kmers_threshold = max_kmers_threshold
		self.temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(output_directory))
		self.populate_database_name()
		self.populate_fofn_name()
		
		self.verbose = verbose
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
		
	def run(self):	
		self.create_file_of_file_names(self.sample.file_of_fastq_files)
		self.logger.warning("Running KMC command" )
		
		subprocess.check_call(self.construct_kmc_command(),shell=True)
	
	def create_file_of_file_names(self, filename):
		with open(filename, 'w') as file_of_sample_fastqs:
			file_of_sample_fastqs.write(self.sample.forward_file + "\n")
			if not self.sample.is_a_fasta():
				file_of_sample_fastqs.write(self.sample.reverse_file + "\n")
			
	def populate_database_name(self):
		self.sample.database_name =  os.path.join(self.temp_working_dir,'kmc_' + self.sample.basename)
		
	def populate_fofn_name(self):
		self.sample.file_of_fastq_files = os.path.join(self.temp_working_dir, 'fofn')
	
	def construct_kmc_command(self):
		redirect_output = ''
		if self.verbose:
			redirect_output = ''
		else:
			redirect_output = '> /dev/null 2>&1'
		
		command_to_run =  " ".join(['kmc', 
			'-t' +  str(self.threads), 
			'-ci' + str(self.min_kmers_value()),
			'-cx' + str(self.max_kmers_threshold),
			'-k' + str(self.kmer),
			self.file_type_option(),
			'@' + self.sample.file_of_fastq_files,
			self.sample.database_name,
			self.temp_working_dir,
			redirect_output
		])
		
		self.logger.warning("Running: "+command_to_run )
		return command_to_run
		
	def min_kmers_value(self):
		if self.sample.is_a_fasta():
			return 1
		else:
			return self.min_kmers_threshold
	
	'''A FASTA file needs a different option for KMC'''
	def file_type_option(self):
		if self.sample.is_a_fasta():
			return '-fm'
		else:
			return '-fq'
				
	def cleanup(self):
		shutil.rmtree(self.temp_working_dir)