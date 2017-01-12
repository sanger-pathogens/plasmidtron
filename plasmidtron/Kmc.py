'''Run KMC over a single sample (2 FASTQs) to get kmers database'''

import os
import tempfile
import subprocess
import logging
import shutil
from plasmidtron.SampleData import SampleData

class Kmc:
	def __init__(self,output_directory, sample, threads, kmer, min_kmers_threshold):
		self.logger = logging.getLogger(__name__)
		self.output_directory = output_directory
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.min_kmers_threshold = min_kmers_threshold
		self.temp_working_dir = tempfile.mkdtemp(dir=output_directory)
		self.populate_database_name()
		self.populate_fofn_name()
		
	def run(self):	
		self.create_file_of_file_names(self.sample.file_of_fastq_files)
		self.logger.info("KMC command: %s" % self.construct_kmc_command())
		
		subprocess.call(self.construct_kmc_command(),shell=True)
	
	def create_file_of_file_names(self, filename):
		with open(filename, 'w') as file_of_sample_fastqs:
			file_of_sample_fastqs.write(self.sample.forward_file + "\n")
			file_of_sample_fastqs.write(self.sample.reverse_file + "\n")
			
	def populate_database_name(self):
		self.sample.database_name =  os.path.join(self.temp_working_dir,'kmc_' + self.sample.basename)
		
	def populate_fofn_name(self):
		self.sample.file_of_fastq_files = os.path.join(self.temp_working_dir, 'fofn')
	
	def construct_kmc_command(self):
		return " ".join['kmc', 
			'-t' +  str(self.threads), 
			'-ci' + str(self.min_kmers_threshold),
			'-k' + str(self.kmer),
			'@' + self.sample.file_of_fastq_files,
			self.sample.database_name,
			self.temp_working_dir
		]
	def cleanup(self):
		shutil.rmtree(self.temp_working_dir)