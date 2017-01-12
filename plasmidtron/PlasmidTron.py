import argparse
import sys
import os
import tempfile
import subprocess
import logging
from plasmidtron.SampleData import SampleData
from plasmidtron.SpreadsheetParser import SpreadsheetParser
from plasmidtron.Kmc import Kmc
from plasmidtron.KmcComplex import KmcComplex

class PlasmidTron:
	def __init__(self,options):
		self.logger = logging.getLogger(__name__)
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
		
		self.logger.info("Generating a kmer database for each sample"))
		for set_of_samples in [trait_samples, nontrait_samples]:
			for sample in set_of_samples:
				kmc_sample = Kmc(self.output_directory, sample, self.threads, self.kmer, self.min_kmers_threshold)
				kmc_sample.run()
		
		self.logger.info("Generating a database of kmers which are in the traits but not in the nontraits set"))
		kmc_complex = KmcComplex(self.output_directory, self.threads, self.min_kmers_threshold, trait_samples, nontrait_samples)
		kmc_complex.run()
		
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
		