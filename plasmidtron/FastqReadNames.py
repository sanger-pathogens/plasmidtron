import os
import logging
import re
from Bio import SeqIO
 
'''Given a FASTQ file, extract all of the read names'''
class FastqReadNames:
	def __init__(self,fastq_file, output_readnames_file, verbose, match_both_pairs):
		self.logger = logging.getLogger(__name__)
		self.fastq_file = fastq_file
		self.output_readnames_file = output_readnames_file
		self.match_both_pairs = match_both_pairs
		
		self.verbose = verbose
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
	
	def extract_readnames_from_fastq(self):
		if not os.path.exists(self.fastq_file):
			self.logger.error('Cannot read the FASTQ file %s', self.fastq_file)
			raise
			
		if self.match_both_pairs:
			self.match_both_pairs_filter()
		else:
			self.match_one_pair_filter()
			
	def match_both_pairs_filter(self):
		self.logger.warning("Extracting read names from FASTQ file where both reads must match")
		regex = re.compile(r'/[12]$')
		base_read_names = {}
		with open(self.fastq_file, "r") as fastq_file_input:
			for record in SeqIO.parse(fastq_file_input, "fastq"):
				base_read_name = regex.sub('', record.id)
				
				if base_read_name in base_read_names:
					base_read_names[base_read_name] = record.id
				else:
					base_read_names[base_read_name] = 1
		
		with open(self.output_readnames_file, "w") as readnames_output:
			for base_name, read_name in base_read_names.items():
				if read_name== 1:
					continue
				readnames_output.write(read_name + '\n')
				
	def match_one_pair_filter(self):
		self.logger.warning("Extracting read names from FASTQ file matching 1 or more read")
		with open(self.fastq_file, "r") as fastq_file_input, open(self.output_readnames_file, "w") as readnames_output:
			for record in SeqIO.parse(fastq_file_input, "fastq"):
				readnames_output.write(record.id + '\n')
