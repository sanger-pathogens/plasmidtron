import os
import logging
from Bio import SeqIO
 
'''Given a FASTQ file, extract all of the read names'''
class FastqReadNames:
	def __init__(self,fastq_file, output_readnames_file):
		self.logger = logging.getLogger(__name__)
		self.fastq_file = fastq_file
		self.output_readnames_file = output_readnames_file
	
	def extract_readnames_from_fastq(self):
		self.logger.info("Extracting read names from '%s' FASTQ file and saving them to '%s'" % (self.fastq_file, self.output_readnames_file))
		with open(self.fastq_file, "r") as fastq_file_input, open(self.output_readnames_file, "w") as readnames_output:
			for record in SeqIO.parse(fastq_file_input, "fastq"):
				readnames_output.write(record.id + '\n')

