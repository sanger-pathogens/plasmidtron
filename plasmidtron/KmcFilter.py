import os
import logging
import tempfile
import subprocess
import shutil
from plasmidtron.FastqReadNames import FastqReadNames
 
'''Given a kmer database extract filter a FASTQ file for a sample'''
class KmcFilter:
	def __init__(self,sample, output_directory, threads,result_database, verbose):
		self.logger = logging.getLogger(__name__)
		self.verbose = verbose
		self.sample = sample
		self.threads = threads
		self.result_database = result_database
		self.temp_working_dir = tempfile.mkdtemp(dir=output_directory)
		self.store_intermediate_files()
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
	
	def store_intermediate_files(self):
		self.intermediate_filtered_fastq = os.path.join(self.temp_working_dir, 'intermediate.fastq')
		self.read_names_file = os.path.join(self.temp_working_dir, 'read_names_file')
		self.sample.filtered_forward_file = os.path.join(self.temp_working_dir, 'sample_1.fastq.gz')
		self.sample.filtered_reverse_file = os.path.join(self.temp_working_dir, 'sample_2.fastq.gz')
		
		if self.verbose:
			self.logger.warning('Filtered FASTQ files for sample %s:\t%s\t%s', self.sample.basename, self.sample.filtered_forward_file, self.sample.filtered_reverse_file)
	
	def redirect_output(self):
		redirect_output_str = ''
		if self.verbose:
			redirect_output_str = ''
		else:
			redirect_output_str = '> /dev/null 2>&1'
		return redirect_output_str
	
	def kmc_filter_command(self):
		return ' '.join(['kmc_tools', '-t'+str(self.threads), 'filter', self.result_database, '@'+self.sample.file_of_fastq_files, self.intermediate_filtered_fastq, self.redirect_output()])
		
	def filtered_fastaq_command(self):
		return ' '.join(['fastaq', 'filter', '--ids_file', self.read_names_file, '--mate_in', self.sample.reverse_file, ' --mate_out', self.sample.filtered_reverse_file, self.sample.forward_file, self.sample.filtered_forward_file, self.redirect_output() ])
	
	def filter_fastq_file_against_kmers(self):
		self.logger.warning("Filter reads against kmer database for sample")
		self.logger.warning('%s',self.kmc_filter_command())
		subprocess.call(self.kmc_filter_command(), shell=True)
	
		# The FASTQ file that comes out of kmc doesnt output all pairs, so we have to refilter it to get all mates.
		fastq_read_names = FastqReadNames(self.intermediate_filtered_fastq, self.read_names_file, self.verbose)
		fastq_read_names.extract_readnames_from_fastq()
	
		# Given a file of read names, pull out the mate paired FASTQ files for the sample
		self.logger.warning('%s',self.filtered_fastaq_command())
		subprocess.call(self.filtered_fastaq_command(), shell=True)
		
	def cleanup(self):
		os.remove(self.intermediate_filtered_fastq)
		os.remove(self.read_names_file)
		shutil.rmtree(self.temp_working_dir)
