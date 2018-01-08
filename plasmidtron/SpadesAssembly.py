import os
import logging
import tempfile
import subprocess
import shutil
import re
from Bio import SeqIO
 
'''Assemble a filtered sample with SPAdes'''
class SpadesAssembly:
	def __init__(self, sample, output_directory, threads, kmer, spades_exec, minimum_length, use_temp_directory,min_spades_contig_coverage, assemble_with_careful, verbose, max_spades_contig_coverage):
		self.logger = logging.getLogger(__name__)
		self.output_directory = output_directory
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.minimum_length = minimum_length
		self.spades_exec = spades_exec
		self.use_temp_directory = use_temp_directory
		self.min_spades_contig_coverage = min_spades_contig_coverage
		self.max_spades_contig_coverage = max_spades_contig_coverage
		self.assemble_with_careful = assemble_with_careful
		self.verbose = verbose
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)

		if not os.path.exists(self.output_directory):
			os.makedirs(self.output_directory)
		self.temp_working_dir = tempfile.mkdtemp(dir=output_directory)
	
	def spades_output_directory(self):
		if self.use_temp_directory:
			return str(self.temp_working_dir)
		else:
			return os.path.join(self.output_directory, 'spades_'+self.sample.basename)

	def spades_command(self):
		redirect_output = ''
		if self.verbose:
			redirect_output = ''
		else:
			redirect_output = '> /dev/null 2>&1'
		
		careful_flag = ''
		if self.assemble_with_careful :
			careful_flag = '--careful'
			
		return ' '.join([self.spades_exec, careful_flag, '--only-assembler','-t', str(self.threads), '-k', str(self.kmer), '-1', self.sample.filtered_forward_file, '-2', self.sample.filtered_reverse_file, '-o', self.spades_output_directory(), redirect_output ])

	def spades_assembly_file(self):
		return os.path.join(self.spades_output_directory(), 'scaffolds.fasta')

	def filtered_spades_assembly_file(self):
		return os.path.join(self.spades_output_directory(), 'filtered_scaffolds.fasta')
	
	def sequence_coverage(self, contig_name):
		# SPAdes add the coverage to the contig name.
		#NODE_447_length_1644_cov_25.5669
		m = re.search('cov_([\d]+)', contig_name)
		if m and m.group(0):
			return int(m.group(1))
		else:
			# return a number big enough that it will always keep the contig
			return self.min_spades_contig_coverage + 1
		
	def remove_small_large_contigs(self,input_file, output_file):
		with open(input_file, "r") as spades_input_file, open(output_file, "w") as spades_output_file:
			sequences = []
			for record in SeqIO.parse(spades_input_file, "fasta"):
				if self.min_spades_contig_coverage != 0 and self.sequence_coverage(record.id) < self.min_spades_contig_coverage:
					self.logger.warning("Excluding contig with low coverage of "+ str(self.sequence_coverage(record.id))+ " from "+ self.spades_assembly_file())
					continue
					
				if self.max_spades_contig_coverage != 0 and self.sequence_coverage(record.id) > self.max_spades_contig_coverage:
					self.logger.warning("Excluding contig with high coverage of "+ str(self.sequence_coverage(record.id))+ " from "+ self.spades_assembly_file())
					continue
				
				if len(record.seq) > self.minimum_length:
					sequences.append(record)
				else:
					self.logger.warning("Excluding contig of length "+ str(len(record.seq))+ " from "+ self.spades_assembly_file() )
			
			SeqIO.write(sequences, spades_output_file, "fasta")
	
	def run(self):
		self.logger.warning("Assembling sample" )
		subprocess.check_call(self.spades_command(), shell=True)
		self.remove_small_large_contigs(self.spades_assembly_file(), self.filtered_spades_assembly_file())

	def cleanup(self):
		self.__remove_spades_subdirectory_if_exists('tmp')
		self.__remove_spades_subdirectory_if_exists('mismatch_corrector' )
		self.__remove_spades_subdirectory_if_exists('misc' )
		self.__remove_spades_subdirectory_if_exists('K'+str(self.kmer) )
		if self.use_temp_directory:
			shutil.rmtree(self.temp_working_dir)
			
	def __remove_spades_subdirectory_if_exists(self,subdirectory):
		if os.path.exists(os.path.join(self.spades_output_directory(), subdirectory )):
			shutil.rmtree(os.path.join(self.spades_output_directory(), subdirectory )) 
			
	
		