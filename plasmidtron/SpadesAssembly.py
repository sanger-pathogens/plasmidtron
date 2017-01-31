import os
import logging
import tempfile
import subprocess
import shutil
from Bio import SeqIO
 
'''Assemble a filtered sample with SPAdes'''
class SpadesAssembly:
	def __init__(self, sample, output_directory, threads, kmer, spades_exec, minimum_length, use_temp_directory):
		self.logger = logging.getLogger(__name__)
		self.output_directory = output_directory
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.minimum_length = minimum_length
		self.spades_exec = spades_exec
		self.use_temp_directory = use_temp_directory
		self.temp_working_dir = tempfile.mkdtemp(dir=output_directory)
	
	def spades_output_directory(self):
		if self.use_temp_directory:
			return str(self.temp_working_dir)
		else:
			return os.path.join(self.output_directory, 'spades_'+sample.basename)

	def spades_command(self):
		return ' '.join([self.spades_exec, '--careful', '--only-assembler','-k', str(self.kmer), '-1', self.sample.filtered_forward_file, '-2', self.sample.filtered_reverse_file, '-o', self.spades_output_directory() ])

	def spades_assembly_file(self):
		return os.path.join(self.spades_output_directory(), 'scaffolds.fasta')

	def filtered_spades_assembly_file(self):
		return os.path.join(self.spades_output_directory(), 'filtered_scaffolds.fasta')
		
	def remove_small_contigs(self,input_file, output_file):
		with open(input_file, "r") as spades_input_file, open(output_file, "w") as spades_output_file:
			sequences = []
			for record in SeqIO.parse(spades_input_file, "fasta"):
				if len(record.seq) > self.minimum_length:
					sequences.append(record)
			
			SeqIO.write(sequences, spades_output_file, "fasta")
	
	def run(self):
		self.logger.info("Assembling sample" )
		subprocess.call(self.spades_command(), shell=True)
		self.remove_small_contigs(self.spades_assembly_file(), self.filtered_spades_assembly_file())

	def cleanup(self):
		shutil.rmtree(os.path.join(self.spades_output_directory(), 'tmp' ))
		shutil.rmtree(os.path.join(self.spades_output_directory(), 'mismatch_corrector' ))
		shutil.rmtree(os.path.join(self.spades_output_directory(), 'misc' ))
		shutil.rmtree(os.path.join(self.spades_output_directory(), 'K'+self.kmer ))
		if self.use_temp_directory:
			shutil.rmtree(self.temp_working_dir)
		