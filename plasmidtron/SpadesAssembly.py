import os
import logging
import tempfile
import subprocess
 
'''Assemble a filtered sample with SPAdes'''
class SpadesAssembly:
	def __init__(self, sample, output_directory, threads, kmer, spades_exec, minimum_length = 400):
		self.logger = logging.getLogger(__name__)
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.minimum_length = minimum_length
		self.spades_exec = spades_exec
		self.spades_output_directory = os.path.join(self.output_directory, 'spades_'+sample.basename)

	def spades_command(self):
		return ' '.join([self.spades_exec, '--careful', '--only-assembler','-k', str(self.kmer), '-1', self.sample.filtered_forward_file, '-2', self.sample.filtered_reverse_file, '-o', self.spades_output_directory])

	def spades_assembly_file(self):
		return os.path.join(self.spades_output_directory,'contigs.fasta')

	def filtered_spades_assembly_file(self):
		return os.path.join(self.spades_output_directory,'filtered_contigs.fasta')
		
	def remove_small_contigs(self):
		with open(self.spades_assembly_file(), "r") as spades_input_file, open(self.filtered_spades_assembly_file, "w") as spades_output_file:
			long_sequences = []
			for record in SeqIO.parse(spades_input_file, "fasta"):
				if len(record.seq) < self.minimum_length:
					long_sequences.append(record)
			spades_output_file.write(short_sequences)
	
	def run(self):
		self.logger.info("Assembling sample" )
		subprocess.call(self.spades_command, shell=True)
