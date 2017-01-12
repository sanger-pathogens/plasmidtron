import os
import logging
import tempfile
 
'''Assemble a filtered sample with SPAdes'''
class SpadesAssembly:
	def __init__(self, sample, output_directory, threads, kmer, spades_exec):
		self.logger = logging.getLogger(__name__)
		self.sample = sample
		self.threads = threads
		self.kmer = kmer
		self.spades_exec = spades_exec
		self.spades_output_directory = os.path.join(self.output_directory, 'spades_'+sample.basename)

	def spades_command(self):
		return ' '.join([self.spades_exec, '--careful', '--only-assembler','-k', str(self.kmer), '-1', self.sample.filtered_forward_file, '-2', self.sample.filtered_reverse_file, '-o', self.spades_output_directory])

	def run(self):
		self.logger.info("Assembling sample : %s" % (self.sample.basename))
		subprocess.call(self.spades_command, shell=True)
