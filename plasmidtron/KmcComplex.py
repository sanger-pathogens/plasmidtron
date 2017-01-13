'''Take in trait and non trait databases and do a complex filter on the kmers'''

import os
import tempfile
import subprocess
import logging
import shutil
from plasmidtron.SampleData import SampleData

class KmcComplex:
	def __init__(self,output_directory, threads, min_kmers_threshold, trait_samples, nontrait_samples):
		self.logger = logging.getLogger(__name__)
		self.output_directory = output_directory
		self.threads = threads
		self.min_kmers_threshold = min_kmers_threshold
		self.trait_samples = trait_samples
		self.nontrait_samples = nontrait_samples
		
		self.temp_working_dir = tempfile.mkdtemp(dir=output_directory)
		self.complex_config_filename = self.generate_complex_config_filename()

	def generate_complex_config_filename(self):
		return os.path.join(self.temp_working_dir, 'complex_config_file')

	def sample_definition_line(self, sample):
		return ' '.join([sample.basename.replace('#','_'), '=', sample.database_name])
		
	def create_config_file(self):
		self.logger.info("Creating config file for 'complex' task")
		with open(self.complex_config_filename, 'w') as complex_config_file:
			complex_config_file.write('INPUT:\n')
			complex_config_file.write(self.sample_definitions_str())
			complex_config_file.write('OUTPUT:\n')
			complex_config_file.write( self.samples_to_set_operation_str() + '\n')
			complex_config_file.write('OUTPUT_PARAMS:\n')
			complex_config_file.write( self.output_parameters_str() + '\n')

	def sample_definitions_str(self):
		sample_definition_lines = ''
		
		for set_of_samples in [self.trait_samples, self.nontrait_samples]:
			for sample in set_of_samples:
				sample_definition_lines += self.sample_definition_line(sample) + "\n"
		return sample_definition_lines
		
	def result_database(self):
		return os.path.join(self.temp_working_dir, 'result')

	def output_parameters_str(self):
		return ' '.join(['-ci'+str(self.min_kmers_threshold)])

	def samples_to_set_operation_str(self):
		trait_basenames = []
		for sample in self.trait_samples:
			trait_basenames.append(sample.basename.replace('#','_'))
			
		nontrait_basenames = []
		for sample in self.nontrait_samples:
			nontrait_basenames.append(sample.basename.replace('#','_'))
	
		set_operation_str  = 'result = '
		set_operation_str += '('+  '+'.join(trait_basenames) +')'
		set_operation_str += '-'
		set_operation_str += '('+  '+'.join(nontrait_basenames) +')'
		return set_operation_str
	
	def kmc_complex_command(self):
		return " ".join(['kmc_tools', 
			'-t' +  str(self.threads),
			'complex',
			self.complex_config_filename ])
	
	def run(self):
		self.create_config_file()
		self.logger.info("Running KMC complex command")
		# kmc_tools doesnt allow for paths in the output database name (even though they say they do) so change working directory
		# to prevent temp files polluting CWD
		original_cwd = os.getcwd()
		os.chdir(self.temp_working_dir)
		subprocess.call(self.kmc_complex_command(), shell=True)
		os.chdir(original_cwd)
		
	def cleanup(self):
		shutil.rmtree(self.temp_working_dir)
		