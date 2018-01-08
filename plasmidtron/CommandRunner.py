import os
import logging
import subprocess
import shutil
import math
import tempfile
 
'''Given a list of commands run single threaded or in parallel'''
class CommandRunner:
	def __init__(self, output_directory, logger, threads):
		self.logger = logger
		self.threads = threads
		self.output_directory = output_directory
	
	def run_list_of_kmc_commands(self, commands_to_run):
		if self.threads > 1:
			self.run_with_parallel(commands_to_run, self.kmc_processes())
		else:
			self.run_sequentially(commands_to_run)
	
	def run_list_of_commands(self, commands_to_run):
		if self.threads > 1:
			self.run_with_parallel(commands_to_run, self.threads)
		else:
			self.run_sequentially(commands_to_run)
			
	def run_sequentially(self, commands_to_run):
		for c in commands_to_run:
			self.logger.warning('GNU parallel command to run %s', c)
			subprocess.check_call(c, shell=True)
	
	'''KMC handles multithreading badly. So give each process 2 threads and limit the overall total independant processes'''
	def kmc_threads(self):
		if self.threads >= 2:
			return 2
		else:
			return 1
		
	def kmc_processes(self):
		if self.threads >= 2:
			return int(math.floor(self.threads/2))
		else:
			return 1
	
	'''Use GNU parallel to manage parallel processing'''	
	def run_with_parallel(self, commands_to_run, processes_in_parallel):
		temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(self.output_directory))
		file_of_commands = os.path.join(temp_working_dir,'commands_to_run')
		with open(file_of_commands, 'w') as commands_file:
			for c in commands_to_run:
				self.logger.warning('Command to run %s', c)
				commands_file.write(c + "\n")
			
		gnu_parallel_command = ' '.join(['parallel', '--gnu', '-j '+ str(processes_in_parallel), '<',file_of_commands])
		self.logger.warning('GNU parallel command to run %s', gnu_parallel_command)
		subprocess.check_call(gnu_parallel_command, shell=True)
		shutil.rmtree(temp_working_dir)
