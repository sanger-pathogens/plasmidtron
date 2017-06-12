import os
import shutil
import tempfile
import subprocess
import logging
import csv
import math
import sys
from collections import OrderedDict

import matplotlib
# Remove the need for a DISPLAY
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
import numpy

from plasmidtron.KmcFasta import KmcFasta
from plasmidtron.KmcVersionDetect import KmcVersionDetect
from plasmidtron.CommandRunner import CommandRunner
 
'''Given a list of assembly files in FASTA format, output a kmer presence/absense plot'''
class PlotKmers:
	def __init__(self,assemblies,output_directory,threads,kmer,max_kmers_threshold, verbose, kmer_plot_filename, max_kmers_to_show = 100000):
		self.kmer_plot_filename = kmer_plot_filename
		self.threads = threads
		self.kmer = kmer 
		self.max_kmers_threshold = max_kmers_threshold
		self.verbose = verbose
		self.max_kmers_to_show = max_kmers_to_show
		self.output_directory = output_directory
		
		if not os.path.exists(self.output_directory):
			os.makedirs(self.output_directory)
		
		self.assemblies = []
		self.assembly_index = {}
		for counter, assembly_file in enumerate(assemblies):
			self.assemblies.append(os.path.abspath(assembly_file))
			self.assembly_index[os.path.abspath(assembly_file)] = counter

		self.temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(self.output_directory))

		self.logger = logging.getLogger(__name__)
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
		self.kmc_major_version = KmcVersionDetect(self.verbose).major_version()
		self.command_runner = CommandRunner( self.output_directory, self.logger, self.threads )
		
			
	def output_filename(self):
		return os.path.join(self.output_directory,self.kmer_plot_filename)
	
	def generate_plot(self):
		self.logger.warning('Using KMC syntax version %s', self.kmc_major_version)
		kmers_to_assemblies = self.get_kmers_to_assemblies()
		kmer_matrix = self.create_matrix_for_plot(kmers_to_assemblies)
		self.plot_kmer_matrix(kmer_matrix)
		self.cleanup()

	def get_kmers_to_assemblies(self):
		self.logger.warning('Extract kmers from assemblies')
		kmers_to_assemblies = OrderedDict()
		# get kmers for each assembly
		kmc_fastas = []
		kmc_fasta_commands = []
		for assembly in self.assemblies:
			self.logger.warning('Finding kmers for assembly %s', assembly)
			kmc_fasta = KmcFasta(self.output_directory, 
								assembly, 
								1, 
								self.kmer,
								1, 
								self.max_kmers_threshold,
								self.verbose)						
			kmc_fastas.append(kmc_fasta)
			kmc_fasta_commands.append(kmc_fasta.kmc_command())
			
		self.command_runner.run_list_of_commands( kmc_fasta_commands )		
			
		for kmc_fasta in kmc_fastas:
			kmers_for_assembly = self.get_kmers_from_db(kmc_fasta.output_database_name())
			kmc_fasta.cleanup()
			
			assembly_idx = self.assembly_index[kmc_fasta.input_filename]
			# read in kmers to ordered hash
			self.logger.warning('read in kmers to ordered dict')
			for kmer in kmers_for_assembly:
				if kmer not in kmers_to_assemblies:
					kmers_to_assemblies[kmer] = {assembly_idx: 1}
				else:
					kmers_to_assemblies[kmer][assembly_idx] = 1
		return kmers_to_assemblies
		
	def create_matrix_for_plot(self,kmers_to_assemblies):
		self.logger.warning('Create matrix for plot')
		kmer_matrix = numpy.zeros((len(self.assemblies),len(kmers_to_assemblies)))
		kmer_counter = 0
		for kmer, assemblies in kmers_to_assemblies.items():
			for assembly_counter, assembly_file in enumerate(self.assemblies):
				if assembly_counter in assemblies and assemblies[assembly_counter] == 1:
					kmer_matrix[assembly_counter][kmer_counter] = 1
				else:
					kmer_matrix[assembly_counter][kmer_counter] = 0
			kmer_counter +=1
		return kmer_matrix
		
	def plot_kmer_matrix(self, kmer_matrix):
		interval_size = 1
		if len(kmer_matrix) > self.max_kmers_to_show:
			interval_size = int(math.ceil(len(kmer_matrix)/self.max_kmers_to_show))
			
		base_sample_names = []
		common_prefix = os.path.commonprefix(self.assemblies)		
		for assembly in self.assemblies:
			stripped_assembly = assembly.replace(common_prefix, '').replace('/filtered_scaffolds.fasta','')
			base_sample_names.append(os.path.basename(stripped_assembly))
		
		fig, ax = plt.subplots()
		ax.matshow(kmer_matrix, cmap=plt.cm.Greys, aspect='auto')
		
		plt.ylabel('Samples')
		plt.xlabel('k-mers')
        
		ax.set_yticks(range(len(base_sample_names)))
		ax.set_yticklabels(base_sample_names)
		make_axes_area_auto_adjustable(ax)
		plt.savefig(self.output_filename())

	def get_kmers_from_db(self,database):
		self.logger.warning('Get kmers from database %s', database)
		dump_file = os.path.join(self.temp_working_dir,'dump.txt')

		command_to_run = ''
		if self.kmc_major_version == 2:
			command_to_run =  ' '.join(['kmc_dump', database, dump_file, self.redirect_output()])
		else:
			command_to_run =  ' '.join(['kmc_tools', '-t'+str(self.threads), 'transform', database, 'dump', dump_file, self.redirect_output()])
		self.logger.warning('KMC dump command: %s', command_to_run)
		subprocess.call(command_to_run, shell=True)
		
		kmers = []
		# read in the kmer dump file - sequence then frequency
		# ACGTAAAAAAAA	45
		# ACGTAAAACCCC	21
		with open(dump_file, "r") as kmer_dump:
			kmer_dump_reader = csv.reader(kmer_dump, delimiter='\t')
			for row in kmer_dump_reader:
				kmers.append(row[0])
				
		os.remove(dump_file)
		return kmers
		
	def redirect_output(self):
		redirect_output_str = ''
		if self.verbose:
			redirect_output_str = ''
		else:
			redirect_output_str = '> /dev/null 2>&1'
		return redirect_output_str

	def cleanup(self):
		shutil.rmtree(self.temp_working_dir)
		