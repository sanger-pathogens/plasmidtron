import os
import shutil
import tempfile
import subprocess
import logging
import csv
import sys
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy
from plasmidtron.KmcFasta import KmcFasta
 
'''Given a list of assembly files in FASTA format, output a kmer presence/absense plot'''
class PlotKmers:
	def __init__(self,assemblies,output_directory,threads,kmer,max_kmers_threshold, verbose):
		self.assemblies = []
		for assembly_file in assemblies:
			self.assemblies.append(os.path.abspath(assembly_file))
		
		self.output_directory = output_directory
		self.threads = threads
		self.kmer = kmer 
		self.max_kmers_threshold = max_kmers_threshold
		self.verbose = verbose
		self.temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(self.output_directory))
		
		self.logger = logging.getLogger(__name__)
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
	
	def generate_plot(self):
		kmers_to_assemblies = self.get_kmers_to_assemblies()
		kmer_matrix = self.create_matrix_for_plot(kmers_to_assemblies)
		self.plot_kmer_matrix(kmer_matrix)

	def get_kmers_to_assemblies(self):
		self.logger.warning('Extract kmers from assemblies')
		kmers_to_assemblies = OrderedDict()
		# get kmers for each assembly
		for assembly in self.assemblies:
			self.logger.warning('Finding kmers for assembly %s', assembly)
			kmc_fasta = KmcFasta(self.output_directory, 
								assembly, 
								self.threads, 
								self.kmer,
								1, 
								self.max_kmers_threshold)			
			kmc_fasta.run()
			kmers_for_assembly = self.get_kmers_from_db(kmc_fasta.output_database_name())
			kmc_fasta.cleanup()
			
			# read in kmers to ordered hash
			self.logger.warning('read in kmers to ordered dict')
			for kmer in kmers_for_assembly:
				if kmer not in kmers_to_assemblies:
					kmers_to_assemblies[kmer] = {assembly: 1}
				else:
					kmers_to_assemblies[kmer][assembly] = 1
		return kmers_to_assemblies
		
	def create_matrix_for_plot(self,kmers_to_assemblies):
		self.logger.warning('Create matrix for plot')
		kmer_matrix = numpy.zeros((len(kmers_to_assemblies),len(self.assemblies)))
		kmer_counter = 0
		for kmer, assemblies in kmers_to_assemblies.items():
			for assembly_counter, assembly_file in enumerate(self.assemblies):
				if assembly_file in assemblies and assemblies[assembly_file] == 1:
					kmer_matrix[kmer_counter][assembly_counter] = 1
				else:
					kmer_matrix[kmer_counter][assembly_counter] = 0
			kmer_counter +=1
		return kmer_matrix
		
	def plot_kmer_matrix(self, kmer_matrix):
		plt.matshow(kmer_matrix, cmap=plt.cm.gray)
		plt.show()

	def get_kmers_from_db(self,database):
		self.logger.warning('Get kmers from database %s', database)
		dump_file = 'dump.txt'
		command_to_run =  ' '.join(['kmc_tools', '-t'+str(self.threads), 'transform', database, 'dump', dump_file])
		subprocess.call(command_to_run, shell=True)
		
		kmers = []
		# read in the kmer dump file - sequence then frequency
		# ACGTAAAAAAAA	45
		# ACGTAAAACCCC	21
		with open(dump_file, "r") as kmer_dump:
			kmer_dump_reader = csv.reader(kmer_dump, delimiter='\t')
			for row in kmer_dump_reader:
				kmers.append(row[0])
				
		####os.remove(dump_file)
		return kmers

	def cleanup(self):
		pass