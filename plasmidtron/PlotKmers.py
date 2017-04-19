import os
import shutil
from plasmidtron.KmcFasta import KmcFasta
 
'''Given a list of assembly files in FASTA format, output a kmer presence/absense plot'''
class PlotKmers:
	def __init__(self,assemblies,output_directory,threads,kmer,max_kmers_threshold, verbose):
		self.assemblies = assemblies
		self.output_directory = output_directory
		self.threads = threads
		self.kmer = kmer 
		self.max_kmers_threshold = max_kmers_threshold
		self.verbose = verbose

	def run(self):
		assembly_kmers = []
		# get kmers for each assembly
		for assembly in self.assemblies:
			kmc_fasta = KmcFasta(self.output_directory, 
								assembly, 
								self.threads, 
								self.kmer,
								1, 
								self.max_kmers_threshold)
			assembly_kmers.append(kmc_fasta)
			
		# read in kmers to ordered hash
		
	
	def get_kmers_from_db(self,database):
				command_to_run =  ' '.join(['kmc_tools', '-t'+str(self.threads), 'transform', database, 'dump', 'dump.txt'])

	def cleanup(self):
		pass