import os
import shutil
import tempfile
import subprocess
import csv
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
		self.temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(self.output_directory))

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
		# Generate a dump file with all kmers
		dump_file = os.path.join(self.temp_working_dir, 'dump.txt')
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
				
		os.remove(dump_file)
		return kmers

	def cleanup(self):
		pass