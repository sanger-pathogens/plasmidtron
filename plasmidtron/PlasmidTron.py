import argparse
import sys
import os
import tempfile
import subprocess
import logging
from plasmidtron.SampleData import SampleData
from plasmidtron.SpreadsheetParser import SpreadsheetParser
from plasmidtron.Kmc import Kmc
from plasmidtron.KmcComplex import KmcComplex
from plasmidtron.FastqReadNames import FastqReadNames
from plasmidtron.KmcFilter import KmcFilter
from plasmidtron.SpadesAssembly import SpadesAssembly

class PlasmidTron:
	def __init__(self,options):
		self.logger = logging.getLogger(__name__)
		self.output_directory        = options.output_directory 
		self.file_of_trait_fastqs    = options.file_of_trait_fastqs
		self.file_of_nontrait_fastqs = options.file_of_nontrait_fastqs
		self.verbose                 = options.verbose
		self.threads                 = options.threads
		self.kmer                    = options.kmer
		self.min_kmers_threshold     = options.min_kmers_threshold
		self.spades_exec             = options.spades_exec
		self.min_contig_len          = options.min_contig_len

	def run(self):
		if not os.path.exists(self.output_directory):
		    os.makedirs(self.output_directory)
		else:
			sys.exit("The output directory already exists")
		
		trait_samples = SpreadsheetParser(self.file_of_trait_fastqs).extract_samples()
		nontrait_samples = SpreadsheetParser(self.file_of_nontrait_fastqs).extract_samples()
		
		self.logger.info("Generating a kmer database for each sample")
		kmc_samples =[]
		for set_of_samples in [trait_samples, nontrait_samples]:
			for sample in set_of_samples:
				kmc_sample = Kmc(self.output_directory, sample, self.threads, self.kmer, self.min_kmers_threshold)
				kmc_sample.run()
				kmc_samples.append(kmc_sample)
		
		self.logger.info("Generating a database of kmers which are in the traits but not in the nontraits set")
		kmc_complex = KmcComplex(self.output_directory, self.threads, self.min_kmers_threshold, trait_samples, nontrait_samples)
		kmc_complex.run()
		
		# Delete all sample temp directories
		self.logger.info("Deleting individual kmer databases for samples")
		for kmc_sample in kmc_samples:
			kmc_sample.cleanup()
		
		for sample in trait_samples:
			kmc_filter = KmcFilter(sample, self.output_directory, self.threads)
			kmc_filter.filter_fastq_file_against_kmers()
		
		for sample in trait_samples:
			spades_assembly = SpadesAssembly( sample, self.output_directory, self.threads, self.kmer, self.spades_exec, self.min_contig_len)
			spades_assembly.run()
			spades_assembly.remove_small_contigs()
			print(spades_assembly.filtered_spades_assembly_file() + '\n')
			sample.cleanup()
		