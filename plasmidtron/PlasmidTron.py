import sys
import os
import logging
import time
from plasmidtron.SampleData import SampleData
from plasmidtron.SpreadsheetParser import SpreadsheetParser
from plasmidtron.Kmc import Kmc
from plasmidtron.KmcComplex import KmcComplex
from plasmidtron.FastqReadNames import FastqReadNames
from plasmidtron.KmcFilter import KmcFilter
from plasmidtron.SpadesAssembly import SpadesAssembly
from plasmidtron.Methods import Methods

class PlasmidTron:
	def __init__(self,options):
		self.start_time = int(time.time())
		self.logger = logging.getLogger(__name__)
		self.output_directory        = options.output_directory 
		self.file_of_trait_fastqs    = options.file_of_trait_fastqs
		self.file_of_nontrait_fastqs = options.file_of_nontrait_fastqs
		self.verbose                 = options.verbose
		self.threads                 = options.threads
		self.kmer                    = options.kmer
		self.min_kmers_threshold     = options.min_kmers_threshold
		self.max_kmers_threshold     = options.max_kmers_threshold
		self.spades_exec             = options.spades_exec
		self.min_contig_len          = options.min_contig_len
		self.action                  = options.action

	def run(self):
		os.makedirs(self.output_directory)
		trait_samples = SpreadsheetParser(self.file_of_trait_fastqs).extract_samples()
		nontrait_samples = SpreadsheetParser(self.file_of_nontrait_fastqs).extract_samples()
		
		self.logger.info("Generating a kmer database for each sample")
		kmc_samples =[]
		for set_of_samples in [trait_samples, nontrait_samples]:
			for sample in set_of_samples:
				kmc_sample = Kmc(self.output_directory, sample, self.threads, self.kmer, self.min_kmers_threshold, self.max_kmers_threshold)
				kmc_sample.run()
				kmc_samples.append(kmc_sample)
		
		self.logger.info("Generating a database of kmers which are in the traits but not in the nontraits set")
		kmc_complex = KmcComplex(self.output_directory, self.threads, self.min_kmers_threshold, trait_samples, nontrait_samples, self.action)
		kmc_complex.run()

		kmc_filters = []
		for sample in trait_samples:
			kmc_filter = KmcFilter(sample, self.output_directory, self.threads,kmc_complex.result_database())
			kmc_filter.filter_fastq_file_against_kmers()
			kmc_filters.append(kmc_filter)
		
		spades_assemblies = []
		for sample in trait_samples:
			spades_assembly = SpadesAssembly( sample, self.output_directory, self.threads, self.kmer, self.spades_exec, self.min_contig_len)
			spades_assembly.run()
			spades_assemblies.append(spades_assembly)
			print(spades_assembly.filtered_spades_assembly_file() + '\n')
			sample.cleanup()
			
		method_file = Methods(os.path.join(self.output_directory, 'methods_summary.txt'), trait_samples, nontrait_samples, self.min_kmers_threshold, self.min_contig_len, self.start_time, self.spades_exec)
		method_file.create_file()
		
		if not self.verbose:
			# Delete all sample temp directories
			self.logger.info("Deleting intermediate files, use --verbose if you wish to keep them")
			for kmc_sample in kmc_samples:
				kmc_sample.cleanup()
				
			kmc_complex.cleanup()
			
			for kmc_filter in kmc_filters:
				kmc_filter.cleanup()
			
			for spades_assembly in spades_assemblies:
				spades_assembly.cleanup()
		