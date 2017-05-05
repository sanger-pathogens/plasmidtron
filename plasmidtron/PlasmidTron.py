import sys
import os
import logging
import subprocess
import time
import shutil
import tempfile
from plasmidtron.Kmc import Kmc
from plasmidtron.KmcComplex import KmcComplex
from plasmidtron.KmcFasta import KmcFasta
from plasmidtron.KmcFilter import KmcFilter
from plasmidtron.Methods import Methods
from plasmidtron.SampleData import SampleData
from plasmidtron.SpadesAssembly import SpadesAssembly
from plasmidtron.SpreadsheetParser import SpreadsheetParser
from plasmidtron.PlotKmers import PlotKmers
from plasmidtron.KmcVersionDetect import KmcVersionDetect
from plasmidtron.CommandRunner import CommandRunner

class PlasmidTron:
	def __init__(self,options):
		self.start_time = int(time.time())
		self.logger = logging.getLogger(__name__)
		self.output_directory           = options.output_directory 
		self.file_of_traits       = options.file_of_traits
		self.file_of_nontraits    = options.file_of_nontraits
		self.verbose                    = options.verbose
		self.threads                    = options.threads
		self.kmer                       = options.kmer
		self.min_kmers_threshold        = options.min_kmers_threshold
		self.max_kmers_threshold        = options.max_kmers_threshold
		self.spades_exec                = options.spades_exec
		self.min_contig_len             = options.min_contig_len
		self.action                     = options.action
		self.min_spades_contig_coverage = options.min_spades_contig_coverage
		self.keep_files                 = options.keep_files
		self.plot_filename              = options.plot_filename
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
		self.kmc_major_version = KmcVersionDetect(self.verbose).major_version()
		self.command_runner = CommandRunner( self.output_directory, self.logger, self.threads )

	def generate_kmer_databases(self, trait_samples, nontrait_samples):
		kmc_samples =[]
		kmc_commands_to_run = []
		for set_of_samples in [trait_samples, nontrait_samples]:
			for sample in set_of_samples:
				self.logger.warning('Generating a kmer database for sample %s', sample.basename)
				kmc_sample = Kmc(self.output_directory, sample, 1, self.kmer, self.min_kmers_threshold, self.max_kmers_threshold, self.verbose)
				kmc_sample.create_file_of_file_names(kmc_sample.sample.file_of_fastq_files)
				kmc_commands_to_run.append(kmc_sample.construct_kmc_command())
				kmc_samples.append(kmc_sample)
		
		self.command_runner.run_list_of_commands( kmc_commands_to_run)	
		return kmc_samples
		
	def filter_data_against_kmers(self,trait_samples, result_database):
		kmc_filters = []
		for sample in trait_samples:
			if sample.is_a_fasta():
				self.logger.warning('Not filtering FASTA file for trait kmers %s', sample.basename)
				continue
				
			self.logger.warning('Filtering reads which contain trait kmers %s', sample.basename)
			kmc_filter = KmcFilter(sample, self.output_directory, 1, result_database, self.verbose)
			kmc_filters.append(kmc_filter)
			
		kmc_filter_commands = [ k.kmc_filter_command() for k in kmc_filters ]
		
		self.command_runner.run_list_of_commands( kmc_filter_commands)
		
		# Convert to parallel
		for k in kmc_filters:
			k.extract_read_names_from_fastq()
		
		kmc_fastaq_commands = [ k.filtered_fastaq_command() for k in kmc_filters ]
		self.command_runner.run_list_of_commands( kmc_fastaq_commands)
		
		return kmc_filters

	def assemble_samples(self,trait_samples,keep_files):
		spades_first_assemblies = []
		spades_commands = []

		for sample in trait_samples:
			if sample.is_a_fasta():
				self.logger.warning('Not assembling sample %s', sample.basename)
				continue
			
			self.logger.warning('First assembly with reads only matching kmers %s', sample.basename)
			spades_assembly = SpadesAssembly(	sample, 
												self.output_directory, 
												1, 
												self.kmer, 
												self.spades_exec, 
												self.min_contig_len,
												True,
												self.min_spades_contig_coverage,
												False,
												self.verbose)
			spades_commands.append(spades_assembly.spades_command())
			spades_first_assemblies.append(spades_assembly)
			
		self.command_runner.run_list_of_commands(spades_commands)
		
		kmc_fastas = []
		kmc_fasta_commands = []
		self.logger.warning('Extract the kmers from the first assembly')
		for spades_assembly in spades_first_assemblies:
			sample = spades_assembly.sample
			
			if not os.path.exists(spades_assembly.spades_assembly_file()):
				continue
			
			spades_assembly.remove_small_contigs(spades_assembly.spades_assembly_file(), spades_assembly.filtered_spades_assembly_file())
		
			if os.path.getsize(spades_assembly.filtered_spades_assembly_file()) <= self.min_contig_len:
				self.logger.warning('Not enough data in the 1st assembly after filtering, skipping the rest of the steps %s', sample.basename)
				continue
			
			self.logger.warning('Rescaffold 1st assembly with all reads %s', sample.basename)
			# Next we want to scaffold by using all of the original reads to join up the small contigs.
			# Extract all of the kmers found in the filtered assembly
			self.logger.warning('Extract kmers from assembly %s', sample.basename)
			kmc_fasta = KmcFasta(	self.output_directory, 
									spades_assembly.filtered_spades_assembly_file(), 
									1, 
									self.kmer,
									1, 
									self.max_kmers_threshold,
									self.verbose)
			kmc_fasta.sample = sample			
			kmc_fastas.append(kmc_fasta)
			kmc_fasta_commands.append(kmc_fasta.kmc_command())
			
		self.command_runner.run_list_of_commands(kmc_fasta_commands)	
			
		for spades_assembly in spades_first_assemblies:	
			if not self.keep_files:
				spades_assembly.cleanup()
		
		kmc_filters = []
		kmc_filter_commands = []
		filtered_fastaq_commands = []
		self.logger.warning('Filter the reads against the kmers from the 1st assembly')
		for kmc_fasta in kmc_fastas:
			sample = kmc_fasta.sample
			
			# Pull out any reads matching the kmers found in the assembly
			self.logger.warning('Pull out reads from original FASTQ files matching assembly kmers %s', sample.basename)
			kmc_filter = KmcFilter(	sample, 
									self.output_directory, 
									1, 
									kmc_fasta.output_database_name(),
									self.verbose)
			#kmc_filter.filter_fastq_file_against_kmers()
			kmc_filters.append(kmc_filter)
			kmc_filter_commands.append(kmc_filter.kmc_filter_command())
			filtered_fastaq_commands.append(kmc_filter.filtered_fastaq_command())
		
		self.command_runner.run_list_of_commands(kmc_filter_commands)	
		for k in kmc_filters:
			k.extract_read_names_from_fastq()
		self.command_runner.run_list_of_commands(filtered_fastaq_commands)	

		for kmc_fasta in kmc_fastas:
			if not self.keep_files:
				kmc_fasta.cleanup()
		
		spades_assemblies = []	
		final_spades_commands = []
		self.logger.warning('Reassemble all filtered reads with SPAdes')
		for kmc_filter in kmc_filters:
			sample = kmc_filter.sample
			self.logger.warning('Reassemble with SPAdes %s', sample.basename)
			final_spades_assembly = SpadesAssembly(	sample, 
												self.output_directory, 
												1, 
												self.kmer, 
												self.spades_exec, 
												self.min_contig_len,
												False,
												self.min_spades_contig_coverage,
												True,
												self.verbose)
			spades_assemblies.append(final_spades_assembly)
			final_spades_commands.append(final_spades_assembly.spades_command())
			
		self.command_runner.run_list_of_commands(final_spades_commands)
			
		self.logger.warning('Cleanup kmc filter output')
		for kmc_filter in kmc_filters:
			if not self.keep_files:
				kmc_filter.cleanup()
			
		self.logger.warning('Cleanup final spades output')
		for final_spades_assembly in spades_assemblies:
			if not os.path.exists(final_spades_assembly.spades_assembly_file()):
				continue
			
			final_spades_assembly.remove_small_contigs(final_spades_assembly.spades_assembly_file(), final_spades_assembly.filtered_spades_assembly_file())
			
			if os.path.exists(final_spades_assembly.filtered_spades_assembly_file()):
				print(final_spades_assembly.filtered_spades_assembly_file())
			else:
				self.logger.error('Final SPAdes output doesnt exist %s', final_spades_assembly.filtered_spades_assembly_file())
			
			if not self.keep_files:
				final_spades_assembly.cleanup()

		return spades_assemblies
			
	def run(self):
		self.logger.warning('Using KMC syntax version %s', self.kmc_major_version)
		os.makedirs(self.output_directory)
		trait_samples = SpreadsheetParser(self.file_of_traits, self.verbose).extract_samples()
		nontrait_samples = SpreadsheetParser(self.file_of_nontraits, self.verbose).extract_samples()

		self.logger.warning('Generating kmer databases for all samples')
		kmc_samples = self.generate_kmer_databases(trait_samples, nontrait_samples)
		
		self.logger.warning("Generating a database of kmers which are in the traits but not in the nontraits set")
		kmc_complex = KmcComplex(self.output_directory, self.threads, self.min_kmers_threshold, trait_samples, nontrait_samples, self.action, self.verbose)
		kmc_complex.run()

		kmc_filters = self.filter_data_against_kmers(trait_samples,kmc_complex.result_database())
		
		self.logger.warning('Assembling all of the trait samples')
		spades_assemblies = self.assemble_samples(trait_samples, self.keep_files)
	
		spades_assembly_files = [s.filtered_spades_assembly_file() for s in spades_assemblies if os.path.exists(s.filtered_spades_assembly_file())]
		plot_kmers = PlotKmers( spades_assembly_files,
								self.output_directory,
								self.threads,
								self.kmer,
								self.max_kmers_threshold, 
								self.verbose, 
								self.plot_filename)
		plot_kmers.generate_plot()
			
		method_file = Methods(
						os.path.join(self.output_directory, 'methods_summary.txt'), 
						trait_samples, 
						nontrait_samples, 
						self.min_kmers_threshold, 
						self.min_contig_len, 
						self.start_time, 
						self.spades_exec, 
						self.verbose)
		method_file.create_file()
		self.cleanup(kmc_samples, kmc_complex, kmc_filters, plot_kmers)
		
	def cleanup(self,kmc_samples, kmc_complex, kmc_filters, plot_kmers):
		if not self.keep_files:
			# Delete all sample temp directories
			self.logger.warning("Deleting intermediate files, use --verbose if you wish to keep them")
			for kmc_sample in kmc_samples:
				kmc_sample.cleanup()
				
			kmc_complex.cleanup()
			
			for kmc_filter in kmc_filters:
				kmc_filter.cleanup()
			
			#plot_kmers.cleanup()
		