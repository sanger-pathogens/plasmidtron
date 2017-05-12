import os
import argparse

class InputTypes:
	
	'''All of the input files listed should exist'''
	def is_input_file_valid(filename):
		if not os.path.exists(filename):
			raise argparse.ArgumentTypeError('Cannot access input file')
		return filename
		
	def is_output_directory_valid(filename):
		if os.path.exists(filename):
			raise argparse.ArgumentTypeError("The output directory already exists")
		return filename
	
	def is_file_of_traits_valid(filename):
		if not os.path.exists(filename):
			raise argparse.ArgumentTypeError("Cannot access the traits input file")
		return filename
		
	def is_file_of_nontraits_valid(filename):
		if not os.path.exists(filename):
			raise argparse.ArgumentTypeError("Cannot access the nontraits input file")
		return filename
		
	def is_kmer_valid(value_str):
		if value_str.isdigit():
			kmer = int(value_str)
			if kmer%2 == 1 and kmer >= 21 and kmer <= 127:
				return kmer
		raise argparse.ArgumentTypeError("Invalid Kmer value, it must be an odd integer between 21 and 127")
		
	def is_min_kmers_per_read(value_str):
		min_kmers = float(value_str)
		if min_kmers > 0 and min_kmers < 1:
				return min_kmers
		raise argparse.ArgumentTypeError("Invalid min percentage kmer coverage, must be between 0 and 1")
		
	def is_min_contig_len_valid(value_str):
		if value_str.isdigit():
			min_contig_len = int(value_str)
			if  min_contig_len >= 0 and min_contig_len <= 1000000:
				return min_contig_len
		raise argparse.ArgumentTypeError("Invalid minimum contig length value, try a resonable value like 600")
	
	def is_min_kmers_threshold_valid(value_str):
		if value_str.isdigit():
			min_kmers_threshold = int(value_str)
			if  min_kmers_threshold >= 0 and min_kmers_threshold <= 255:
				return min_kmers_threshold
		raise argparse.ArgumentTypeError("Invalid minimum kmers threshold, it must be between 0 and 255, but ideally between 20-30.")
		
	def is_threads_valid(value_str):
		if value_str.isdigit():
			threads = int(value_str)
			if  threads > 0 and threads <= 32:
				return threads
		raise argparse.ArgumentTypeError("Invalid number of threads, it must at least 1 and less than 32.")
		
	def is_max_kmers_threshold_valid(value_str):
		if value_str.isdigit():
			max_kmers_threshold = int(value_str)
			if  max_kmers_threshold >= 10 and max_kmers_threshold <= 255:
				return max_kmers_threshold
		raise argparse.ArgumentTypeError("Invalid maximum kmers threshold, it must be between 10 and 255, and greater than the minimum kmer value, but ideally greater than the coverage.")
			
	def is_min_spades_contig_coverage_valid(value_str):
		if value_str.isdigit():
			min_coverage = int(value_str)
			if  min_coverage >= 0 and min_coverage <= 1000000:
				return min_coverage
		raise argparse.ArgumentTypeError("Invalid minimum SPAdes contig coverage, try between 20 and 30.")
		
	def is_max_spades_contig_coverage_valid(value_str):
		if value_str.isdigit():
			max_coverage = int(value_str)
			if  max_coverage >= 0 and max_coverage <= 1000000:
				return max_coverage
		raise argparse.ArgumentTypeError("Invalid maximum SPAdes contig coverage, try between 1000.")

