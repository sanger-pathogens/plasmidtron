import os
import argparse

class InputTypes:
	
	def is_output_directory_valid(filename):
		if os.path.exists(filename):
			raise argparse.ArgumentTypeError("The output directory already exists")
		return filename
	
	def is_file_of_trait_fastqs_valid(filename):
		if not os.path.exists(filename):
			raise argparse.ArgumentTypeError("Cannot access the traits input file")
		return filename
		
	def is_file_of_nontrait_fastqs_valid(filename):
		if not os.path.exists(filename):
			raise argparse.ArgumentTypeError("Cannot access the nontraits input file")
		return filename
		
	def is_kmer_valid(value_str):
		if isinstance( value_str, int ):
			kmer = int(value_str)
			if kmer%2 == 1 and kmer >= 21 and kmer <= 127:
				return kmer
		raise argparse.ArgumentTypeError("Invalid Kmer value, it must be an odd integer between 21 and 127")
		
	def is_min_contig_len_valid(value_str):
		if isinstance( value_str, int ):
			min_contig_len = int(value_str)
			if  min_contig_len >= 0 and min_contig_len <= 1000000:
				return min_contig_len
		raise argparse.ArgumentTypeError("Invalid minimum contig length value, try a resonable value like 600")
	
	def is_min_kmers_threshold_valid(value_str):
		if isinstance( value_str, int ):
			min_kmers_threshold = int(value_str)
			if  min_kmers_threshold >= 0 and min_kmers_threshold <= 255:
				return min_kmers_threshold
		raise argparse.ArgumentTypeError("Invalid minimum kmers threshold,t must be between 0 and 255, but ideally between 20-30.")
		
	def is_threads_valid(value_str):
		if isinstance( value_str, int ):
			threads = int(value_str)
			if  threads > 0 and threads <= 512:
				return threads
		raise argparse.ArgumentTypeError("Invalid number of threads, it must at least 1 and less than the No. of CPUs")
