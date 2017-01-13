import unittest
import os
from plasmidtron.PlasmidTron import PlasmidTron

class Options:
	def __init__(self,output_directory, file_of_trait_fastqs, file_of_nontrait_fastqs, verbose, threads, kmer, min_kmers_threshold, spades_exec, min_contig_len):
		self.output_directory        = output_directory 
		self.file_of_trait_fastqs    = file_of_trait_fastqs
		self.file_of_nontrait_fastqs = file_of_nontrait_fastqs
		self.verbose                 = verbose
		self.threads                 = threads
		self.kmer                    = kmer
		self.min_kmers_threshold     = min_kmers_threshold
		self.spades_exec             = spades_exec
		self.min_contig_len          = min_contig_len

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plasmidtron')

class TestPlasmidTron(unittest.TestCase):
	
	def test_small_valid_chrom_plasmid(self):
		'''Given small FASTQS of simulated reads, with a chromosome in 1 and chromosome+plasmid in the other run the whole pipeline'''
		options = Options(os.path.join(data_dir,'out'), os.path.join(data_dir,'traits.csv'), os.path.join(data_dir,'nontraits.csv'),False, 1, 81, 20, 'spades-3.9.0.py', 100)
		
		plasmid_tron = PlasmidTron(options)
		plasmid_tron.run()
	
