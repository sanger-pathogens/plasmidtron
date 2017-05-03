import unittest
import os
import shutil
from plasmidtron.PlasmidTron import PlasmidTron

class Options:
	def __init__(self,output_directory, file_of_traits, file_of_nontraits, verbose, threads, kmer, min_kmers_threshold,max_kmers_threshold, spades_exec, min_contig_len, action, min_spades_contig_coverage, keep_files,plot_filename):
		self.output_directory           = output_directory 
		self.file_of_traits       = file_of_traits
		self.file_of_nontraits    = file_of_nontraits
		self.verbose                    = verbose
		self.threads                    = threads
		self.kmer                       = kmer
		self.min_kmers_threshold        = min_kmers_threshold
		self.spades_exec                = spades_exec
		self.min_contig_len             = min_contig_len
		self.max_kmers_threshold        = max_kmers_threshold
		self.action                     = action
		self.min_spades_contig_coverage = min_spades_contig_coverage
		self.keep_files = keep_files
		self.plot_filename = plot_filename

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plasmidtron')

class TestPlasmidTron(unittest.TestCase):
	
	def test_small_valid_chrom_plasmid(self):
		'''Given small FASTQS of simulated reads, with a chromosome in 1 and chromosome+plasmid in the other run the whole pipeline'''
		if os.path.exists(os.path.join(data_dir,'out')):
			shutil.rmtree(os.path.join(data_dir,'out'))
		options = Options(os.path.join(data_dir,'out'), os.path.join(data_dir,'traits.csv'), os.path.join(data_dir,'nontraits.csv'),False, 1, 61, 20,200, 'spades.py', 100,'union', 1, False, 'plot.png')
		
		plasmid_tron = PlasmidTron(options)
		plasmid_tron.run()
		
		final_assembly = os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/filtered_scaffolds.fasta')
		
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/scaffolds.fasta')))
		self.assertTrue(os.path.isfile(final_assembly))
		'''The final assembly should be about 6k so leave some margin for variation in SPAdes'''
		self.assertTrue(os.path.getsize(final_assembly) > 5000)
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/plot.png')))
		shutil.rmtree(os.path.join(data_dir,'out'))
		
		
	def test_small_multithreaded(self):
		'''Run in parallel'''
		if os.path.exists(os.path.join(data_dir,'out')):
			shutil.rmtree(os.path.join(data_dir,'out'))
		options = Options(os.path.join(data_dir,'out'), os.path.join(data_dir,'traits.csv'), os.path.join(data_dir,'nontraits.csv'),False, 2, 61, 20,200, 'spades.py', 100,'union', 1, False, 'plot.png')
		
		plasmid_tron = PlasmidTron(options)
		plasmid_tron.run()
		
		final_assembly = os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/filtered_scaffolds.fasta')
		
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/scaffolds.fasta')))
		self.assertTrue(os.path.isfile(final_assembly))
		'''The final assembly should be about 6k so leave some margin for variation in SPAdes'''
		self.assertTrue(os.path.getsize(final_assembly) > 5000)
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/plot.png')))
		shutil.rmtree(os.path.join(data_dir,'out'))
		

	def test_non_trait_fasta(self):
		'''Given a nontrait file consisting only of a FASTA file, run the whole pipeline'''
		if os.path.exists(os.path.join(data_dir,'out')):
			shutil.rmtree(os.path.join(data_dir,'out'))
		options = Options(os.path.join(data_dir,'out'), os.path.join(data_dir,'traits.csv'), os.path.join(data_dir,'fasta_nontraits.csv'),False, 1, 61, 20,200, 'spades.py', 100,'union', 1, False, 'plot.png')
		
		plasmid_tron = PlasmidTron(options)
		plasmid_tron.run()
		
		final_assembly = os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/filtered_scaffolds.fasta')
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/spades_S_typhi_CT18_chromosome_pHCM2/scaffolds.fasta')))
		self.assertTrue(os.path.isfile(final_assembly))
		'''The final assembly should be about 6k so leave some margin for variation in SPAdes'''
		self.assertTrue(os.path.getsize(final_assembly) > 5000)
		self.assertTrue(os.path.isfile(os.path.join(data_dir,'out/plot.png')))
		shutil.rmtree(os.path.join(data_dir,'out'))