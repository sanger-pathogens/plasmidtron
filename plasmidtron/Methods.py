import os
import logging
import subprocess
import re
import time
 
'''Produce a file with a paragraph of text and references which can easily be used in a publication'''
class Methods:
	def __init__(self, filename,trait_samples, nontrait_samples, min_kmers_threshold, min_contig_length, start_time, spades_exec ):
		self.logger = logging.getLogger(__name__)
		self.filename              = filename           
		self.trait_samples         = trait_samples      
		self.nontrait_samples      = nontrait_samples   
		self.min_kmers_threshold   = min_kmers_threshold
		self.min_contig_length     = min_contig_length  
		self.start_time            = start_time         
		self.spades_exec           = spades_exec
		
	def plasmidtron_version(self):
		current_version= 'X'
		try:
		    current_version = get_distribution('plasmidtron').version
		except:
		    current_version= 'X'
		return str(current_version)
			
	def decode(self,x):
	    try:
	        s = x.decode()
	    except:
	        return x
	    return s
		
	def kmc_version(self):
		#K-Mer Counter (KMC) ver. 2.3.0 (2015-08-21)
		regex = re.compile('ver. (\d+\.\d+\.\d+) ')
		cmd_output = subprocess.Popen('kmc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		cmd_output = self.decode(cmd_output[0]).split('\n')[:-1] + self.decode(cmd_output[1]).split('\n')[:-1]
		
		hits = regex.search(cmd_output[0])
		if hits:
			return str(hits.group(1))
		else:
			return 'X'
	
	# TODO remove duplication
	def spades_version(self):
		#SPAdes v3.9.0
		regex = re.compile('SPAdes v(\d+\.\d+\.\d+)')
		cmd_output = subprocess.Popen(self.spades_exec + ' -v', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		cmd_output = self.decode(cmd_output[0]).split('\n')[:-1] + self.decode(cmd_output[1]).split('\n')[:-1]

		hits = regex.search(cmd_output[0])
		if hits:
			return str(hits.group(1))
		else:
			return 'X'
			
	def running_time(self):
		return str(int(time.time()) - self.start_time)

	def create_file(self):
		self.logger.info("Creating file with methods")
		with open(self.filename, 'w') as output_fileh:
			output_fileh.write(self.methods_paragraph(len(self.trait_samples),len(self.nontrait_samples), self.plasmidtron_version(),self.kmc_version(), self.min_kmers_threshold, self.spades_version(), self.min_contig_length, self.running_time() ))
			output_fileh.write('\n\nReferences\n')
			output_fileh.write(self.references_paragraph())

	def methods_paragraph(self, num_trait_samples, num_non_trait_samples, plasmidtron_version,kmc_version, min_kmers_threshold, spades_version, min_contig_length, running_time ):
		methods_text =  'Method\n\n'
		methods_text += 'A set of '+ str(num_trait_samples) +' samples displayed a known phenotype/trait.\n'
		methods_text += 'To investigate this, these samples were sequenced using Illumina to produced paired ended reads in FASTQ format.\n'
		methods_text += 'These were provided to PlasmidTron (v'+ plasmidtron_version +')(Page et. al., 2017) along with '+ str(num_non_trait_samples)+' control samples which do not display the phenotype/trait.\n'
		methods_text += 'A k-mer analysis of the sets was performed using KMC (v'+ kmc_version +') (Deorowicz et. al., 2015) where a k-mer must be seen at least '+ str(min_kmers_threshold) +' times in a sample to be considered.\n'
		methods_text += 'The paired reads, which contained k-mers only found in the phenotype/trait set, were extracted and de novo assembled with SPAdes (v'+ spades_version +') (Bankevich et. al., 2012).\n'
		methods_text += 'Sequences in the resulting assemblies less than '+str(min_contig_length)+' bases long were excluded. The total running time was '+str(running_time)+' seconds.\n'
		return methods_text
		
	def references_paragraph(self):
		references_text =  'Anton Bankevich, Sergey Nurk, et. al. , and Pavel A. Pevzner. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19(5) (2012), 455-477. doi:10.1089/cmb.2012.0021\n'
		references_text += 'Andrew J. Page, Alexander Wailan,Yan Shao, Nicholas R. Thomson, Jacqueline A. Keane, "PlasmidTron: kmer based de novo assembly of genome sequences based on phenotypes", in preparation (2017).\n'
		references_text += 'Deorowicz, S., Kokot, M., Grabowski, Sz., Debudaj-Grabysz, A., KMC 2: Fast and resource-frugal k-mer counting, Bioinformatics, 2015; doi: 10.1093/bioinformatics/btv022.\n'
		return references_text
		